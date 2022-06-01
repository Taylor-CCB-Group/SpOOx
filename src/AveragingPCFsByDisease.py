from matplotlib import path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from utils_alt import dataset,dataset_filterSampleID, crossPCF, getAnnulusAreasAroundPoints, plotPCFWithBootstrappedConfidenceInterval, getPCFContributionsWithinGrid
from scipy.spatial.distance import cdist
import os
import pickle
import argparse
import json
sns.set_style('white')
sns.set(font_scale=2)

from spatialstats import preprocessing

def main():     
        parser = argparse.ArgumentParser(
        )
        parser.add_argument(
                '-i', '--pathToData', 
                help = 'File path to data table (must be a tab-delimited file).',
                required = True
        ) 
        parser.add_argument(
                '-o', '--output', 
		help = 'Path to write all the outputs to.',
                required = True
        )
        parser.add_argument(
                '-c', '--cluster_annotations', 
                help = 'Input file with cluster annotations.',
                required = True
        )
        parser.add_argument(
                '-cl', '--clusteringToUse', 
                help = 'Name of column in the file specified by pathToData to use as the clustering label. Should correspond with the annotations in cluster_annotations',
                required = True
        )
        parser.add_argument(
                '-j', '--config', 
                help = 'json config',
                required = True
        )

        args = parser.parse_args()

        
        pathToData = args.pathToData
        clusteringToUse = args.clusteringToUse
        cluster_annotations = args.cluster_annotations
        
        pathToSaveFigures = args.output
        os.makedirs(pathToSaveFigures,exist_ok=True)

        clf= open(cluster_annotations).read().split("\n")
        clf= [x.strip().split("\t") for x in clf[1:]]
        mappings = {}
        for i in clf:
                if len(i) <2:
                        continue
                clname= "cl"+i[0]
                if len(i[0])==1:
                        clname= "cl0"+i[0]
                mappings[i[1]]=clname

        conf = json.loads(open(args.config).read())

        allClustersToCompare = []

        ca= conf["clusters_to_run"]["cell_type_A"]
        cb = conf["clusters_to_run"]["cell_type_B"]
        for a in ca:
                for b in cb:
                        if a==b:
                                continue
                        allClustersToCompare.append([mappings[a],mappings[b]])
        

        diseases= conf["disease_state"]
        diseasesToAverage = list(diseases.keys())
        for disease in diseasesToAverage:
                rois = diseases[disease]

                df_annotations, datasets = preprocessing(pathToData, cluster_annotations,rois)
                clusterNames = {df_annotations.ClusterNumber.iloc[v] : df_annotations.Annotation.iloc[v].strip() for v in range(len(df_annotations))}

                for clustersToCompare in allClustersToCompare:
                        
                        allContributions = []
                        allPCFs = []
                        all_nRectangles = []
                        all_rectContributions = []
                        all_rectNs = []
                        
                        for ds in datasets:
                                domainX = ds.domainX
                                domainY = ds.domainY
                                
                                dr_mum = 10
                                maxR_mum = 300      
                                
                                # First we pre-calculate the area around each point (within the domain)
                                areas = {}
                                allPoints = {}
                                has_data =True
                                for i, cluster in enumerate(clustersToCompare):
                                        print('Cluster',cluster)
                                
                                        p0 = ds.points[ds.df[clusteringToUse] == cluster]
                                        if np.shape(p0)[0]>0:
                                                p0_areas = getAnnulusAreasAroundPoints(p0, dr_mum, maxR_mum, ds.domainX, ds.domainY)
                                                allPoints[cluster] = p0
                                                areas[cluster] = p0_areas
                                        else:
                                                has_data=False
                                if has_data==False:
                                        print("No cells for {} {}".format(ds.name,clustersToCompare))
                                        continue
                                # Now we use that to find the PCFs
                                a = 0
                                b = 0
                                clusterA = clustersToCompare[0]
                                clusterB = clustersToCompare[1]
                                pair = [clusterA, clusterB]
                                print(pair)
                                
                                p_A = ds.points[ds.df[clusteringToUse] == clusterA]
                                p_B = ds.points[ds.df[clusteringToUse] == clusterB]
                                
                                density_B = np.shape(p_B)[0]/(ds.domainX*ds.domainY)
                                
                                areas_A = areas[clusterA]
                                areas_B = areas[clusterB]
                                
                                distances_AtoB = cdist(p_A, p_B, metric='euclidean')
                                radii, g, contributions = crossPCF(distances_AtoB, areas_A, areas_B, density_B, dr_mum, maxR_mum)
                                g = g.transpose()[0]
                                allContributions.append(contributions)
                                allPCFs.append(g)
                                
                                plt.figure(figsize=(12,9))
                                plotPCFWithBootstrappedConfidenceInterval(plt.gca(), radii, g, contributions, p_A, ds.domainX, ds.domainY, label=ds.indication, includeZero=True)
                                plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                                
                                print("Pair correlation function completed.")
                                
                                os.makedirs(pathToSaveFigures + disease + '/Samples',exist_ok=True)
                                plt.savefig(pathToSaveFigures + disease + '/Samples/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.png')
                                plt.close()
                                
                                nRectangles, rectContributions, rectNs = getPCFContributionsWithinGrid(contributions, domainX, domainY, p_A)
                                all_nRectangles.append(nRectangles)
                                all_rectContributions.append(rectContributions)
                                all_rectNs.append(rectNs)
                                
                                # Now dump the data too
                                toSave = {'disease':disease,
                                'name':ds.name,
                                'PCF_radii':radii,
                                'contributions':contributions
                                }
                                os.makedirs(pathToSaveFigures + disease + '/Pickles/',exist_ok=True)
                                with open(pathToSaveFigures + disease + '/Pickles/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.p',"wb") as fid:
                                        pickle.dump(toSave,fid)
                        
                       
                        if len(allContributions)==0:
                                print ("No data for {} in {}".format(clustersToCompare,disease))
                                continue
                        #%% Bootstrapping
                        # Each bootstrap sample, we select totalNrectangles quadrats and construct a PCF from them
                        totalNrectangles = np.sum(all_nRectangles)
                        # Since all rectangles are the same size, we can treat the whole lot as one big sample
                        rectContributions = np.concatenate(all_rectContributions, axis=0)
                        rectNs = np.concatenate(all_rectNs, axis=0)
                        
                        numBootstrapSims = 999
                        samplePCFs = np.zeros(shape=(numBootstrapSims, np.shape(rectContributions)[1]))
                        toSample = np.random.choice(totalNrectangles, size=(totalNrectangles,numBootstrapSims))

                        sample = np.sum(rectContributions[toSample,:],axis=0)
                        Ns = np.sum(rectNs[toSample],axis=0)
                    
                        samplePCFs = sample / Ns[:,np.newaxis]
                    
                        # Average according to mean over cells, rather than ROIs
                        allContributions_concat = np.concatenate(allContributions, axis=0)
                        PCF_mean_fromContributions = np.mean(allContributions_concat,axis=0)
                        PCF_max_fromContributions = 2*PCF_mean_fromContributions - np.percentile(samplePCFs, 97.5, axis=0)
                        PCF_min_fromContributions = 2*PCF_mean_fromContributions - np.percentile(samplePCFs, 2.5, axis=0)

                        plt.figure(figsize=(12,9))
                        plt.plot(radii,PCF_mean_fromContributions)
                        plt.gca().fill_between(radii, PCF_min_fromContributions, PCF_max_fromContributions, alpha=0.4)
                        plt.gca().axhline(1,c='k',linestyle=':')
                        plt.ylim([0,5])
                        #plt.gca().set_ylim(ymin=0)
                        plt.gca().set_xlabel('r ($\\mu$m)')
                        plt.gca().set_ylabel('$g(r)$')
                        #plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                        plt.savefig(pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_averageWithCI95.png')
                        plt.close()
                        


if __name__ == "__main__":
    main()   
    
















