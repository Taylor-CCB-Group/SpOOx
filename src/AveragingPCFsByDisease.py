from matplotlib import path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from utils_alt import dataset,dataset_filterSampleID, crossPCF, getAnnulusAreasAroundPoints, plotPCFWithBootstrappedConfidenceInterval, CalculateBootstrapAroundCSRForPValues
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

                # 1. Mat Neutrophils vs endothelial cells
                # cl05 vs cl37
                
                # 2. MAC2 vs Myofibroblast
                # cl06 vs cl33
                
                # 3. MAC1 vs MAC1 self interaction
                # cl17 vs cl17
                
                # 4. MAC2 vs MAC3
                # cl06 vs cl02
                
                # 5. Mat Neutrophils vs RAGEAlveolarE
                # cl05 vs cl43
                
                # 6. CCR2hi monocytes vs RAGEAlveolarE
                # cl07 vs cl43
                
                # 7. endothelial cell vs endothelial cell
                # cl37 vs cl37
                
                # 8. CCR2lo monocytes vs CD45 RO CD4 T cells 10μm
                # cl09 vs cl30
                
                # 9. CCR2 mid mono and Endothelial cells
                # cl18 vs cl37
        

                for clustersToCompare in allClustersToCompare:
                        
                        allContributions = []
                        allBootstraps = []
                        for ds in datasets:
                                domainX = ds.domainX
                                domainY = ds.domainY
                                
                                dr_mum = 20
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
                                
                                numBootstrapSamples = 1000
                                N_A = np.shape(p_A)[0]
                                N_B = np.shape(p_B)[0]
                                PCF_radii, bootstrapSamples = CalculateBootstrapAroundCSRForPValues(N_A, N_B, domainX, domainY, dr_mum, maxR_mum, numBootstrapSamples)
                                
                                plt.figure(figsize=(12,9))
                                plotPCFWithBootstrappedConfidenceInterval(plt.gca(), radii, g, contributions, p_A, ds.domainX, ds.domainY, label=ds.indication, includeZero=True)
                                plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                                
                                # Plot halos
                                plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),5,axis=1),color=[1,0,0],linestyle=':')
                                plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),95,axis=1),color=[1,0,0],linestyle=':')
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),2.5,axis=1),color=[1,0,1],linestyle=':')
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),97.5,axis=1),color=[1,0,1],linestyle=':')
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),1,axis=1),color=[1,0.6,0.6])
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),99,axis=1),color=[1,0.6,0.6])
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),0.5,axis=1),color=[1,0.8,0.8])
                                # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),99.5,axis=1),color=[1,0.8,0.8])
                                # plt.gca().axhline(1,color=[0,0,0],linestyle=':')
                                allBootstraps.append(bootstrapSamples)
                                print("Pair correlation function completed.")
                                
                                os.makedirs(pathToSaveFigures + disease + '/Samples',exist_ok=True)
                                plt.savefig(pathToSaveFigures + disease + '/Samples/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.png')
                                plt.close()
                                
                                # Now dump the data too
                                toSave = {'disease':disease,
                                'name':ds.name,
                                'bootstrapSamples':bootstrapSamples,
                                'PCF_radii':PCF_radii,
                                'contributions':contributions
                                }
                                os.makedirs(pathToSaveFigures + disease + '/Pickles/',exist_ok=True)
                                with open(pathToSaveFigures + disease + '/Pickles/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.p',"wb") as fid:
                                        pickle.dump(toSave,fid)
                        
                       
                        if len(allContributions)==0:
                                print ("No data for {} in {}".format(clustersToCompare,disease))
                                continue
                        stackedContributions = []
                        pcfMeans = []
                        ns = []
                        
                        for v in range(len(allContributions)):
                                stackedContributions.extend(allContributions[v])
                                pcfMeans.append(np.mean(allContributions[v],axis=0))
                                ns.append(np.shape(allContributions[v])[0])
                        
                        # Step 1: Equally weighted ROIs
                        # PCF - just average the PCFs
                        # Bootstrap: average bootstrap n for all ROIs to get numBootstrapSamples things, then take p-value
                        
                        plt.figure(figsize=(12,9))
                        plt.plot(radii,np.mean(pcfMeans,axis=0),label='Avg (ROIs equally weighted)')
                        bs = np.sum(allBootstraps,axis=0)/np.shape(allBootstraps)[0]
                        # Plot halos
                        plt.plot(radii,np.percentile(bs.transpose(),5,axis=1),color=[1,0,0],linestyle=':')
                        plt.plot(radii,np.percentile(bs.transpose(),95,axis=1),color=[1,0,0],linestyle=':')
                        plt.gca().set_xlabel('r ($\\mu$m)')
                        plt.gca().set_ylabel('g(r)')
                        plt.gca().axhline(y=1, c=[0.5, 0.5, 0.5], linestyle=':')
                        plt.title('Weight ROIs equally')
                        plt.savefig(pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_averageByROIs.png')
                        # plt.savefig(pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_averageByROIs.eps')
                        plt.close()
                        
                        
                        
                        # Step 2: Equally weighted points
                        # PCF - weighted average according to number of points
                        # Bootstrap - same thing
                        plt.figure(figsize=(12,9))
                        weightedAvg = np.sum([ns[v]*pcfMeans[v] for v in range(len(ns))],axis=0)/sum(ns)
                        bs = np.sum([ns[v]*allBootstraps[v] for v in range(len(ns))],axis=0)/sum(ns)
                        plt.plot(radii,weightedAvg,label='Avg (cells equally weighted)')
                        # Plot halos
                        plt.plot(radii,np.percentile(bs.transpose(),5,axis=1),color=[1,0,0],linestyle=':')
                        plt.plot(radii,np.percentile(bs.transpose(),95,axis=1),color=[1,0,0],linestyle=':')
                        plt.gca().set_xlabel('r ($\\mu$m)')
                        plt.gca().set_ylabel('g(r)')
                        plt.gca().axhline(y=1, c=[0.5, 0.5, 0.5], linestyle=':')
                        plt.title('Weight cells equally')
                        plt.savefig(pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_averageByCells.png')
                        # plt.savefig(pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_averageByROIs.eps')
                        plt.close()

if __name__ == "__main__":
    main()   
    




















