
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils_alt import  crossPCF, getAnnulusAreasAroundPoints, plotPCFWithBootstrappedConfidenceInterval, getPCFContributionsWithinGrid
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

        parser.add_argument(
                '-s', '--save_roi_data', 
                help = 'save individual ROI data',
                required = False,
                default = False
        )

        args = parser.parse_args()

        
        pathToData = args.pathToData
        clusteringToUse = args.clusteringToUse
        cluster_annotations = args.cluster_annotations
        
        pathToSaveFigures = args.output+"/"
        os.makedirs(pathToSaveFigures,exist_ok=True)

        clf= open(cluster_annotations).read().split("\n")
        clf= [x.strip().split("\t") for x in clf[1:]]
        mappings = {}
        ca = []
        cb = []
        for i in clf:
                if len(i) <2:
                        continue
                clname= "cl"+i[0]
                if len(i[0])==1:
                        clname= "cl0"+i[0]
                mappings[i[1]]=clname
                ca.append(i[1])
                cb.append(i[1])
        
        conf = json.loads(open(args.config).read())

        allClustersToCompare = []
        #only certain interactions specified
        if conf.get("clusters_to_run"):
            ca= conf["clusters_to_run"]["cell_type_A"]
            cb = conf["clusters_to_run"]["cell_type_B"]
        for a in ca:
            for b in cb:
                allClustersToCompare.append([mappings[a],mappings[b]])


        diseases= conf["conditions"]
        diseasesToAverage = list(diseases.keys())
        summary_data=[]
        for disease in diseasesToAverage:
                rois = diseases[disease]
                df_annotations, datasets = preprocessing(pathToData, cluster_annotations,rois)
                clusterNames = {df_annotations.ClusterNumber.iloc[v] : df_annotations.Annotation.iloc[v].strip() for v in range(len(df_annotations))}
                #calculate average number of cells
                cell_totals= {k:0  for k in clusterNames}
                for ds in datasets:
                    vs = dict(ds.df[clusteringToUse].value_counts())
                    for k in vs:
                        cell_totals[k]+=vs[k]
                cell_averages = {k:int(v/len(datasets)) for k,v in cell_totals.items()}
                for clustersToCompare in allClustersToCompare:
                        
                        allContributions = []
                        allPCFs = []
                        all_nRectangles = []
                        all_rectContributions = []
                        all_rectNs = []
                        
                        #summary line
                        cl1=clustersToCompare[0]
                        cl2= clustersToCompare[1]
                        line  = "{}\t{}\t{}\t{}\t{}".format(disease,clusterNames[cl1],clusterNames[cl2],cell_averages[cl1],cell_averages[cl2])

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
                                radii, g, contributions = crossPCF(distances_AtoB, areas_A, density_B, dr_mum, maxR_mum)
                                g = g.transpose()[0]
                                allContributions.append(contributions)
                                allPCFs.append(g)
                                
                                if args.save_roi_data:
                                    plt.figure(figsize=(12,9))
                                    plotPCFWithBootstrappedConfidenceInterval(plt.gca(), radii, g, contributions, p_A, ds.domainX, ds.domainY, includeZero=True)
                                    plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                                    os.makedirs(pathToSaveFigures + disease + '/Samples',exist_ok=True)
                                    plt.savefig(pathToSaveFigures + disease + '/Samples/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.png')
                                    plt.close()
                                
                                print("Pair correlation function completed.")
                                nRectangles, rectContributions, rectNs = getPCFContributionsWithinGrid(contributions, domainX, domainY, p_A)
                                all_nRectangles.append(nRectangles)
                                all_rectContributions.append(rectContributions)
                                all_rectNs.append(rectNs)
                                
                                # Now dump the data too
                               
                                os.makedirs(pathToSaveFigures + disease + '/Pickles/',exist_ok=True)
                                if args.save_roi_data:
                                    toSave = {'disease':disease,
                                    'name':ds.name,
                                    'PCF_radii':radii,
                                    'contributions':contributions
                                    }
                                    with open(pathToSaveFigures + disease + '/Pickles/' + ds.name + '_' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.p',"wb") as fid:
                                            pickle.dump(toSave,fid)
                        
                       
                        if len(allContributions)==0:
                                print ("No data for {} in {}".format(clustersToCompare,disease))
                                summary_data.append(line+"\tND\tND\tND\tND\tND\tND\tND\tND\tND\n")
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
                        stub = pathToSaveFigures + disease + '/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '_'
                        plt.savefig(stub+'averageWithCI95.png')
                        

                        plt.close()
                        to_save={
                            "radii":radii,
                            "PCF_mean":PCF_mean_fromContributions,
                            "PCF_min":PCF_min_fromContributions,
                            "PCF_max":PCF_max_fromContributions
                        }
                        with open(pathToSaveFigures + disease + '/Pickles/' + clusterNames[pair[0]] + '-to-' + clusterNames[pair[1]] + '.p',"wb") as fid:
                                        pickle.dump(to_save,fid)

                        #store line to add to summary file
                        me = to_save["PCF_mean"]
                        up = to_save["PCF_min"]
                        lo = to_save["PCF_max"]
                       
                        line += "\t{}\t{}\t{}".format( max(to_save["PCF_mean"]),to_save["PCF_mean"][1],to_save["PCF_mean"][2])
                        line += "\t{}\t{}\t{}\t{}".format(up[1],lo[1],up[2],lo[2])
                        r_min="NA"
                        r_mean="NA"
                        if me[0]>1:
                            for i,n in enumerate(to_save["radii"]):
                                if lo[i]<=1 and r_min =="NA":
                                    r_min=str(n)
                                    
                                if me[i]<=1 and r_mean=="NA":
                                    r_mean=str(n)
                                    break
                       
                        line+="\t{}\t{}\n".format(r_min,r_mean)
                        summary_data.append(line)
                        
        sf = open(pathToSaveFigures+"summary.txt","w")
        sf.write("state\tCell Type 1 1\tCell Type 2\tmean cell 1 number\tmean cell 2 number")
        sf.write("\tg(r)max\tgr10\tgr20")
        sf.write("\tgr10 PCF upper\tgr10 PCF lower\tgr20 PCF upper\tgr20 PCF lower")
        sf.write ("\tPCF min intersect\tPCF max intersect\n")
        for line in summary_data:
            sf.write(line)
        sf.close()            
                        


if __name__ == "__main__":
    main()   
    
















