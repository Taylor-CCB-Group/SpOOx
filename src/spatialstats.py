#!/usr/bin/env python
# coding: utf-8

DESC = ("Apply spatial statistics pipeline with command line inputs.")

# Imports
import os
import os.path
import sys
import argparse
from matplotlib import path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils_alt import dataset,dataset_filterSampleID, crossPCF, getAnnulusAreasAroundPoints, plotPCFWithBootstrappedConfidenceInterval, changeSomeElements
import skimage
from pathlib import Path
from scipy.spatial.distance import cdist
import pickle
from statsmodels.stats.multitest import fdrcorrection
import skimage.io

sns.set_style('white')
sns.set(font_scale=2)


#name of current program
prog = os.path.basename(__file__)
description = DESC
rois=[]

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
                '-r', '--rois',
                nargs = '+',
                help = 'Name of roi in the file, must be in the sample_id column',
                required = False
        ) 
        parser.add_argument(
                '-f', '--functions',
                nargs = '+', 
                help = 'Input the functions that are to be carried out (celllocationmap, contourplots, quadratcounts, quadratcelldistributions, paircorrelationfunction, morueta-holme). If no argument is given, all functions will be carried out.',
                required = False
        ) 
        parser.add_argument(
                '-v', '--verbose', 
                action = 'store_true', 
                required = False
        )

        parser.add_argument(
                '-q', '--quadrat_size', 
                default=100,
                required = False,
                type=int
        )

   





        args = parser.parse_args()

        if (os.path.exists(args.cluster_annotations) == False):
                print("Annotation file does not exist.")
                sys.exit()



        pathToData = args.pathToData#'/project/covidhyperion/shared/data/panel2/tree/HEALTHY/SAMPLE_1/ROI_1/clustering/cellData.tab'#
        pathToWriteOutput = args.output#'/Filers/home/j/jbull/Temp1/'#
        cluster_annotations = args.cluster_annotations#'/project/covidhyperion/shared/data/panel2//config//exampleclusterannotation.tab'#
        clusteringToUse = args.clusteringToUse
        rois = args.rois

        # pathToData = '/t1-data/project/covidhyperion/shared/data/panel2/tree/clustering/COVID_HC/COVID_HC_cellData_Harmonyclustered_k30.txt'#'/t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/SAMPLE_14/ROI_2/clustering/COVID_SAMPLE_14_ROI_2_cellData_clustered.txt'
        # cluster_annotations = '/t1-data/project/covidhyperion/shared/data/panel2/config/clusterannotation.tab'
        # clusteringToUse = 'harmony_phenograph_exprs'#'phenograph_cluster_scaledtrim_k30'
        # pathToWriteOutput = '/Filers/home/j/jbull/Temp_2222_2/'
         
        
        functions = args.functions
     
        if bool(functions) == False:
                functions = ['celllocationmap', 'contourplots', 'quadratcounts', 'quadratcelldistributions','paircorrelationfunction','morueta-holme','networkstatistics','localclusteringheatmaps']
                print("Running all functions:", functions)
        
        print("Checking what samples and rois to analyse...")

        df_annotations, datasets = preprocessing(pathToData, cluster_annotations,rois)
        #todo should not need to add a / 
        
        for ds in datasets:
            ds.pathToWriteOutput = pathToWriteOutput+"/quadratMethods/"
            if not os.path.exists(ds.pathToWriteOutput):
                os.mkdir(ds.pathToWriteOutput)
            #clusteringToUse = 'phenoGraph_cluster' ###ADD to cmd line?
            if len(df_annotations.ClusterNumber) < 21:
                colors = [plt.cm.tab20(v) for v in range(len(df_annotations.ClusterNumber))]
            else:
                # We have more than 20 colors, so this is going to be annoying. Generate a new colormap
                nNeeded = len(df_annotations.ClusterNumber) - 20
                colors = [plt.cm.tab20(v) for v in range(20)]
                newColors = [plt.cm.Set3(v) for v in range(nNeeded)]
                colors.extend(newColors)
            
            clusterNames = {df_annotations.ClusterNumber.iloc[v] : df_annotations.Annotation.iloc[v].strip() for v in range(len(df_annotations))}
    
            #make an output directory
            Path(pathToWriteOutput).mkdir(parents=True, exist_ok=True)
            print("Output file directory created.")
    
            #Calling the functions.
            functions = [x.casefold() for x in functions]
            if 'quadratcounts' in functions or 'quadratcelldistributions' in functions or 'morueta-holme' in functions:
                    r, VMRs, counts = quadratMethods(ds, df_annotations, clusteringToUse,args.quadrat_size)
    
            for i in functions:
                    #get directory according to function and create the dirctory
                    #if one does not exist
                    ds.pathToWriteOutput=pathToWriteOutput+"/"+i+"/"
                    if not os.path.exists(ds.pathToWriteOutput):
                            os.mkdir(ds.pathToWriteOutput)
                    if i == 'celllocationmap':
                            cellLocationMap(ds, df_annotations, clusteringToUse, clusterNames, colors)
                    if i == 'contourplots':
                            contourPlots(ds, df_annotations, clusteringToUse, clusterNames)
                    if i == 'quadratcounts':
                            quadratCounts(ds, df_annotations, colors, clusterNames, r, VMRs)
                    if i == 'quadratcelldistributions':
                            quadratCellDistributions(ds, df_annotations, colors, clusterNames, counts)
                    if i == 'morueta-holme':
                            moruetaHolmeAssociationMatrix(ds, df_annotations, colors, clusterNames, counts)
                    if i == 'paircorrelationfunction':
                            pairCorrelationFunction(ds, df_annotations, clusteringToUse, clusterNames)
                    if i == 'networkstatistics':
                            # Infer path to label matrix from file structure - this is easily breakable (and hopefully easily fixable too)
                            #pathToLabelMatrix = '/project/covidhyperion/shared/data/panel2/tree/'+ ds.df.condition + '/'+ds.df.sample_id+'/'+ds.df.ROI+'/deepcell/COVIDPANEL2_'+ds.df.sample_id+'_'+ds.df.ROI+'.tif'
                            sub_path = str(ds.df.sample_id)
                            sub_path = sub_path.split()[1]
                            print("split path: ", sub_path)
                            sub_path = sub_path.replace('_SAMPLE', '/SAMPLE')
                            sub_path = sub_path.replace('_ROI', '/ROI')
                            print("sub_path" + sub_path)
                            pathToLabelMatrix = '/project/covidhyperion/shared/data/panel2/tree/' + sub_path + '/deepcell/deepcell.tif'
                            labels = skimage.io.imread(pathToLabelMatrix)
                            networkStatistics(ds, df_annotations, clusteringToUse, clusterNames, labels, colors)
                    if i == 'localclusteringheatmaps':
                            localClusteringHeatmaps(ds, df_annotations, clusteringToUse, clusterNames)


def preprocessing(pathToData, cluster_annotations,rois):
        # Filter dataframe to ensure spatial stats are only applied to the same sample   
        df = pd.read_csv(pathToData,delimiter='\t')        
        uniqueSamples = np.unique(df['sample_id'])
        print("Unique samples: ", uniqueSamples)

        # check to see if rois specified in cmd line are specified in the dataframe (sample_id)
	# if correct then run only those rois
        if (bool(rois) == True):
                print("ROIs selected =",rois)
                check =  all(item in uniqueSamples for item in rois)
                if (check == True):           
                        uniqueSamples = rois
                else:                                           
                        print("Some ROIs do not match samples in sample_id column. Please check and try again.")
                        sys.exit()
	
        datasets = []   

        print('Analysing',uniqueSamples)
        
        #if len(uniqueSamples) == 1:
        #    ds = dataset(pathToData)
        #    datasets.append(ds)
        if len(uniqueSamples) >= 1:
            # Split it up, return multiple datasets
            for sample_id in uniqueSamples:
                print('Preparing',sample_id)
                ds = dataset_filterSampleID(pathToData,sample_id)
                datasets.append(ds)
        else:
            print('No samples found')
            assert(1==2)
            #todo Throw a proper error dealing with this
        
        
        #PART 1 - Preprocess main data table & make path to write output.
        #ds = dataset(pathToData)
        print("Main data tables processed.")
        
        #PART 2 - Preprocess annotation table.
        df_annotations = pd.read_csv(cluster_annotations, delimiter='\t')
        # Put into annoying format
        df_annotations.ClusterNumber = ['cl' + f"{v:02d}" for v in  df_annotations.ClusterNumber]        # Find the correct ROI.
        #df_annotations = df_clusterAnnotations[(df_clusterAnnotations.Indication.str.casefold()==ds.indication.casefold()) & (df_clusterAnnotations.Sample.str.casefold()==ds.sample.casefold()) & (df_clusterAnnotations.ROI.str.casefold()==ds.roi.casefold())]
        print("Annotation table processed.")
        
        return df_annotations, datasets
        

def cellLocationMap(ds, df_annotations, clusteringToUse, clusterNames, colors):
        #PART 4 - Generate full map of cell locations
        plt.figure(figsize=(12,12))
        plt.title(ds.name)
        for i, cl in enumerate(df_annotations.ClusterNumber):
                p = ds.points[ds.df[clusteringToUse] == cl]
                plt.scatter(p[:,0],p[:,1],color=colors[i],s=10,label=clusterNames[cl] + ' (cluster ' + str(cl) + ')')
        plt.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)
        plt.axis('square')
        plt.xlim([0,ds.domainX])
        plt.ylim([0,ds.domainY])
        plt.savefig(ds.pathToWriteOutput + ds.name + '__Overview.png',bbox_inches='tight')
        plt.close()
        print("Cell location map completed.")


def contourPlots(ds, df_annotations, clusteringToUse, clusterNames):
        #PART 5
        for i, cl in enumerate(df_annotations.ClusterNumber):
                plt.figure(figsize=(12, 9))
                p = ds.points[ds.df[clusteringToUse] == cl]
                if np.shape(p)[0]>20:
                    plotting_points = pd.DataFrame(np.array(p), columns=['x', 'y'])
                    ax = sns.kdeplot(data=plotting_points,x='x',y='y', shade=True, palette="PuRd", cut=0, levels=15)#200
                    ax.axis('square')
                    plt.title(clusterNames[cl] + ' (cluster ' + str(cl) + ')')
                    plt.xlim([0, ds.domainX])
                    plt.ylim([0, ds.domainY])
                    plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[cl] + ' (cluster ' + str(cl) + ')' + '_KDEplot.png',bbox_inches='tight')
                    plt.close()
        print("Contour plots completed.")


def quadratMethods(ds, df_annotations, clusteringToUse, quadratEdgeLength=100):
        #PART 6 - Time to implement some quadrat methods
        
        x = np.arange(0,ds.domainX,quadratEdgeLength)
        y = np.arange(0,ds.domainY,quadratEdgeLength)
        x = np.append(x,ds.domainX)
        y = np.append(y,ds.domainY)

        xv, yv = np.meshgrid(x,y)
        nx = len(x)
        ny = len(y)
        nQuadrats = (nx-1)*(ny-1)

        quadrats = []
        for i in range(nx-1):
                for j in range(ny-1):
                        vertices = [(x[i],y[j]), (x[i],y[j+1]), (x[i+1],y[j+1]), (x[i+1],y[j])]
                        quadrats.append(vertices)
        assert(nQuadrats == len(quadrats))

        nSpecies = len(df_annotations.ClusterNumber)
        counts = np.zeros(shape=(nQuadrats,nSpecies))
        # totals = np.zeros(shape=(nSpecies))
        
        i = 0
        celllabels = ds.df[clusteringToUse]
        for xi in range(nx-1):
            for yi in range(ny-1):
                px = (ds.points[:,0] >= x[xi]) & (ds.points[:,0] < x[xi+1])
                py = (ds.points[:,1] >= y[yi]) & (ds.points[:,1] < y[yi+1])
                mask = px & py
                labs = celllabels[mask]
                for ind, key in enumerate(df_annotations.ClusterNumber):          
                    counts[i,ind] = sum(labs == key)
                i = i + 1

        VMRs = []
        for ind, cl in enumerate(df_annotations.ClusterNumber):
                mean = np.nanmean(counts[:,ind])
                variance = np.nanstd(counts[:,ind])**2
                VMR = variance / mean
                VMRs.append(VMR)            
            
        r = np.corrcoef(counts.transpose())
        
        saveFile = ds.pathToWriteOutput + ds.name + '_quadratCounts_data.p'
        toSave = {'PickleFileGeneratedBy':'spatialstats.py',
                  'quadratEdgeLength':quadratEdgeLength,
                  'counts':counts,
                  'r':r,
                  'VMRs':VMRs
                  }
                
        with open(saveFile,"wb") as fid:
            pickle.dump(toSave,fid)


        print("Quadrat methods completed.")

        return r, VMRs, counts


def quadratCounts(ds, df_annotations, colors, clusterNames, r, VMRs):
        #PART 7
        plt.figure(figsize=(12,9))
        g = sns.heatmap(r, cmap='RdBu_r', vmin=-1, vmax=1, annot=False, xticklabels=clusterNames.values(), yticklabels=clusterNames.values())
        plt.title('Quadrat counts: correlation coefficient')
        plt.savefig(ds.pathToWriteOutput + ds.name + '__QuadratCountsCorrelationCoefficient.png',bbox_inches='tight')

        plt.figure(figsize=(12,12))
        plt.scatter(VMRs,df_annotations.ClusterNumber,c=colors,s=400)
        # for i, cl in enumerate(clustersToCompare):
        #     plt.plot([cl,cl],[lowVMR99[i],highVMR99[i]],color=colors[i])
        plt.ylabel('Cluster ID')
        plt.xlabel('VMR')
        plt.gca().set_yticks(list(clusterNames.keys()))
        plt.gca().set_yticklabels(clusterNames.values())
        plt.plot([1,1],[np.sort(list(clusterNames.keys()))[0],np.sort(list(clusterNames.keys()))[-1]],'k:')
        plt.savefig(ds.pathToWriteOutput + ds.name + '__VarianceToMeanRatio.png',bbox_inches='tight')
        
        print("Quadrat counts completed.")


def quadratCellDistributions(ds, df_annotations, colors, clusterNames, counts):
        #PART 8 - Visualise distributions of cells in each quadrat
        import scipy.stats as ss

        sns.set(font_scale=1)

        for i, cl in enumerate(df_annotations.ClusterNumber):
                plt.figure(figsize=(6,6))
                b = np.arange(0,30,1)

                mu = sum(counts[:,i])/len(counts[:,i])
                #print(sum(counts[:,i]), mu, len(counts[:,i]))
                #print(ss.poisson.pmf(b, mu))
                #print(sum(counts[:,i])*ss.poisson.pmf(b, mu))
                #plt.figure()
                #plt.plot(b,ss.poisson.pmf(b, mu))
                f_exp = ss.poisson.pmf(b, mu)

                plt.plot(b, f_exp, 'bo', ms=2, label='poisson pmf')
                plt.hist(counts[:,i],color=colors[i], bins=b, density=True)
                plt.title(clusterNames[cl])

                plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[cl] + ' (cluster ' + str(cl) + ')' + '_HistogramOfQuadratCounts.png',bbox_inches='tight')
                plt.close()

        print("Quadrat cell distributions completed.")
        
def moruetaHolmeAssociationMatrix(ds, df_annotations, colors, clusterNames, counts):
        print('Calculating MH Association matrix for',ds.df.sample_id.iloc[0])
        # Calculate cell association matrix following the method from Morueta-Holme et al 2016 DOI: 10.1111/ecog.01892
        O = counts.transpose()
        clusterNames_MH = list(clusterNames.values())
       
        # If any species is not represented in this sample, O will be singular and hence non-invertible
        # So we get rid of any cell types that aren't present in the sample
        if len(np.where(np.sum(counts,axis=0) == 0)) > 0:
            toExclude = list(np.where(np.sum(counts,axis=0) == 0)[0])
            O = np.delete(O,toExclude,0)
            clusterNames_MH = np.delete(clusterNames_MH,toExclude,0)
        
        nSpecies = len(clusterNames_MH)
        # STEP ONE:
        # Calculate the observed species x site community matrix, O (nSpecies x nQuadrats)
        
        assert(np.shape(O)[0] == nSpecies)
        
        # STEP TWO:
        # Calculate the expected abundance matrix, N (nSpecies x nQuadrats)
        # Randomise O to create N, keeping row sums (total species) and column sums (site abundances) fixed
        # First, shuffle O 10000 times to avoid transient effects
        matrix = O.copy()
        i = 0
        while i < 10000:
            matrix, success = changeSomeElements(matrix)
            if success:
                i = i + 1
        
        # Now we get 1000 matrices N
        i = 0
        Ns = []
        print('Shuffling matrix')
        singularCount = 0
        while i < 1000:
            # Do 500 shuffles in between each N
            if i % 25 == 0:
                print('\t\t', i, 'of',1000)
            j = 0
            while j < 500:
                matrix, success = changeSomeElements(matrix)
                if success:
                    j = j + 1
            if np.linalg.det(np.cov(matrix)) != 0:
                # Only add matrices with non-zero determinant in covariance matrix
                out = matrix.copy()
                Ns.append(out)
                i = i + 1
            else:
                singularCount = singularCount + 1
                print('Found singular covariance matrix, reshuffling (' + str(singularCount) + ' occurrences so far)')
        Ns = np.asarray(Ns)
        
        # OK, we now have 1000 matrices with randomly shuffled values such that the row and column sums are the same!
        # That was a faff.
        
        # STEP THREE:
        # Calculate the covariance matrix Sigma for each of O and N
        # Use Sigma to calculate the partial correlation matrix C(O) and C(N) (nSpecies x nSpecies)
        # C_ij (M) = -Sigma^-1_ij (M) / sqrt(  Sigma^-1_ii (M) . Sigma^-1_jj (M)   )
        # Repeat for MANY resamples of N
        
        def getCFromSigma_inv(Sigma_inv):
            C = np.zeros(np.shape(Sigma_inv))
            for i in range(len(Sigma_inv)):
                for j in range(len(Sigma_inv)):
                    C[i,j] = -Sigma_inv[i,j] / np.sqrt(Sigma_inv[i,i] * Sigma_inv[j,j])
            return C
        
        Sigma_O = np.cov(O)
        Sigma_O_inv = np.linalg.inv(Sigma_O)
        C_O = getCFromSigma_inv(Sigma_O_inv)
        # Fix the main diagonal for plotting
        C_O_plot = C_O + 2*np.identity(nSpecies)
        plt.figure(figsize=(12,12))
        plt.imshow(C_O_plot,cmap='RdBu_r',vmin=-1,vmax=1)
        plt.colorbar()
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH)
        plt.savefig(ds.pathToWriteOutput + ds.name + '_PartialCorrelation.png',bbox_inches='tight')
        plt.close()

        C_ns = []
        for i in range(np.shape(Ns)[0]):
            Sigma = np.cov(Ns[i])
            Sigma_inv = np.linalg.inv(Sigma)
            C_n = getCFromSigma_inv(Sigma_inv)
            C_ns.append(C_n)
        C_ns = np.asarray(C_ns)      
        
        # STEP FIVE:
        # Calculate the standard effect size by rescaling the mean and SD of the null distributions
        SES = np.zeros(shape=(nSpecies,nSpecies))
        for i in range(nSpecies):
            for j in range(nSpecies):
                SES[i,j] = (C_O[i,j] - np.mean(C_ns[:,i,j])) / np.std(C_ns[:,i,j])
        
        plt.figure(figsize=(12,12))
        lim = np.ceil(max(np.abs(np.nanmin(SES)),np.nanmax(SES)))
        plt.imshow(SES,cmap='RdBu_r',vmin=-lim,vmax=lim)
        plt.colorbar()
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH)
        plt.savefig(ds.pathToWriteOutput + ds.name + '_StandardEffectSize.png',bbox_inches='tight')
        plt.close()

        # STEP SIX: 
        # Calculate a 2-tailed p-value for each pair of species
        p = np.zeros(shape=(nSpecies,nSpecies))
        for i in range(nSpecies):
            for j in range(nSpecies):
                p[i,j] = sum(np.abs(C_O[i,j]) < np.abs(C_ns[:,i,j]))/np.shape(C_ns)[0]
        
        # Now do Benjamini-Hochberg correction
        
        alpha = 0.05
        
        rejected, p_star = fdrcorrection(p.reshape((1,-1))[0], alpha, method='indep')
        
        rejected = rejected.reshape(nSpecies,nSpecies)
        p_star = p_star.reshape(nSpecies,nSpecies)
        
        A = np.zeros(shape=(nSpecies,nSpecies))
        for i in range(nSpecies):
            for j in range(nSpecies):
                if p_star[i,j] < alpha:
                    A[i,j] = SES[i,j]
                else:
                    A[i,j] = np.nan
                    
        plt.figure(figsize=(12,12))
        plt.imshow(A,cmap='RdBu_r',vmin=-15,vmax=15)
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames_MH)
        plt.colorbar()
        plt.savefig(ds.pathToWriteOutput + ds.name + '_AdjacencyMatrix.png',bbox_inches='tight')
        
        plt.close('all')
        
        saveFile = ds.pathToWriteOutput + ds.name + '_MoruetaHolme_data.p'
        toSave = {'PickleFileGeneratedBy':'spatialstats.py',
                  'A':A,
                  'C_O':C_O,
                  'C_ns':C_ns,
                  'counts':counts,
                  'clusterNames_MH':clusterNames_MH,
                  'Ns':Ns,
                  'nSpecies':nSpecies,
                  #'nx':nx,
                  #'ny':ny,
                  'O':O,
                  #'outfolder':outfolder,
                  'p':p,
                  'p_star':p_star,
                  #'pathToData':pathToData,
                  #'quadratEdgeLength':quadratEdgeLength,
                  #'r':r,
                  'rejected':rejected,
                  'ds.name':ds.name,
                  'Sigma_O':Sigma_O,
                  'Sigma_O_inv':Sigma_O_inv,
                  #'x':x,
                  #'y':y,
                  #'xlim':xlim,
                  #'ylim':ylim,
                  'ds':ds, 
                  'df_annotations':df_annotations, 
                  'colors':colors, 
                  'clusterNames':clusterNames}
                
        with open(saveFile,"wb") as fid:
            pickle.dump(toSave,fid)


        print("Quadrat neighbourhood analysis (Morueta-Holme) completed.")

def pairCorrelationFunction(ds, df_annotations, clusteringToUse, clusterNames):
        #PART 9 - Cell-cell methods
        # Compute the crossPCF between each pairwise combination of clusters
        maxR_mum = 300
        dr_mum = 10

        # First we pre-calculate the area around each point (within the domain)
        areas = {}
        allPoints = {}
        for i, cluster in enumerate(df_annotations.ClusterNumber):
                print('Cluster',cluster)

                p0 = ds.points[ds.df[clusteringToUse] == cluster]
                if np.shape(p0)[0]>0:
                    p0_areas = getAnnulusAreasAroundPoints(p0, dr_mum, maxR_mum, ds.domainX, ds.domainY)
                    allPoints[cluster] = p0
                    areas[cluster] = p0_areas
      
        # Now we use that to find the PCFs
        gs = np.zeros(shape=(len(df_annotations.ClusterNumber), len(df_annotations.ClusterNumber), len(np.arange(0,maxR_mum,dr_mum))))
        for a, clusterA in enumerate(df_annotations.ClusterNumber):
                for b, clusterB in enumerate(df_annotations.ClusterNumber):
                        pair = [clusterA, clusterB]
                        print(pair)
                        
                        p_A = ds.points[ds.df[clusteringToUse] == clusterA]#allPoints[clusterA]# ds.points[ds.df[clusteringToUse] == int(pair[0])]    
                        p_B = ds.points[ds.df[clusteringToUse] == clusterB]#allPoints[clusterB]#ds.points[ds.df[clusteringToUse] == int(pair[1])]
                        if np.shape(p_A)[0]>0 and np.shape(p_B)[0]>0:
                            density_B = np.shape(p_B)[0]/(ds.domainX*ds.domainY)
                            
                            areas_A = areas[clusterA]
                            areas_B = areas[clusterB]
                            
                            distances_AtoB = cdist(p_A, p_B, metric='euclidean')
                            radii, g, contributions = crossPCF(distances_AtoB, areas_A, areas_B, density_B, dr_mum, maxR_mum)
                            gs[a, b, :] = g.transpose()[0]
                    
                            plt.figure(figsize=(12,9))
                            plotPCFWithBootstrappedConfidenceInterval(plt.gca(), radii, g, contributions, p_A, ds.domainX, ds.domainY, label=ds.indication, includeZero=True)
                            plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                            plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]] + '_PCF.png',bbox_inches='tight')
                            plt.close()               
        print("Pair correlation function completed.")   
    

        #PART 10 - Make a heatmap
        #plt.imshow(gs[:,:,1],cmap='RdBu_r',vmin=0,vmax=2)
        for radius in [1,2]:
                fig = plt.figure(figsize=(11, 9))
                labs = clusterNames.values()
                g = sns.heatmap(gs[:, :, radius], cmap='RdBu_r', vmin=0, vmax=2, annot=False, xticklabels=labs, yticklabels=labs, cbar_kws={'label': '$g(r=' + str(radii[radius]) + '\mu m)$'})
                # plt.title(indications[i] + ' - ' + str(radii[radius]) + '$\mu$m')
                plt.title(ds.name + ' - ' + str(radii[radius]) + '$\mu$m')
                # fig.autofmt_xdate()
                plt.savefig(ds.pathToWriteOutput + ds.name + '_PCF-Heatmap' + str(radii[radius]) + '-micrometers.png',bbox_inches='tight')
                # Make graph image
                # Save everything
                plt.close()
           
        saveFile = ds.pathToWriteOutput + ds.name + '_PCF_data.p'
        toSave = {'PickleFileGeneratedBy':'spatialstats.py',
                  'maxR_mum':maxR_mum,
                  'dr_mum':dr_mum,
                  'gs':gs,
                  'radii':radii}
                
        with open(saveFile,"wb") as fid:
            pickle.dump(toSave,fid)
        print("PCF heatmap completed.")


def networkStatistics(ds, df_annotations, clusteringToUse, clusterNames, labels, colors):
        
        label2CellType = {}
        annotationDict = {df_annotations.iloc[v]['ClusterNumber']:df_annotations.iloc[v]['Annotation'] for v in range(len(df_annotations))}
        for i in range(len(ds.df)):
            label = ds.df.iloc[i].label
            cellType = ds.df.iloc[i][clusteringToUse] 
            label2CellType[label] = cellType
            
        import scipy.ndimage as ndimage
        # Dilate the label matrix https://stackoverflow.com/questions/12747319/scipy-label-dilation
        dilationValue = 5
        dilation = ndimage.maximum_filter(labels,dilationValue)
        dilation[labels != 0] = labels[labels != 0]

        # Get adjaceny matrix https://stackoverflow.com/questions/26486898/matrix-of-labels-to-adjacency-matrix
        G = np.zeros([dilation.max() + 1]*2)
        
        # left-right pairs
        G[dilation[:, :-1], dilation[:, 1:]] = 1
        # right-left pairs
        G[dilation[:, 1:], dilation[:, :-1]] = 1
        # top-bottom pairs
        G[dilation[:-1, :], dilation[1:, :]] = 1
        # bottom-top pairs
        G[dilation[1:, :], dilation[:-1, :]] = 1
        
        nx, ny = np.shape(G)
        
        # Remove these asserts as they only hold when no labels have been removed in filtering steps earlier
        #assert(np.max(labels) == np.max(ds.df.label))
        #assert(len(G) == np.max(ds.df.label)+1)
        
        # Quantify numbers of connections between different labels
        
        celltypes = np.asarray(ds.df[clusteringToUse])
        outputs = {}
        for ct_i, celltype in enumerate(annotationDict):
            posInds = np.where(ds.df[clusteringToUse] == celltype)[0]
            print(ct_i, celltype)
            print(len(posInds))
            if len(posInds) > 0:
                
                # G is indexed by label value, not by dataframe index
                posLabels = ds.df.label.iloc[posInds]
            
                degree = []
                neighbourTypes = []
                for i in posLabels:
                    neighbours = np.where(G[i,:])[0]
                    if neighbours[0] == 0:
                        neighbours = np.delete(neighbours,0)
                    unfilteredNeighbours = []
                    for v in neighbours:
                        if v in label2CellType:
                            unfilteredNeighbours.append(v)
                        
                    degree.append(len(unfilteredNeighbours))
                    neighbourTypes.extend([label2CellType[v] for v in unfilteredNeighbours])
                
                # OK, now we need to convert cluster labels into integers so that we can visualise
                # No, I don't know why they're strings either
                # Use the ordering of annotationDict to provide the integers
                celltype2Int = {list(annotationDict.keys())[v] : v for v in range(len(annotationDict))}
                neighbourTypesAsInts = [celltype2Int[v] for v in neighbourTypes]
                
                hist,bins = np.histogram(neighbourTypesAsInts,bins=np.arange(0,len(annotationDict)+1))
                fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15,6))
                ax[1].pie(hist/len(neighbourTypesAsInts),labeldistance=1.15,colors=colors)
                my_circle=plt.Circle( (0,0), 0.6, color='white')
                ax[1].add_artist(my_circle)
                my_circle=plt.Circle( (0,0), 0.55, color=colors[ct_i])
                ax[1].add_artist(my_circle)
                ax[1].set_title('n = ' + str(len(posInds)))
                
                ax[0].hist(degree,bins=np.arange(0,19),color=colors[ct_i])
                ax[0].set_title(annotationDict[celltype] + ', mean degree: ' + str(np.round(np.mean(degree),2)))
    
                plt.savefig(ds.pathToWriteOutput + ds.name + annotationDict[celltype] + '_adjacentCells.png',bbox_inches='tight')
                    
                plt.close()
                outputs[celltype] = {'degree':degree,
                                     'hist':hist,
                                     'neighbourTypesAsInts':neighbourTypesAsInts}
                
        saveFile = ds.pathToWriteOutput + ds.name + '_networkAnalysis_data.p'
        toSave = {'PickleFileGeneratedBy':'spatialstats.py',
                  'dilationValue':dilationValue,
                  'outputs':outputs}
                
        with open(saveFile,"wb") as fid:
            pickle.dump(toSave,fid)
        print('Network statistics completed')

def localClusteringHeatmaps(ds, df_annotations, clusteringToUse, clusterNames):
        print("Starting clustering heatmaps")   
        from utils_alt import returnAreaOfCircleInDomainAroundPoint
        # Identify where contributions to clustering/exclusion between two cell types are strongest or weakest, and save the resulting heatmap
        radiusOfInterest = 100 # microns
        # First we pre-calculate the area around each point (within the domain)
        vfunc_returnAreaOfCircleInDomainAroundPoint = np.vectorize(returnAreaOfCircleInDomainAroundPoint,excluded=['points'])
        areas = {}
        allPoints = {}
        for i, cluster in enumerate(df_annotations.ClusterNumber):
                print('Cluster',cluster)

                p0 = ds.points[ds.df[clusteringToUse] == cluster]
                if np.shape(p0)[0]>0:
                    p0_areas = vfunc_returnAreaOfCircleInDomainAroundPoint(index=np.arange(len(p0)), points=p0, r=radiusOfInterest, domainX=ds.domainX, domainY=ds.domainY)
                    allPoints[cluster] = p0
                    areas[cluster] = p0_areas
      
        # Now we use that to find the PCFs
        # gs = np.zeros(shape=(len(df_annotations.ClusterNumber), len(df_annotations.ClusterNumber), len(np.arange(0,maxR_mum,dr_mum))))
        outputs = {}
        for a, clusterA in enumerate(df_annotations.ClusterNumber):
                subOutputs = {}
                for b, clusterB in enumerate(df_annotations.ClusterNumber):
                        pair = [clusterA, clusterB]
                        print(pair)
                        
                        p_A = ds.points[ds.df[clusteringToUse] == clusterA]#allPoints[clusterA]# ds.points[ds.df[clusteringToUse] == int(pair[0])]    
                        p_B = ds.points[ds.df[clusteringToUse] == clusterB]#allPoints[clusterB]#ds.points[ds.df[clusteringToUse] == int(pair[1])]
                        if np.shape(p_A)[0]>0 and np.shape(p_B)[0]>0:
                            density_B = np.shape(p_B)[0]/(ds.domainX*ds.domainY)
                            
                            areas_A = areas[clusterA]
                            
                            distances_AtoB = cdist(p_A, p_B, metric='euclidean')
                            contributions = distances_AtoB <= radiusOfInterest
                            
                            BnearA_observed = np.sum(contributions,axis=1)/areas_A # observed per unit area
                            localPCFcontributions = BnearA_observed/density_B
                            
                            # Map PCF interpretation to [-1,1]
                            maxPCFthreshold = 5.0
                            minPCFthreshold = 1/maxPCFthreshold
                            
                            transformedLocalPCFcontributions = np.copy(localPCFcontributions)
                            transformedLocalPCFcontributions[transformedLocalPCFcontributions < minPCFthreshold] = minPCFthreshold
                            transformedLocalPCFcontributions[transformedLocalPCFcontributions > maxPCFthreshold] = maxPCFthreshold
                            
                            transformedLocalPCFcontributions[transformedLocalPCFcontributions<1] = -1/transformedLocalPCFcontributions[transformedLocalPCFcontributions<1]
                            # That gives us values in [-maxPCFthreshold,-1] U [1,maxPCFthreshold]
                            # Now map to [-1,1]
                            transformedLocalPCFcontributions[transformedLocalPCFcontributions<0] = (transformedLocalPCFcontributions[transformedLocalPCFcontributions<0]+1)/(maxPCFthreshold-1)
                            transformedLocalPCFcontributions[transformedLocalPCFcontributions>0] = (transformedLocalPCFcontributions[transformedLocalPCFcontributions>0]-1)/(maxPCFthreshold-1)
                            
                            # plt.figure()
                            # plt.hist(transformedPs,bins=100)
                            
                            # plt.figure(figsize=(12,9))
                            # plt.scatter(p_A[:,0],p_A[:,1],c=transformedLocalPCFcontributions,cmap='RdBu_r',vmin=-1,vmax=1)
                            # plt.colorbar()
                            # plt.gca().axis('equal')
                          
                            
                            savePngHeatmaps = True
                            if savePngHeatmaps:
                                
                                # Now we stick a Gaussian on top of each point 
                                # This turns out to be a bit annoying. 
                                # Solution: make a prototypical bell curve with height 1
                                # Add a new one centred at (x,y), scaled by transformedPs
                                containingCircleRadius = 150 # For bell curve, in microns
                                x, y = np.meshgrid(np.arange(-containingCircleRadius,containingCircleRadius+0.1,1),np.arange(-containingCircleRadius,containingCircleRadius+0.1,1))
                                dst = np.sqrt(x*x + y*y)
                                sigma = 50
                                bellCurve = np.exp(-( dst**2 / ( 2.0 * sigma**2 ) ) )
                                # plt.figure()
                                # plt.imshow(bellCurve)
                                
                                xrange = [-containingCircleRadius, ds.domainX+1+containingCircleRadius]
                                yrange = [-containingCircleRadius, ds.domainY+1+containingCircleRadius]
                                heatmap = np.zeros(shape=(xrange[1]-xrange[0],yrange[1]-yrange[0]))
                                
                                def addWeightedContribution(heatmap, weight, coordinate, xrange, yrange, bellCurve,containingCircleRadius):
                                    x0 = int(coordinate[0]) - xrange[0] - containingCircleRadius
                                    x1 = x0 + 2*containingCircleRadius + 1
                                    y0 = int(coordinate[1]) - yrange[0] - containingCircleRadius
                                    y1 = y0 + 2*containingCircleRadius + 1
                                    heatmap[x0:x1,y0:y1] = heatmap[x0:x1,y0:y1] + bellCurve*weight
                                    return heatmap
                                
                                for i in range(len(p_A)):
                                    coordinate = p_A[i,:]
                                    weight = transformedLocalPCFcontributions[i]
                                    heatmap = addWeightedContribution(heatmap, weight, coordinate, xrange, yrange, bellCurve, containingCircleRadius)
                                
                                #%
                                plt.figure(figsize=(12,9))
                                l = int(np.ceil(np.max([heatmap.min(),heatmap.max()])))
                                plt.imshow(heatmap.transpose(),cmap='RdBu_r',vmin=-l,vmax=l,origin='lower')
                                plt.xlim([0,ds.domainX])
                                plt.ylim([0,ds.domainY])
                                plt.colorbar()
                                ax = plt.gca()
                                ax.grid(False)
                                plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                                plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]] + '_localClusteringHeatmap.png',bbox_inches='tight')
    
                                plt.figure(figsize=(12,9))
                                plt.scatter(p_A[:,0],p_A[:,1],c='b')
                                plt.scatter(p_B[:,0],p_B[:,1],c='r')
                                plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                                plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]] + '_localClusteringHeatmapScatterplot.png',bbox_inches='tight')
                                plt.close('all')
                            subOutputs[clusterB] = transformedLocalPCFcontributions
                outputs[clusterA] = subOutputs
                        
                        
                
        saveFile = ds.pathToWriteOutput + ds.name + '_localClusteringHeatmap_data.p'
        toSave = {'PickleFileGeneratedBy':'spatialstats.py',
                  'transformedLocalPCFcontributions':outputs,
                  'radiusOfInterest':radiusOfInterest,
                  'maxPCFthreshold':maxPCFthreshold,
                  'minPCFthreshold':minPCFthreshold
                  }
        if savePngHeatmaps:
            toSave['bellcurve_containingCircleRadius'] =  containingCircleRadius
            toSave['bellcurve_sigma'] = sigma
                
        with open(saveFile,"wb") as fid:
            pickle.dump(toSave,fid)
        
        print("Local clustering heatmaps completed")   


if __name__ == "__main__":
    main()
