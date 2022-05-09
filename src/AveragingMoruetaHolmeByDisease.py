
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from utils_alt import dataset_filterSampleID, changeSomeElements
from statsmodels.stats.multitest import fdrcorrection
import os
import pickle
import argparse
import json

from spatialstats import preprocessing

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-i', '--pathToData', 
            help = 'Path to pickles',
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
            '-j', '--config', 
            help = 'json config',
            required = True
    )
    parser.add_argument(
            '-s', '--save_pickle', 
            help = 'json config',
            required = False,
            default=False
    )

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    pathToWriteOutput= args.output
    pathToData = args.pathToData
    cluster_annotations= args.cluster_annotations

    conf = json.loads(open(args.config).read())
    diseases= conf["disease_state"]
    diseasesToAverage = list(diseases.keys())
    

    clf= open(cluster_annotations).read().split("\n")
    allClusterNames= [x.strip().split("\t")[1] for x in clf[1:] if x != ""]
  

    conf = json.loads(open(args.config).read())


    for disease in diseasesToAverage:
        
        def getCFromSigma_inv(Sigma_inv):
            C = np.zeros(np.shape(Sigma_inv))
            for i in range(len(Sigma_inv)):
                for j in range(len(Sigma_inv)):
                    C[i,j] = -Sigma_inv[i,j] / np.sqrt(Sigma_inv[i,i] * Sigma_inv[j,j])
            return C
        
        Os = np.empty((0,0))
        Ns = np.empty((0,0,0))
        nSpecies = len(allClusterNames)
        for roi in diseases[disease]:
            arr= roi.split("_")
            #hack
            if arr[0]=="HEALTHY":
                arr[0]="HC"
            f= arr[0]+"_"+roi+"_"+arr[-2]+"_"+arr[-1]+"_MoruetaHolme_data.p"
            with open(os.path.join(pathToData, f),"rb") as fid:
                result = pickle.load(fid)

            O = result['O']
            N = result['Ns']
            clusterNames_MH = list(result['clusterNames_MH'])

            lenO = np.shape(O)[1]
            lenN = np.shape(N)[0]
            # Check saved names against full list of names to see what's missing and add the rows of zeros back in
            missing = [v for v in allClusterNames if v not in clusterNames_MH]
            if len(missing) > 0:
            # We need to insert the previously removed missing cells - removed for inverting on
            # a per ROI basis, but can be added to the big old list
                for missingVal in missing:
                    ind = allClusterNames.index(missingVal)
                    O = np.insert(O, ind, np.zeros((1,lenO)), 0)
                    N = np.insert(N, ind, np.zeros((1,lenO)), 1)
            if len(Os) == 0:
                # If we're at the start, initialise Os and Ns to the right size
                Os = np.empty((nSpecies,0))
                Ns = np.empty((lenN,nSpecies,0))

            Os = np.append(Os,O,axis=1)
            Ns = np.append(Ns,N,axis=2)

        clusterNames=allClusterNames
        if len(np.where(np.sum(Os,axis=1) == 1)) > 0:
            toExclude = list(np.where(np.sum(Os,axis=1)==0)[0])
            Os = np.delete(Os,toExclude,0)
            Ns = np.delete(Ns,toExclude,1)
            clusterNames= np.delete(np.asarray(clusterNames),toExclude,0)
            nSpecies = len(clusterNames)


        #%
        # Now continue as before from...
        # STEP THREE:
        # Calculate the covariance matrix Sigma for each of O and N
        # Use Sigma to calculate the partial correlation matrix C(O) and C(N) (nSpecies x nSpecies)
        # C_ij (M) = -Sigma^-1_ij (M) / sqrt(  Sigma^-1_ii (M) . Sigma^-1_jj (M)   )
        # Repeat for MANY resamples of N
        stub = os.path.join(pathToWriteOutput,disease+"-Average")
        O = Os
        Sigma_O = np.cov(O)
        Sigma_O_inv = np.linalg.inv(Sigma_O)
        C_O = getCFromSigma_inv(Sigma_O_inv)
        # Fix the main diagonal for plotting
        C_O_plot = C_O + 2*np.identity(nSpecies)
        plt.figure(figsize=(12,12))
        plt.imshow(C_O_plot,cmap='RdBu_r',vmin=-1,vmax=1)
        plt.colorbar()
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.savefig(stub + '_PartialCorrelation.png',bbox_inches='tight')
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
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.savefig(stub + '_StandardEffectSize.png',bbox_inches='tight')
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
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.colorbar()
        plt.savefig(stub + '_AdjacencyMatrix.png',bbox_inches='tight')
        
        plt.close('all')
        
        print("Quadrat neighbourhood analysis (Morueta-Holme) completed.")
        if (args.save_pickle):
            toSave = {'A':A,
                        'C_O':C_O,
                        'C_ns':C_ns,
                        #'counts':counts,
                        'clusterNames_MH':clusterNames_MH,
                        'Ns':Ns,
                        'nSpecies':nSpecies,
                        'O':O,
                        'p':p,
                        'p_star':p_star,
                        'rejected':rejected,
                        #'ds.name':ds.name,
                        'Sigma_O':Sigma_O,
                        'Sigma_O_inv':Sigma_O_inv,
                        'clusterNames':clusterNames}
            saveFile = stub + '-Average_Data.p'
            with open(saveFile,"wb") as fid:
                pickle.dump(toSave,fid)
        
if __name__ == "__main__":
    main()        