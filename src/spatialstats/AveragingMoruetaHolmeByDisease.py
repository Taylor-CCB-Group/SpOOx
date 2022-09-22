
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
            help = 'save the pickle file',
            required = False,
            default=False
    )

    parser.add_argument(
            '-f', '--fdr', 
            help = 'fdr',
            required = False,
            type=float,
            default=0.05
    )

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    pathToWriteOutput= args.output
    pathToData = args.pathToData
    cluster_annotations= args.cluster_annotations

    conf = json.loads(open(args.config).read())
    diseases= conf["conditions"]
    diseasesToAverage = list(diseases.keys())

    if not os.path.exists(pathToWriteOutput):
        os.mkdir(pathToWriteOutput)
    

    clf= open(cluster_annotations).read().split("\n")
    allClusterNames= [x.strip().split("\t")[1] for x in clf[1:] if x != ""]
  

    conf = json.loads(open(args.config).read())

  
    main_data={}
    for disease in diseasesToAverage:
        ddata={}
        for i1,n1 in enumerate(allClusterNames):
            for i2,n2 in enumerate(allClusterNames):
                ddata["{}|{}".format(n1,n2)]={}
        main_data[disease]=ddata


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
            f= roi+"_MoruetaHolme_data.p"
            fl = os.path.join(pathToData, f)
            if not os.path.exists(fl):
                continue
            with open(fl,"rb") as fid:
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
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.savefig(stub + '_PartialCorrelation.png',bbox_inches='tight')
        plt.close()

        
        o = open(stub+"_correlation_matrix.txt","w")

        for i1,n1 in enumerate(clusterNames):
            for i2,n2 in enumerate(clusterNames):
                v= C_O_plot[i1,i2]
                o.write("{}|{}".format(n1,n2)+"\t")
                o.write(str(v)+"\n")              
                main_data[disease]["{}|{}".format(n1,n2)]["PC"]=v

        o.close()
        
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
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.savefig(stub + '_StandardEffectSize.png',bbox_inches='tight')
        plt.close()

        for i1,n1 in enumerate(clusterNames):
            for i2,n2 in enumerate(clusterNames):   
                    main_data[disease]["{}|{}".format(n1,n2)]["SES"]=SES[i1,i2]


        
        # STEP SIX: 
        # Calculate a 2-tailed p-value for each pair of species
        p = np.zeros(shape=(nSpecies,nSpecies))
        p_list = []
        p_index = []
        for i in range(nSpecies):
            for j in range(nSpecies):
                p[i,j] = sum(np.abs(C_O[i,j]) < np.abs(C_ns[:,i,j]))/np.shape(C_ns)[0]
                if j < i:     
                    p_list.append(p[i,j])
                    p_index.append([i,j])

        
      
        
        alpha = args.fdr
        rejected, p_star = fdrcorrection(p_list, alpha, method='indep')
       
        A = np.zeros(shape=(nSpecies,nSpecies))
        A[:,:] = np.nan # Initialise A with Nans
        for k in range(len(p_index)):
            i = p_index[k][0]
            j = p_index[k][1]
            cn1= clusterNames[i]
            cn2 = clusterNames[j]
            main_data[disease]["{}|{}".format(cn1,cn2)]["FDR"]=p_star[k]
            main_data[disease]["{}|{}".format(cn2,cn1)]["FDR"]=p_star[k]
            
            if p_star[k] < alpha:  
                A[i,j] = SES[i,j]
        
        # create other half of matrix
        # simplest A= A + A.T (but adding value to NaN is NaN)
        B=A.T
        na= np.isnan(A)
        nb = np.isnan(B)
        A[na]=0
        B[nb]=0
        A +=B
        na &= nb
        A[na] = np.nan


        #write out txt file
        o = open(stub+"_adj_matrix_{}.txt".format(alpha),"w")
        o.write("\t"+"\t".join(clusterNames)+"\n")
        for i,n in enumerate(clusterNames):
            line = [str(x) for x in A[i]]
            o.write(n+"\t"+"\t".join(line)+"\n")
        o.close()

        o = open(stub+"_adj_groups_{}.txt".format(alpha),"w")
        for i1,n in enumerate(clusterNames):
            o.write(n+"\n")
            for i2,n in enumerate(clusterNames):
                v= A[i1,i2]
                if not np.isnan(v):
                    o.write(n+"\t"+str(v)+"\n")
            o.write("\n")
        o.close()


      

        



        plt.figure(figsize=(12,12))
        plt.imshow(A,cmap='RdBu_r',vmin=-15,vmax=15)
        plt.xticks(ticks=np.arange(0,nSpecies),labels=clusterNames,rotation=90)
        plt.yticks(ticks=np.arange(0,nSpecies),labels=clusterNames)
        plt.colorbar()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.savefig(stub + '_AdjacencyMatrix_{}.png'.format(alpha),bbox_inches='tight')
        
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
    main_table = open(os.path.join(pathToWriteOutput,"all_data.txt"),"w")
    main_table.write("state\tcell1\tcell2\tpartial_correlation\tstandard_effect_size\tFDR\n")
    for d in main_data:
        for intera in main_data[d]:
            ns = intera.split("|")
            info = main_data[d][intera]
            main_table.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(d,ns[0],ns[1],info.get("PC","NA"),info.get("SES","NA"),info.get("FDR","NA")))
    main_table.close()    
if __name__ == "__main__":
    main()