
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils_alt import  getPCFContributionsWithinGrid
from statsmodels.stats.multitest import fdrcorrection
import os
import pickle
import argparse
import json
import pandas as pd
import zipfile
from scipy.stats import combine_pvalues,norm
import math



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
                '-cl', '--clusteringIdColumn', 
                help = 'Name of column in the file specified by pathToData to use as the clustering label. Should correspond with the annotations in cluster_annotations',
                required = True
        )

        parser.add_argument(
                '-p', '--pathToSS', 
                help = 'The path to the spatial stats data',
                required = True
        )
        parser.add_argument(
                '-j', '--config', 
                help = 'json config',
                required = True
        )

        parser.add_argument(
                '-pi', '--processImages', 
                help = 'rename and zip pcf images for upload to mdv',
                required = False,
                default=False,
                type=bool
        )

        args = parser.parse_args()

        pathToData = args.pathToData
        clusteringIdColumn = args.clusteringIdColumn

        annotations = pd.read_csv(os.path.join(args.pathToSS,"annotations.tsv"),sep="\t")
     
        
        pathToSaveFigures = args.output+"/"
        os.makedirs(pathToSaveFigures,exist_ok=True)

       
        conf = json.loads(open(args.config).read())

        diseases= conf["conditions"]
        state_name = conf.get("name","state")
        print("state_name:"+state_name)
        if (os.path.exists(os.path.join(args.pathToSS,"contactingcellnetwork"))):
            ccn_folder = os.path.join(pathToSaveFigures,"contactingcellnetwork")
            os.makedirs(ccn_folder,exist_ok=True)
            average_ccn(ccn_folder,diseases,args.pathToSS,annotations)
 
        
        pc_folder = os.path.join(pathToSaveFigures,"paircorrelationfunction")
        os.makedirs(pc_folder,exist_ok=True)
        average_pcf(pc_folder,diseases,args.pathToSS,pathToData,clusteringIdColumn, annotations)

        mh_folder = os.path.join(pathToSaveFigures,"morueta-holme")
        os.makedirs(mh_folder,exist_ok=True)
        average_mh(mh_folder,diseases,args.pathToSS,annotations)

        nt_folder = os.path.join(pathToSaveFigures,"networkstatistics")
        os.makedirs(nt_folder,exist_ok=True)
        average_network(nt_folder,diseases,args.pathToSS,annotations)

        if args.processImages:
            im_folder = os.path.join(pathToSaveFigures,"pcf_images")
            os.makedirs(im_folder,exist_ok=True)
            process_mdv_images(im_folder,diseases,pc_folder)

        make_summary_table(pathToSaveFigures,diseases,args.processImages,state_name)



def average_ccn(ccn_folder,diseases,pathToSS,annotations):
    summary_file = os.path.join(pathToSS,"summary.tsv")
    df= pd.read_csv(summary_file,sep="\t")


    cell_numbers= {f"{k1}|{k2}":v1 for k1,k2,v1 in zip(df["sample_id"],df["Cell Type 1"],df["cell 1 number"]) }
    cell_inters = [f"{k1}|{k2}" for k1,k2 in zip(df["Cell Type 1"],df["Cell Type 2"]) ]

    cell_inters=set(cell_inters)

    z_scores = {f"{k1}|{k2}|{k3}":v1 for k1,k2,k3,v1 in zip(df["sample_id"],df["Cell Type 1"],df["Cell Type 2"],df["contact_zscore"])}

    conds = diseases
    res={}
    for c in conds:
        res_to_index={}
        all_ps=[]
        all_zs=[]
        index=0
        for i in cell_inters:
            zsc=[]
            weights=[]
            pvals=[]
            for s in conds[c]:
                zs = z_scores.get(f"{s}|{i}")
                if zs==None or zs== "None" or math.isnan(zs):
                    continue
                pv = norm.sf(abs(zs)) 
                pvals.append(pv)
                zsc.append(zs)
                c1= i.split("|")[0]
                weights.append(cell_numbers[f"{s}|{c1}"])
            if len(pvals)==0:
                continue
            stat,pval= combine_pvalues(pvals,"stouffer",weights)
            ws = sum(weights)
            wzc=0
            for n in range(len(zsc)):
                wzc+=(zsc[n]*weights[n])/ws
            all_ps.append(pval)
            all_zs.append(wzc)
            res_to_index[i]=index
            index+=1
        bl,fdr = fdrcorrection(all_ps)
        for r in res_to_index:
            ind= res_to_index[r]
            res[f"{c}|{r}"]=[fdr[ind],all_zs[ind]]        
    o=open(f"{ccn_folder}/summary.tsv","w")
    o.write("state\tCell Type 1\tCell Type 2\tcontact_zscore\tcontact_pval_fdr\n")

    for r in res:
        state,ct1,ct2 = r.split("|")
        o.write(f"{state}\t{ct1}\t{ct2}\t{res[r][1]}\t{res[r][0]}\n")
    o.close()


def make_summary_table(dir,diseases,add_pcf_images,state_name):
    def get_df(f):
        dframe= pd.read_csv(f,sep="\t")
        dframe = dframe.rename(columns={'state': state_name})
        return dframe.set_index([state_name,"Cell Type 1","Cell Type 2"])

    dfs=[ get_df(os.path.join(dir,"paircorrelationfunction",x,"summary.tsv")) for x in diseases]

    pc_df= pd.concat(dfs)

    mh_df = get_df(os.path.join(dir,"morueta-holme","summary.tsv"))
    nt_df=  get_df(os.path.join(dir,"networkstatistics","summary.tsv"))
   

    all_tables= [pc_df,mh_df,nt_df]
    ccn_file =os.path.join(dir,"contactingcellnetwork","summary.tsv")
    if os.path.exists(ccn_file):
        ccn_df  = get_df(os.path.join(dir,"contactingcellnetwork","summary.tsv"))
        all_tables.append(ccn_df)

    if add_pcf_images:
        all_tables.append(get_df(os.path.join(dir,"pcf_images","summary.tsv")))

    table = pd.concat(all_tables,axis=1,join="outer")

    table.to_csv(os.path.join(dir,"summary.tsv"),sep="\t")

   


def process_mdv_images(im_folder,diseases,pathToSS):
    images_zip  = os.path.join(im_folder,"av_pcf_images.zip")
    zip_file = zipfile.ZipFile(images_zip, 'w')
    index=0
    pathToPCF= os.path.join(pathToSS,"")
    o=open(os.path.join(im_folder,"summary.tsv"),"w")
    o.write("state\tCell Type 1\tCell Type 2\tpcf_image_key\n")
    for state in diseases:
        
        
        ifolder  = os.path.join(pathToPCF,state)
        ims= os.listdir(ifolder)
        for im in ims:
            if im.endswith(".png"):
                t= im.replace("_averageWithCI95.png","")
                cells = t.split("-to-")     
                f= os.path.join(ifolder,im)
                zip_file.write(f,arcname="im{}.png".format(index), compress_type=zipfile.ZIP_DEFLATED)
                o.write("{}\t{}\t{}\t{}\n".format(state,cells[0],cells[1],index))
                index+=1
                    
    zip_file.close()
    o.close()

def average_mh(mh_folder,diseases,pathToSS,annotations):
    
    pathToWriteOutput= mh_folder
    mhPath = os.path.join(pathToSS,"morueta-holme")
    diseasesToAverage = list(diseases.keys())
    allClusterNames= list(annotations["Annotation"])
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
            fl = os.path.join(mhPath, f)
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

        
      
        
        alpha = 0.05
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
      
    main_table = open(os.path.join(pathToWriteOutput,"summary.tsv"),"w")
    main_table.write("state\tCell Type 1\tCell Type 2\tMH_PC\tMH_SES\tMH_FDR\n")
    for d in main_data:
        for intera in main_data[d]:
            ns = intera.split("|")
            info = main_data[d][intera]
            main_table.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(d,ns[0],ns[1],info.get("PC","NA"),info.get("SES","NA"),info.get("FDR","NA")))
    main_table.close() 


def average_network(nt_folder,diseases,pathToSS,annotations):
    dir = os.path.join(pathToSS,"networkstatistics")
    fs= os.listdir(dir)
   
    cl_length = len(annotations["Annotation"])
    id_to_clus= {x[0]:x[1]  for x in zip(annotations["ClusterNumber"],annotations["Annotation"])}
    clus_li = list(annotations["Annotation"])
    data={}
    for f in fs:
        if f.endswith(".p"):
            sample_id=f.replace("_networkAnalysis_data.p","")
            fi = os.path.join(dir,f)
            with open (fi,"rb") as fo:
                x=pickle.load(fo)
                for cl in x["outputs"]:
                    index=0
                    cellsn=[0] * cl_length
                    #degree is number of cells in contact with cell each cellA
                    #hence length is number of A cells
                    n_cells=  len(x["outputs"][cl]["degree"])
                    #total number of cells surrounding all cell As
                    n_surrounding_cells= sum(x["outputs"][cl]["degree"])
                    for ncells in x["outputs"][cl]["degree"]:
                        neighbours = x["outputs"][cl]["neighbourTypesAsInts"][index:index+ncells]
                        #start with no cells in conatct
                        incontact=[False]*cl_length
                        for n in neighbours:
                            #at least one in contact
                            incontact[n]=True
                        for i,n in enumerate(incontact):
                            if n:
                                cellsn[i]+=1
                        index+=ncells
                    #md = [x/n_cells for x in cellsn]
                    cl1name =id_to_clus[cl]
                    for i,v in enumerate(cellsn):
                        ntw= x["outputs"][cl]["hist"][i]
                        cl2name = clus_li[i]
                        data["{}|{}|{}".format(sample_id,cl1name,cl2name)] = [v,(v/n_cells)*100,ntw,(ntw/n_surrounding_cells)*100,ntw/n_cells]
    sf= open(os.path.join(nt_folder,"summary.tsv"),"w")
    sf.write("state\tCell Type 1\tCell Type 2\tcontacts\t%contacts\tNetwork\tNetwork(%)\tmean degree\n")
    for cond in diseases:
        for name1 in clus_li:
            for name2 in clus_li:
                totals=[0,0,0,0,0]
                for sid in diseases[cond]:
                    info = data.get("{}|{}|{}".format(sid,name1,name2))
                    if info:
                        for n in range(0,5):
                            totals[n]+=info[n]
                
                for n in range(0,5):
                    totals[n]= round(totals[n]/len(diseases[cond]),3)
                sf.write("{}\t{}\t{}\t".format(cond,name1,name2))
                sf.write("\t".join([str(x) for x in totals]))
                sf.write("\n")
    
    sf.close()


def average_pcf(pc_folder,diseases,pathToSS,pathToData,clusteringIdColumn, annotations):
    diseasesToAverage = list(diseases.keys())
    for d in diseasesToAverage:
        os.makedirs(os.path.join(pc_folder,d),exist_ok=True)
        os.makedirs(os.path.join(pc_folder,d,"Pickles"),exist_ok=True)


    
    for disease in diseasesToAverage:
        summary_data=[]
        rois = diseases[disease]
        datasets = preprocessing(pathToData, rois)
        clusterNames = {annotations.ClusterNumber.iloc[v] : annotations.Annotation.iloc[v].strip() for v in range(len(annotations))}
        #calculate average number of cells
        cell_totals= {k:0  for k in clusterNames}
        for ds in datasets:
            vs = dict(ds.df[clusteringIdColumn].value_counts())
            for k in vs:
                cell_totals[k]+=vs[k]
        cell_averages = {k:int(v/len(datasets)) for k,v in cell_totals.items()}

        pfiles={}
        for roi in rois:
            pfile = os.path.join(pathToSS,"paircorrelationfunction",f'{roi}_PCF_data.p')     
            with open(pfile,"rb") as fid:
                pfiles[roi] = pickle.load(fid)


        for a,cla in enumerate(annotations.Annotation):
            for b,clb in  enumerate(annotations.Annotation):
                
                allContributions = []
                allPCFs = []
                all_nRectangles = []
                all_rectContributions = []
                all_rectNs = []
                
                #summary line
                cla=annotations.ClusterNumber.iloc[a]
                clb=annotations.ClusterNumber.iloc[b]
                pair = [cla, clb]
                line  = "{}\t{}\t{}\t{}\t{}".format(disease,clusterNames[cla],clusterNames[clb],cell_averages[cla],cell_averages[clb])

                for ds in datasets:
                    domainX = ds.domainX
                    domainY = ds.domainY
                                
                    # First we pre-calculate the area around each point (within the domain)
                    p_A = ds.points[ds.df[clusteringIdColumn] == cla]
            
                    #get the pcf data
                    pdata= pfiles[ds.sample]

                    radii = pdata["radii"]
                    g = pdata["gs"][a,b]
                    contributions= pdata["cs"][a][b]
                    if  contributions is None:
                        continue
                                 
                    allContributions.append(contributions)
                    allPCFs.append(g)   
                    nRectangles, rectContributions, rectNs = getPCFContributionsWithinGrid(contributions, domainX, domainY, p_A)
                    all_nRectangles.append(nRectangles)
                    all_rectContributions.append(rectContributions)
                    all_rectNs.append(rectNs)
                        
                        
                if len(allContributions)==0:
                        print ("No data for {} in {}".format(pair,disease))
                        summary_data.append(line+"\tND\tND\tND\tND\tND\tND\tND\tND\tND\n")
                        continue
            
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
                figpath= os.path.join(pc_folder,disease,f'{clusterNames[pair[0]]}-to-{clusterNames[pair[1]]}_averageWithCI95.png')
                
                plt.savefig(figpath)
                

                plt.close()
                to_save={
                    "radii":radii,
                    "PCF_mean":PCF_mean_fromContributions,
                    "PCF_min":PCF_min_fromContributions,
                    "PCF_max":PCF_max_fromContributions
                }
                p_file= os.path.join(pc_folder,disease,"Pickles",f'{clusterNames[pair[0]]}-to-{clusterNames[pair[1]]}.p')
                with open(p_file,"wb") as fid:
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

                print(f'PCF {disease}:{pair} complete')
                    
        sf = open(os.path.join(pc_folder,disease,"summary.tsv"),"w")
        sf.write("state\tCell Type 1\tCell Type 2\tmean cell 1 number\tmean cell 2 number")
        sf.write("\tg(r)max\tgr10\tgr20")
        sf.write("\tgr10 PCF upper\tgr10 PCF lower\tgr20 PCF upper\tgr20 PCF lower")
        sf.write ("\tPCF min intersect\tPCF max intersect\n")
        for line in summary_data:
            sf.write(line)
        sf.close()            
                        


if __name__ == "__main__":
    main()   
    
















