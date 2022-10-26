import argparse
import pickle
import os
import pandas
import json

#content of the pickle files
#for each cluster
# hist - the number of cells of each type surrounding this type
# degree - the number of cells surrounding each cell
# neighbourTypesASInts - the actual type of cell syrrounding each (length of this array = sum (degree))

#dir = the directory containing the network stats
#the cluster ids to names
# 


def summarize_networksats(dir, id_to_clus, cond_to_sample):
    fs= os.listdir(dir)
    cl_length = len(id_to_clus)
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
                    
                    n_cells=  len(x["outputs"][cl]["degree"])
                    n_surrounding_cells= sum(x["outputs"][cl]["degree"])
                    for ncells in x["outputs"][cl]["degree"]:
                        neighbours = x["outputs"][cl]["neighbourTypesAsInts"][index:index+ncells]
                        incontact=[False]*cl_length
                        for n in neighbours:
                            incontact[n]=True
                        for i,n in enumerate(incontact):
                            if n:
                                cellsn[i]+=1
                        index+=ncells
                    #md = [x/n_cells for x in cellsn]
                    cl1name =id_to_clus[cl]
                    for i,v in enumerate(cellsn):
                        ntw= x["outputs"][cl]["hist"][i]
                        cl2name = id_to_clus["cl0"+str(i+1) if i<9 else "cl"+str(i+1)]
                        data["{}|{}|{}".format(sample_id,cl1name,cl2name)] = [v,(v/n_cells)*100,ntw,(ntw/n_surrounding_cells)*100,ntw/n_cells]
    if not cond_to_sample:
        return data
    final_data={}
    for cond in cond_to_sample:
        for name1 in id_to_clus.values():
            for name2 in id_to_clus.values():
                totals=[0,0,0,0,0]
                for sid in cond_to_sample[cond]:
                    info = data.get("{}|{}|{}".format(sid,name1,name2))
                    if info:
                        for n in range(0,5):
                            totals[n]+=info[n]
                
                for n in range(0,5):
                    totals[n]= round(totals[n]/len(cond_to_sample[cond]),3)
                final_data["{}|{}|{}".format(cond,name1,name2)]=totals
    
    return final_data




def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
            '-p', '--pathToPCF', 
            help = 'Path to pickles',
            required = True
    )
    parser.add_argument(
            '-n', '--pathToNS', 
            help = 'Path to pickles',
            required = True
    )
    
    parser.add_argument(
            '-m', '--pathToMH', 
    help = 'Path to write all the outputs to.',
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
   
    args = parser.parse_args()

    anno = pandas.read_csv(args.cluster_annotations,sep="\t")
    id_to_clus= {"cl0"+str(x[0]) if x[0]<10 else "cl"+str(x[0]):x[1] for x in zip(anno["ClusterNumber"],anno["Annotation"])}
    cond_to_sample  = json.loads(open(args.config).read())["conditions"]
    sample_to_cond={}
    for c in cond_to_sample:
        for s in cond_to_sample[c]:
            sample_to_cond[s]=c

    net_data = summarize_networksats(args.pathToNS,id_to_clus,cond_to_sample)
   
    mh_data={}
    with  open(os.path.join(args.pathToMH,"all_data.txt")) as f:
        first=True
        for line in f:
            if first:
                first=False
                continue
            arr= line.strip().split("\t")
            mh_data["{}|{}|{}".format(arr[0],arr[1],arr[2])]=arr[3:]
    o= open(args.output,"w")
    with open(os.path.join(args.pathToPCF,"summary.txt")) as f:
        first=True
        for line in f:
            line=line.strip()
            if first:
                o.write(line+"\tMH_FDR\tMH_PC\tMH_SES\tcontacts\t%contacts\tNetwork\tNetwork(%)\tmean degree\n")
                first= False
                continue
            o.write(line)
            arr=line.split("\t")
            key =  "{}|{}|{}".format(arr[0],arr[1],arr[2])
            o.write("\t"+"\t".join(mh_data[key]))
            o.write("\t"+"\t".join([str(x) for x in net_data[key]]))
            o.write("\n")
    o.close()





 


if __name__ == "__main__":
    main()