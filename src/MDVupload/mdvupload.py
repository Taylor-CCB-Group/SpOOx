import requests
import json
import os
import zipfile
import gzip
import sys
import pandas
import argparse
import tarfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--name',help="The name of the project",required=True)
    parser.add_argument('-u','--url', help ="the url of the MDV website",default="https://mdv.molbiol.ox.ac.uk")
    parser.add_argument('-t','--token' , help="The token obtained from the website",required=True)
    parser.add_argument('-d','--directory', help="the director where SpOOx was run", default =".")
    parser.add_argument('-g','--genome', help="genome build - not used but required - default hg19", default ="hg19")
    parser.add_argument('-i','--images',nargs = '+', help="specify which images to use")
    parser.add_argument('-dr','--dry_run', help="only create files-do not upload",default =False,type=bool)
    parser.add_argument("-td","--tempdir",help="the temporary directory to put the tar ball to upload. Default is _temp in the run folder")
    
    

    args = parser.parse_args()
    basedir=args.directory
    #config directory is where temporary files will be written
    conf_dir = os.path.join(basedir,"_temp")
    if args.tempdir:
        conf_dir = args.tempdir
    if not os.path.exists(conf_dir):
        os.mkdir(conf_dir)
    #all files will be added to single tar file that will be streamed
    tar_file = os.path.join(conf_dir,"data.tar")
    tar = tarfile.open(tar_file,"w")

    #get the url and token 
    url = args.url +"/hyperion/create_hyperion_project"
    token = args.token
  
    #get the marker and metadata files
    marker_file=os.path.join(basedir,"signalextraction","markers.tsv")
      
    if not os.path.exists(marker_file):
        print ("Can't find marker panel")
        sys.exit()

    sample_file = os.path.join(basedir,"metadata.tsv")
    if not os.path.exists(sample_file):
        print ("Can't find sample file")
        sys.exit()

    #work out the markers
    conf={
        "dimension_reduction":{
            "UMAP":["UMAP_1","UMAP_2","UMAP_3"]
        },
        "name":args.name,
        "images_to_include":["DNA1_Ir191"],
        "default_image":"cellmask",
        "clustering":["harmony_pheno_cluster"],
        "position_fields":["x","y"],
        "create_files_only":args.dry_run,
        "micro_per_pixel":1,
        "cell_identifier":"cellID",
        "sample_identifier":"sample_id",
        "sample_fields":["sample_name","ROI","condition"],
        "cell_data":["area","eccentricity","major_axis_length","minor_axis_length","perimeter"]
    }

    if args.images:
        conf["images_to_include"]=args.images

    
   
   
    mdf =pandas.read_csv(marker_file,sep="\t")
    conf["cluster_markers"]= list(mdf[mdf.clustering==1]["marker_name"])
    conf["all_markers"]= mdf["marker_name"].tolist()
    conf["nucleus_markers"]=list(mdf[mdf.nucleus==1]["marker_name"])
    

    
    data_file = os.path.join(basedir,"clustering","clusters.txt")
    if not os.path.exists(data_file):
        print ("Can't find data file - {}".format(data_file))
        sys.exit()

   

    #get all the fields required
    dr_fields=[]
    for dr in list(conf["dimension_reduction"].values()):
        dr_fields+=dr
    all =conf["all_markers"]+conf["position_fields"]+conf["clustering"] + dr_fields +\
        [conf["sample_identifier"]]+[conf["cell_identifier"]]+conf["sample_fields"]+conf["cell_data"]
    if conf.get("other_fields"):
        for item in conf["other_fields"]:
            all.append(item["name"])

    
    t_f= os.path.join(conf_dir,"t_data.gz")
    o= gzip.open(t_f,'wt')
   
  

    #parse the data file
    sample_id_index =None
    samples= set()
    headers=None
    header_to_index={}
    missing_fields=[]
    annotation_index = None
    annotation_info=None
   

    if conf.get("annotations"):
        anno_file = conf["annotations"]["file"]
        if not os.path.exists(anno_file):
            print ("Can't find annotation file:{}".format(anno_file))
            sys.exit()
        annotation_info={}
        first=True
        with open(anno_file) as f:
            for line in f:
                if first:
                    first=False
                    continue
                arr=line.strip().split("\t")
                #convert to cl.. format
                if len(arr[0])==1:
                    arr[0]= "cl0"+arr[0]
                else:
                    arr[0]="cl"+arr[0]      
                annotation_info[arr[0]]=arr[1].strip()
        if not conf.get("other_fields"):
            conf["other_fields"]=[]
        conf["other_fields"].append({
            "name":conf["annotations"]["name"],
            "type":"text"
        })


        
    print ("parsing and compressing data file")
    sti = conf.get("samples_to_include")
    if sti:
        sti= set(sti)
    read=0
    with open(data_file) as f:
        first=True
        for line in f:
            arr = line.strip().split("\t")
            if first:
                headers=arr

                first =False
                for field in all:
                    try:
                        #r changes - to .
                        if "-" in field:
                            index = headers.index(field.replace("-","."))
                        else:
                            index= headers.index(field)
                        header_to_index[field]=index
                    except:
                        missing_fields.append(field)
                if len(missing_fields)>0:
                     print("The following fields are missing in the data file: {}".format(",".join(missing_fields)))
                     sys.exit()
                sample_id_index= headers.index(conf["sample_identifier"])
                o.write("\t".join(all))
                if annotation_info:
                    try:
                        annotation_index=arr.index(conf["annotations"]["based_on"]) 
                    except:
                        print("annotation field not present")
                        sys.exit()
                    o.write("\t{}".format(conf["annotations"]["name"]))
                o.write("\n")    
                
                
                continue

            
            read+=1
            #only include certain samples
            if sti and not arr[sample_id_index] in sti:
                if read%5000==0:
                    print ("processed {} lines".format(read))
                continue
           
            to_write=[]
            for field in all:
                to_write.append(arr[header_to_index[field]])
            if annotation_info:
                anno = annotation_info.get(arr[annotation_index],"NA")
                to_write.append(anno)
            o.write("\t".join(to_write)+"\n")
           
            if read%5000==0:
                print ("processed {} lines".format(read))
                      
    o.close()
    tar.add(t_f,arcname="data.gz")
    os.remove(t_f)


    #the samples file
    sdf = pandas.read_csv(sample_file,sep="\t")
    if sti:
        sdf=sdf[sdf["sample_id"].isin(sti)]
    samples= list(sdf["sample_id"])
    sample_file = os.path.join(conf_dir,"t_sample.tsv")
    sdf.to_csv(sample_file,index=False,sep="\t")
    tar.add(sample_file,arcname="samples.tsv")
    os.remove(sample_file)
   
    #zip all the images
    print ("zipping images")
   
   
    images_zip = os.path.join(conf_dir,"temp_images.zip")
    zip_file = zipfile.ZipFile(images_zip, 'w')

    #add the cell masks
    for sample in samples:
        f = os.path.join(basedir,"deepcell",sample,"bandwmask.png")
        imname = sample+"_"+"cellmask.png"
        zip_file.write(f,arcname=imname, compress_type=zipfile.ZIP_DEFLATED)


    iti = conf.get("images_to_include")
    imgfolder = os.path.join(basedir,"pngs","img")
    if os.path.exists(imgfolder):
        images = os.listdir(imgfolder)
        for im in images:        
            if iti:
                accept=False
                for iname in iti:
                    if iname in im:
                        accept=True
                        break
                if not accept:
                    continue
            if im.endswith("_Mask.png"):
                continue
            for sample in samples:        
                if im.startswith(sample):
                    zip_file.write(os.path.join(imgfolder,im),arcname=im, compress_type=zipfile.ZIP_DEFLATED)
                    break
    zip_file.close()

    tar.add(images_zip,arcname="main_images.zip")
    os.remove(images_zip)


    #spatial stats
    ss_dir = os.path.join(basedir,"spatialstats")
    ss_sum = os.path.join(ss_dir,"summary.tsv")
    if os.path.exists(ss_sum):
        tar.add(ss_sum,arcname="spatialstats.tsv")
        tar.add(os.path.join(ss_dir,"pcf_images.zip"),arcname="pcf_images.zip")
        tar.add(os.path.join(ss_dir,"annotations.tsv"),arcname="annotations.tsv")
        lchm_file= os.path.join(ss_dir,"lhm_images.zip")
        if os.path.exists(lchm_file):
            tar.add(lchm_file,arcname="lchm_images.zip")

    #average spatial stats
    ssa_dir = os.path.join(basedir,"spatialstats_average")
    ssa_sum = os.path.join(ssa_dir,"summary.tsv")
    if os.path.exists(ssa_sum):
        tar.add(ssa_sum,arcname="spatialstats_av.tsv")
        tar.add(os.path.join(ssa_dir,"pcf_images","av_pcf_images.zip"),arcname="av_pcf_images.zip")
        tar.add(os.path.join(basedir,"cond_to_sampleid.json"),arcname="cond_to_sampleid.json")

    #ome-tiffs
    ot_dir = os.path.join(basedir,"pyramidal_ometiff")
    if os.path.exists(ot_dir):
        for f in os.listdir(ot_dir):
            if f.endswith("ome.tiff"):
                tar.add(os.path.join(ot_dir,f),f)

    #save the config
    mlvconfig= os.path.join(conf_dir,"temp_config.json")
    o = open(mlvconfig,"w")
    o.write(json.dumps(conf,indent=2))
    o.close()

    tar.add(mlvconfig,arcname="config.json")
    tar.close()
     
    if not conf.get("create_files_only"):
        print ("uploading files")
        headers= {
            "x-access-tokens":token,
            "mdv-name":args.name,
            "mdv-genome":args.genome
        }
        with open(tar_file, 'rb') as f:
            r = requests.post(url, headers= headers,data=f)
        
        if not r.status_code==200:
            print ("there was a problem connecting to the server")
        if r.text == "ok":
            print("uploaded")
        elif r.text=="permission denied":
            print("Permission denied - do you have the correct token?")
        elif r.text=="error":
            print("There was an error on the server")
        os.remove(tar_file)

      
    else:
        print ("Files generated")









if __name__ == "__main__":
    main()
    