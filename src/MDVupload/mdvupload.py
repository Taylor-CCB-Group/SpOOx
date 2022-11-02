import requests,json,os,zipfile,gzip,yaml,sys,pandas,re
from tifffile import TiffFile
import math
import subprocess
import random,string

def main(config):
    with open(config) as file:
        try:
            conf = yaml.safe_load(file)
        except Exception as e:
            print("Can't open config -{}".format(config))
            sys.exit()
  
    #config directory is where temporary files will be written
    conf_dir = conf.get("output_dir","'")

    #get the url and token then remove them
    url = conf["url"]
    token = conf["token"]
    del conf["url"]
    del conf["token"]

    #get the marker and metadata files
    marker_file=conf["marker_file"]
 

      
    if not os.path.exists(marker_file):
        print ("Can't find marker panel")
        sys.exit()

    #work out the markers
   
    mdf =pandas.read_csv(marker_file,sep="\t")
    conf["cluster_markers"]= list(mdf[mdf.clustering==1]["marker_name"])
    conf["all_markers"]= mdf["marker_name"].tolist()
    conf["nucleus_markers"]=list(mdf[mdf.nucleus==1]["marker_name"])
    

    
    data_file = conf["data_file"]
    if not os.path.exists(data_file):
        print ("Can't find data file - {}".format(data_file))
        sys.exit()

    if conf.get("image_name_only"):
        df = pandas.read_csv(data_file,sep="\t")
        df["sample_id"]=df.apply (lambda row: "_".join(row["Image Name"].split("_")[:5]), axis=1)
        df["sample_name"]= df.apply(lambda row: "_".join(row["Image Name"].split("_")[:3]),axis=1)
        df["condition"]= df.apply(lambda row: row["Image Name"].split("_")[0],axis=1)
        df["ROI"] = df.apply(lambda row: "_".join(row["Image Name"].split("_")[3:5]),axis=1)
        t_data_file = "temp_data_file.tsv"
        df.to_csv(t_data_file,sep="\t",index=False)
        data_file=t_data_file

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
    #o=open(t_f,"w")
  

    #parse the data file
    sample_id_index =None
    samples= set()
    headers=None
    header_to_index={}
    missing_fields=[]
    annotation_index = None
    annotation_info=None
    sample_info={}
    

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
                        index= headers.index(field)
                        header_to_index[field]=index
                    except:
                        if field=="Va7-2":
                            header_to_index[field]=headers.index("Va7.2")
                        else:
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

            #keep track of samples
            read+=1
            if sti and not arr[sample_id_index] in sti:
                if read%5000==0:
                    print ("processed {} lines".format(read))
                continue
            sa = arr[sample_id_index]
            samples.add(sa)
            si = sample_info.get(sa,{"n_cells":0})
            sample_info[sa]=si
            si["n_cells"]+=1
            for f in conf["sample_fields"]:
                si[f] = arr[header_to_index[f]]
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

    #if pca fileexists add the data to the config file
    pca_file = conf.get("pca_file")
    if pca_file:
        if not os.path.exists(pca_file):
                print("cannot find specified pca file: {}".format(pca_file))
                sys.exit()
        else:
            pc = pandas.read_csv(conf.get("pca_file"),sep="\t")
            for r in pc.itertuples():
                sample_info[r[1]]["PCS"]=list(r[2:6])
    
    
   
    #zip all the images
    print ("zipping images")
    iti = conf.get("images_to_include")
   
    images_zip = os.path.join(conf_dir,"temp_images.zip")
    zip_file = zipfile.ZipFile(images_zip, 'w')
    for imgfolder in conf["image_folders"]:
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

            for sample in samples:        
                if im.startswith(sample):
                    zip_file.write(os.path.join(imgfolder,im),arcname=im, compress_type=zipfile.ZIP_DEFLATED)
                    break
    zip_file.close()

    #no longer needed
    del conf["data_file"]
    del conf["image_folders"]

    conf["sample_info"]=sample_info

    #save the config
    mlvconfig= os.path.join(conf_dir,"temp_config.json")
    o = open(mlvconfig,"w")
    o.write(json.dumps(conf,indent=2))
    o.close()
    

    files={
        "config":open(mlvconfig,"rb"),
        "data_file":open(t_f,"rb"),
        "images":open(images_zip,"rb")
    }

    data= {"token":token}
    if not conf.get("create_files_only"):
        print ("uploading files")
        r= requests.post(url+"/hyperion/create_hyperion_project",files=files,data=data)
        if not r.status_code==200:
            print ("there was a problem connecting to the server")

        if r.text == "ok":
            print("uploaded")
        elif r.text=="permission denied":
            print("Permission denied - do you have the correct token?")
        elif r.text=="error":
            print("There was an error on the server")

        if conf.get("delete_files_after_upload"):
            os.remove(mlvconfig)
            os.remove(t_f)
            os.remove(images_zip)
    else:
        print ("Files generated")




def make_pyramid_ome_tiff(folder,out_dir,tile_size=512):

    file_names={}
    for d in os.listdir(folder):
        for f in os.listdir(os.path.join(folder,d)):
            if f.endswith(".ome.tiff"):
                sample_id = f.replace(".ome.tiff","")
                
                in_file = os.path.join(folder,d,f)
                #work out dimensions for p_res 
                t=TiffFile(in_file)
                t_dim =t.pages[0].shape
                p_res =str(math.ceil(math.log2(max(t_dim)/tile_size)))
                r= ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5))
                out_name= f'{r}.ome.tiff'
                out_file= os.path.join(out_dir,out_name)
                params= [
                    "bfconvert",
                    "-tilex",str(tile_size),
                    "-tiley",str(tile_size),
                    "-pyramid-scale","2",
                    "-pyramid-resolutions",p_res,
                    "-compression","LZW",
                    in_file,
                    out_file
                ]
                subprocess.run(params)
                new_name= out_name.replace("tiff","png")
                new_file = out_file.replace(".tiff",".png")
                os.rename(out_file,new_file)
                file_names[sample_id]=new_name
                print(f'{sample_id}:{new_name}')
    o=open(os.path.join(out_dir,"names.json"),"w")
    o.write(json.dumps(file_names))
    o.close()





if __name__ == "__main__":
    conf= "config.yaml"
    if len(sys.argv)>1:
        conf= sys.argv[1]
    #main(conf)
    make_pyramid_ome_tiff(
        "/t1-data/project/BM_hyperion/shared/data/panel1/ometiff",
        "/t1-data/project/covidhyperion/sergeant/bm_ss/ometiff"
    )
