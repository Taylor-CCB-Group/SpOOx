#!/usr/bin/env python
# coding: utf-8

import cv2
from os import listdir
from os.path import isfile, join
import os
import argparse
import re
import sys

src_image = '/t1-data/user/staylor/lho/hyperion/data/cun1/ImcSegmentationPipeline/Data/output/histocat/CUN6USS3COVIDHyperion2runA_s0_p4_r1_a1_ac/ACE2_Ho165.tiff'
dest_image = '/t1-data/user/staylor/lho/hyperion/data/cun1/tmp/test.png'

def tiff2png(src_image,dest_image):
    img = cv2.imread(src_image, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    normed = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    contrast = 10
    brightness = 0
    out = cv2.addWeighted(normed, contrast, normed, 0, brightness)
    cv2.imwrite(dest_image,out)


def roi2images(src_dir,dest_dir,zegami_tsv):
    image_dest_dir = dest_dir + "/img"
   
    #make directory 
    if not os.path.exists(image_dest_dir):
        os.mkdir(image_dest_dir)
    
    onlyfiles = [f for f in listdir(src_dir) if isfile(join(src_dir, f))]
    
    
    for file in onlyfiles:    
        print(file)
        src_image = src_dir + "/" + file
        basename = os.path.splitext(os.path.basename(src_image))[0];
        dirname = os.path.dirname(src_image)
        sample_name = os.path.basename(dirname);
        dest_image = sample_name + "_" + basename + "." + output_format
        full_path_dest_image = image_dest_dir + '/' + dest_image
        print("dest_image ", dest_image)
        print("sample_name ", sample_name)
        
        #full condition / sample / id
        x =  re.search(r"^(\w+_SAMPLE_\d+_ROI_\d+)", sample_name)
        full_sample_name = x.group(1)
        print("full_sample_name = "+ full_sample_name)

        #disease name
        x =  re.search(r"^(\S+?)_", sample_name)
        disease = x.group(1);
        print("disease = "+ disease)
        
        # get sample id
        x =  re.search(r"_SAMPLE_(\d+)_", sample_name)
        sample_id = x.group(1);
        print("sample_id = "+ sample_id)
                
        # get roi from directory name
        x = re.search(r"_ROI_(\d+)",dirname)
        roi = x.group(1);
        print("roi = " + roi)

        replicate = disease+"_"+sample_id
        print("replicate = " + replicate)

        
        # protein and metal name from image name
        x = re.search(r"^(\w+?)_(\w+)", basename)
        print("basename = "+basename)
    
        tab = "\t"
    
        if x:
            print("found")
            protein = x.group(1)
            metal = x.group(2)
            print(protein, metal)
            
            # write the line of to the tsv
            fields = (dest_image,full_sample_name,sample_id,replicate,roi,protein, metal)
            print("Writing to",fields,"to",output_tsv)
            f.write(tab.join(fields))
            f.write("\n")
        else:
            print("Missing fields (protein/metal) in " + basename)
     
        print(">>>" + full_path_dest_image)  
        tiff2png(src_image,full_path_dest_image)
    

def MakeYaml(yamlFile, name, description,datasetPath, imagesetPath):
    f = open(yamlFile,'w')
    yamlText = f'''# zeg create collections --project XY6PsWre --config zegami.yaml
# The name of the collection
name: {name}
description: {description}
# The type of data set. For now this needs to be set to 'file'. (optional)
dataset_type: file
# Config for the file data set type
imageset_type: file
# Config for the file image set type
file_config:
# Whether to recursively scan any directories. (optional)
    recursive: False
# If provided, the mime-type to use when uploading images. (optional)
    mime_type: image/png
# Path to the dataset file. (optional)
    path: /{datasetPath}
# A collection of paths to image files. Paths can be to both images and directories
    paths:
        - {imagesetPath}
# Name of the column in the dataset that contains the image name. (optional)
dataset_column: ImageFile'''
    f.write(yamlText)

project_dir = '/stopgap/hyperion/msalio/hodgkin/FOXP3/histocat'
dest_dir = '/stopgap/hyperion/msalio/hodgkin/FOXP3/zegami'

parser = argparse.ArgumentParser(description='''Generate an image normalised zegami collection of ROIs in one view.
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='indir',
                    help='input directory of histocat files.')
parser.add_argument('--outdir', dest='outdir',
                    help='output directory for zegami files')
parser.add_argument('--yaml', dest='yaml',
                    help='output yaml file to upload to zegami via cli')
parser.add_argument('--name', dest='name',
                    help='name of zegami collection')

args = parser.parse_args()
		    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   

project_dir = os.path.abspath(args.indir)
dest_dir = os.path.abspath(args.outdir)

output_format = "png"

output_tsv = dest_dir + '/' + 'zegami.tsv'
f = open(output_tsv,'w')
f.write("ImageFile\tfull_sample_name\tsample_name\treplicate\troi\tprotein\tmetal\n")

onlydirs = [d for d in listdir(project_dir) if os.path.isdir(join(project_dir, d))]
for dir in onlydirs:   
    sample_dir = project_dir + "/" + dir 
    print("Processing "+sample_dir+"...")
    roi2images(sample_dir,dest_dir,output_tsv)

MakeYaml(args.yaml, args.name,"", output_tsv, dest_dir + "/img")



