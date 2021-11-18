#!/usr/bin/env python
# coding: utf-8

DESC = ("flatzegami2stackzegami - convert stack of flat Zegami images into a Zegami stacked collection")
from ast import literal_eval
import argparse
import os
from glob import glob
from shutil import copy
import re
import readconfig

#
#input_dir = "/t1-data/user/staylor/lho/hyperion/data/cun1/ImcSegmentationPipeline/Data/renamed/zegami/img"
#output_dir = "/t1-data/user/staylor/lho/hyperion/data/cun1/ImcSegmentationPipeline/Data/renamed/zegami/stacked"
#yaml = "/t1-data/user/staylor/lho/hyperion/data/cun1/ImcSegmentationPipeline/Data/renamed/zegami/stacked/stacked.yaml"
#zegami_out = "/t1-data/user/staylor/lho/hyperion/data/cun1/ImcSegmentationPipeline/Data/renamed/zegami/stacked/zegami.tsv"

def IsValidFile(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        #return open(arg, 'r')  # return an open file handle
        return arg

def File2List(file):
    print("file=",file)
    fileobj=open(file,"r")
    lines=[]
    for line in fileobj:
        lines.append(line.strip())
        print(line)
    return lines


parser = argparse.ArgumentParser(
        prog = "",
        description = DESC
)
parser.add_argument(
        '-i', '--input_dir', 
        help = 'Input directory of images e.g. named PNGs/JPGs',
        required = True
) 
parser.add_argument(
        '-o', '--output_dir', 
        help = 'directory to put images ready to be uploaded as a zegami stacked collection',
        required = True
)
parser.add_argument(
        '-y', '--yaml', 
        help = 'Output YAML file for Zegami upload',
        required = True
)
parser.add_argument(
        '-t', '--zegami_tsv', 
        help = 'Output TSV file for Zegami upload',
        required = True
)

parser.add_argument(
        '-n', '--name', 
        help = 'Name of Zegami collection',
        required = True
)

parser.add_argument(
        '-d', '--description', 
        help = 'Description of Zegami collection',
        required = True        
)

parser.add_argument('-m', '--markers',
        metavar='FILE', type=lambda x: IsValidFile(parser, x),
        help='file containing list of tiff files defining cytoplasm',
        required = False
)



args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

y = open(args.yaml,"w");
z = open(args.zegami_tsv,"w")

z.write("ImageName\tSample Name\tDisease\tSample\tROI\n")

a=("""name: {}
description: {}
dataset_type: file
enable_clustering: false
file_config:
  path: {}
collection_version: 2
image_sources:\n""").format(args.name,args.description,args.zegami_tsv)
        
y.write(a)

#markers= File2List(args.markers)
markers = readconfig.GetMarkerList(args.markers,"clustering")
print("MARKERS",markers)
#add deepcell to the list of markers
# todo check exists!
markers.append("deepcell_mask")

# goes though all the markers and identifies the name in zegami 'flat' image directory

rois = set()

for marker_name in markers:
    print(marker_name)
    # make the output directories
    marker_output_dir = args.output_dir + "/" + marker_name
    #print(marker_output_dir)
    if not os.path.exists(marker_output_dir):
        os.mkdir(marker_output_dir)
    glob_input = args.input_dir+'/*'+ marker_name +'*'
    
    print(">>>>" + glob_input + "<<<<<");
    
    #write yaml file
    y.write("  - paths:\n")
    y.write("     - "+marker_output_dir+"\n")
    y.write("    source_name: "+ marker_name+"\n")
    y.write("    dataset_column: ImageName\n")
    y.write("    imageset_type: file\n")
    
    for file in glob(glob_input):
        #print(">>>" + file + "<<<")
        x =  re.search(r"(\w+_SAMPLE_\d+_ROI_\d+)", file)
        new_file_name = x.group(1) + ".png";
        sample_name = x.group(1)
        
        x = re.search(r"(\w+)_SAMPLE_(\d+)_ROI_(\d+)", file)
        disease = x.group(1)
        sample = x.group(2)
        roi = x.group(3)
      
        #print(">>>" + new_file_name + "<<<")
        copy(file,marker_output_dir+"/" + new_file_name)
        # add to set to avoid duplications
        rois.add(new_file_name + "\t" + sample_name + "\t" + disease + "\t" + sample + "\t" + roi + "\n")


#convert set to string
# from https://stackoverflow.com/questions/17528374/python-convert-set-to-string-and-vice-versa
str=''.join(list(map(str, rois)))
z.write(str)

y.close()
z.close()
