#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
from tifffile import TiffFile


parser = argparse.ArgumentParser(description='''Generate a metadata file describing the data analysed. Output will be tab-delimited: sample_id, sample_name, condition, ROI, path.
''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir', dest='indir', default = "deepcell",
                    help='the location of files used to generate this metadata')		
parser.add_argument('--outfile', dest='outfile', default = "metadata.tsv",
                    help='output file')		

args = parser.parse_args()

indir=args.indir
outfile=args.outfile

cwd = os.getcwd()

filedata = "sample_id\tsample_name\tcondition\tROI\tpath\tregion_width\tregion_height\n"
for d in os.listdir(indir):
    pattern = '_ROI_'
    if (re.search(pattern, d)):
        (sample_name, ROI) = re.split(pattern, d)
        condition= d.split("_SAMPLE_")[0]
        full_path = cwd + '/' + 'signalextraction' + '/' + d
        metadata = f'{d}\t {sample_name}\t{condition}\tROI_{ROI}\t{full_path}'
        #try an get the deep cell tiff
        tiff = os.path.join(indir,d,"deepcell.tif")
        if os.path.exists(tiff):
            #extract width and height
            t=TiffFile(tiff)
            t_dim =t.pages[0].shape
            metadata +=f'\t{t_dim[0]}\t{t_dim[1]}\n'
        else:
            metadata+="\tNA\tNA\n"


        filedata += metadata


        
metadata_file = cwd + '/' + outfile
f = open(metadata_file,'w')  
f.write(filedata)
f.close()
