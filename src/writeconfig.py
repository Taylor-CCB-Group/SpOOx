#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import sys
import pandas as pd
import re

# Command line script to:
# INPUTS:
#       Histocat directory
#       Config file (this should into the config directory)

columns = ["marker_name","nucleus","cytoplasm","clustering"]

parser = argparse.ArgumentParser(description='''Generate a table used for specificing the makers to be used in various parts of the pipeline
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='indir',
                    help='input directory of histocat files.')
parser.add_argument('--outfile', dest='outfile',
                    help='output config file')
parser.add_argument('--columns', dest='columns',
                    help='names of columns that will be output as part of config in addition to makers')

args = parser.parse_args()
		    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   


def GetAllFiles(path):
# get all the files in subdirectroes ending in .tiff    
    tiffList = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if(file.endswith(".tiff")):
                #print((os.path.join(root,file)))
                marker = re.sub("_.*tiff","",file)
                tiffList.append(marker)
    return tiffList

markerSet = set(GetAllFiles(args.indir))
markerList = list((sorted(markerSet))) 
# make empty dataframe
df = pd.DataFrame(columns = columns, dtype=object)
df.style.hide_index()
df['marker_name'] = pd.Series(markerList)
df.fillna(0, inplace=True)
print("Marker table looks like the following:")
print(df)
df.to_csv(args.outfile,sep="\t",index=False)


#print(set().union(GetAllFiles(args.indir)))


