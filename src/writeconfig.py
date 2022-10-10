#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import argparse
import os
import sys
import pandas as pd
import re

# Command line script to:
# INPUTS:
#       Histocat directory
#       Config file (this should into the config directory)


logger = logging.getLogger(__name__)

columns = ["marker_name","nucleus","cytoplasm","clustering"]

parser = argparse.ArgumentParser(description='''Outputs a table (markers_panel.tsv) for specificing the makers to be used in various parts of the pipeline and an overview showing all the markers across all ROIs in the project (marker_overview.tsv). This is useful to identify misspellings and other errors in the marker names.
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='indir',
                    help='input directory of histocat files.')
parser.add_argument('--outfile', dest='outfile', default="markers_panel.tsv",
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
                # typical file name looks like this:
                # CD11b_Sm149.tiff
                marker = re.sub("_.*tiff","",file)
                tiffList.append(marker)

                # get parent dir of file
                parent = os.path.basename(os.path.dirname(os.path.join(root,file)))
                # construct dictionary with parent as key and marker as value
    return tiffList

def GetMarkersPerRoi(path):
# get a list of top level directories in histocat dir which correspond to ROIs
    roiList =[]
    # get a list of durectories in the histocat dir
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            roiList.append(dir)
    logger.error(roiList)
    #print(roiList)

    # put list of files in each directory into a dictionary
    markerDict = {}
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            markerDict[dir] = GetAllFiles(os.path.join(path,dir))  
    #logger.error(markerDict)
    return markerDict

def WriteDebugFile(markerDict,outfile):
# write a config file with the markers in each roi
    with open(outfile, 'w') as f:
        for roi in markerDict:
            f.write(roi + "," + ",".join(sorted(markerDict[roi])) + "\n")


# check there are files ending in .tiff
tiffList = GetAllFiles(args.indir)
if len(tiffList) == 0:
    logger.error("No tiff files found in input directory")
    sys.exit()


markerSet = set(GetAllFiles(args.indir))
markerList = list((sorted(markerSet))) 
# make empty dataframe
df = pd.DataFrame(columns = columns, dtype=object)
df.style.hide_index()
df['marker_name'] = pd.Series(markerList)
df.fillna(0, inplace=True)
logging.info("Marker config table '" + args.outfile + "' generated.")
#logging.info(df)
df.to_csv(args.outfile,sep="\t",index=False)

# write a debug file to look for naming problems and different markers in each ROI
WriteDebugFile(GetMarkersPerRoi(args.indir),"marker_overview.csv")
logging.info("Marker table 'marker_overview.tsv' generated.")




