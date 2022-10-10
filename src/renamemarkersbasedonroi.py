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


parser = argparse.ArgumentParser(description='''Using a specified ROI as reference, will compare to  all other ROIs to see if has the same the antibody / metal ion combination.''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='indir',
                    help='reference directory (ROI) to base renaming all other directories on.')
parser.add_argument('--histocatdir', dest='histocatdir',
                    help='output directory of histocat files.')
parser.add_argument('--rename', dest='rename', action='store_true',
                    help='rename all the files based on the reference naming scheme ')

args = parser.parse_args()
		    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   


def GetAllFiles(path):
# get all the files in subdirectories ending in .tiff    
    tiffList = []
    for root, dirs, files in os.walk(path):
        for file in files:
            #print(files)
            if(file.endswith(".tiff")):
                # typical file name looks like this:
                # CD11b_Sm149.tiff
                tiffList.append(file)

                # get parent dir of file
                parent = os.path.basename(os.path.dirname(os.path.join(root,file)))
                # construct dictionary with parent as key and marker as value
    return tiffList


# exmple
# CD11c_Sm154.tiff
# dictionary key is CD11c and value is Sm154
# from the tiff file name get the marker and the channel
def GetMarkerChannelFromTiff(tiff):
    marker = re.sub("_.*tiff","",tiff)
    channel = re.sub(".*_","",tiff)
    channel = re.sub("\.tiff","",channel)
    return marker,channel

def GetRoiMarkersDictFromTiff(tiffList):
# create a dictionary with the marker as key and the channel as value
    roiMarkersDict = {}
    for tiff in tiffList:
        marker,channel = GetMarkerChannelFromTiff(tiff)
        roiMarkersDict[channel] = marker
    return roiMarkersDict

# rename all the files in the histocat directory based on the roiMarkersDict
def RenameFiles(path,roiMarkersDict):

    for root, dirs, files in os.walk(path):
        for file in files:
            if(file.endswith(".tiff")):
                marker,channel = GetMarkerChannelFromTiff(file)
                # get parent directory
                parent = os.path.basename(os.path.dirname(os.path.join(root,file)))
                #print(marker,channel)
                # if the channel is in the roiMarkersDict and if it is different than the marker
                # stored in the roiMarkersDict then rename the file
                # check if KeyError
                try:
                    if(roiMarkersDict[channel] != marker):
                        # rename the file
                        newFileName = roiMarkersDict[channel] + "_" + channel + ".tiff"
                        print(parent + " old file name "+ file + " will be converted to new file called ",newFileName)
                        if (args.rename):
                            os.rename(os.path.join(root,file),os.path.join(root,newFileName))
                except KeyError:
                    print("Channel not found in roiMarkersDict so removing file",file)
                    #remove the file
                    if (args.rename):
                        os.remove(os.path.join(root,file))

def main():
    tiffList = GetAllFiles(args.indir)
    roiMarkersDict = GetRoiMarkersDictFromTiff(tiffList)
    # get the parent of the parent of args.indir
    #parent = os.path.dirname(os.path.dirname(args.indir))
    #print("parent",parent)
    RenameFiles(args.histocatdir,roiMarkersDict)


main()

