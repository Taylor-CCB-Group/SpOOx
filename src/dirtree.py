#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 06:48:45 2021

@author: staylor
"""
import os
import re
import shutil
import argparse
import sys

#inDir = "/t1-data/project/covidhyperion/staylor/covid2/histocat"
#outDir = "/t1-data/project/covidhyperion/staylor/covid2/tree"

sampleSource = []
sampleIds = []
roiIds = []

def MakeDir(dir):
    try:
        print("Making ",dir)
        os.mkdir(dir)
    except:
        print("Directory "+ dir + " exists.")

def HistocatToTree(inDir,outDir):
# read through all the files in histocat directory 
# now we know which rois have been successfully extracted
    dirs = os.listdir(inDir)
    for file in dirs:
        x =  re.search(r"(\w+)_(SAMPLE_\d+)_(ROI_\d+)", file)
        sampleSource = x.group(1);
        sampleId = x.group(2)
        roiId = x.group(3)    
        print("sampleSource = "+ sampleSource)
        print("sampleId = "+ sampleId)
        print("roiId = "+ roiId)
        
        #make directory tree
        MakeDir(outDir)
        sampleSourceDir = os.path.join(outDir,sampleSource)
        MakeDir(sampleSourceDir)
        sampleDir = os.path.join(sampleSourceDir,sampleId)
        MakeDir(sampleDir)
        roiDir = os.path.join(sampleDir,roiId)
        MakeDir(roiDir)
        pathToHistocatROI=os.path.join(inDir,file)
        pathToHistocatTreeLocation=os.path.join(outDir,sampleSourceDir,sampleDir,roiDir,"histocat")

        #make a symlink to this folder from the histocat dir
        #os.symlink(pathToHistocatROI,pathToHistocatTreeLocation)
        #shutil.copytree(pathToHistocatROI,pathToHistocatTreeLocation)
        # only copy across histocate files if they have been update compared to destination
        #print("cp -R -u -p "+pathToHistocatROI+" "+pathToHistocatTreeLocation);
        os.system("cp -R -u -p "+pathToHistocatROI+" "+pathToHistocatTreeLocation)

def ParseName(file, outDir):
# parses the name of a file or directory and returns the relevant path in the tree
    file = os.path.basename(file)
    x =  re.search(r"(\w+)_(SAMPLE_\d+)_(ROI_\d+)", file)
    treePath = os.path.join(outDir,x.group(1),x.group(2),x.group(3))
    try:
        os.path.isdir(treePath)
        return(treePath)
    except:
        print(treePath, "does not exist so cannot write the file there.")



parser = argparser = argparse.ArgumentParser(description='''Generate a directory tree for spatial omics analysis
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='inDir',
                    help='input directory of histocat files.')
parser.add_argument('--outdir', dest='outDir',
                    help='output directory for analysis files')
args = parser.parse_args()
		    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)	

#print(ParseName("./abc/COVID_SAMPLE_11_ROI_1", outDir))
HistocatToTree(args.inDir,args.outDir)
          