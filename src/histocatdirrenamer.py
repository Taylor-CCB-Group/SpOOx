#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:55:08 2021

@author: staylor
"""

import argparse
import sys
import os
import re

parser = argparse.ArgumentParser(description='''Rename folders labelled imc style to WIMM format
Input: Histocat directory and disease sample name e.g. DISEASE_SAMPLE_1
Output All subfolders renamed
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='indir',
                    help='histocat folder')
parser.add_argument('--pattern', dest='pattern', default='',
                    help='histocat folder')
parser.add_argument('--disease_sample', dest='disease_sample',default='',
                    help='stub for the disease_sample renaming e.g. covid_sample_1')
parser.add_argument('--roi', dest='roi', action='store_true', default=False,
                    help='only rename the roi of in the filename')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   

args = parser.parse_args()

def RenameDirs(rootdir):
    for it in os.scandir(rootdir):
        # check if the directory or symlink (for nextflow) exists and contains the pattern specified
        if (it.is_dir() or it.is_symlink()) and (args.pattern in it.name) or (args.roi == True):
            x =  re.search(r"_s0_a(\d+)_ac$", it.path)
            roi = x.group(1)
            print("ROI = "+ roi)
            print("it.path = "+ it.path)
            fullPath = os.path.abspath(it.path)
            basePath = os.path.abspath(it.path + "/..")
            if (args.roi == True):
                x = re.search(r"(\w+_SAMPLE_\d+)", it.path)
                args.disease_sample = x.group(1)
            newPath = basePath + "/" + args.disease_sample + "_ROI_"+roi
            print("Renaming " + fullPath + " to " + newPath)
            #os.rename(fullPath, newPath)
            os.system("mv " + '"' + fullPath + '"' + " " + newPath)
	    
RenameDirs(args.indir)
