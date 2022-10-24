# read in file names from directory

import os
import glob
from tabnanny import filename_only
import pandas as pd
import numpy as np
import re
import sys
import argparse
import logging
import skimage.io 
from shutil import move


logger = logging.getLogger(__name__)

# set logging verbosity to 5
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Do a QC check on histocat directoty and remove poor quality images or non biological markers or excude list of rois')
parser.add_argument('-i', '--indir', help='input directory', default="histocat")
parser.add_argument('-o', '--outdir', help='output directory', default='badimages')
parser.add_argument('-r', '--roi', nargs='+', help='ROI to remove', required=False)
parser.add_argument('-x','--minImageHeight', required=False,
                       type=int,
                       help='Minimum image height (in pixels) to be considered for analysis')
parser.add_argument('-y','--minImageWidth', required=False,
                       type=int,
                       help='Minimum image width (in pixels) to be considered for analysis')
# filter out images whose marker begins with a number argument
parser.add_argument('-n','--filterNonMarkers', 
                        action='store_true', 
                        default= False, 
                        help='Remove images that are just used to calibrate the machine' )
args = parser.parse_args()

# function to compare two strings are the same even if characters are in random order
def CompareStrings(str1, str2):
    a = sorted(str1)
    b = sorted(str2)
    if(sorted(str1)== sorted(str2)): 
        return True
    else: 
        return False


# exit if no arguments are given
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

# find all tiff files in directory that have smaller than 25x25 pixels (default)
files = glob.glob(args.indir + '/*/*.tiff')
if len(files) == 0:
    logger.error("No tiff files found in input directory")
    sys.exit()

# if outdir doesn't exist make it
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

if (args.roi):
# check if list of rois to remove are files that exists
    for roi in args.roi:
        if not os.path.exists(args.indir +'/'+ roi):
            logging.error("ROI " + roi + " does not exist in "+ args.indir)
        else:
            logging.info("Moving to bad images directory: " + args.indir +'/'+ roi)
            #force move to bad images directory
            os.rename(args.indir +'/'+ roi, args.outdir + '/' + roi)


# check minium size of minwidth and minheight arguments are given
if args.minImageWidth or args.minImageHeight:
    for file in files:
        # if file doesn't exist skip
        if not os.path.exists(file):
            continue

        # get parent dir
        parentDir = os.path.basename(os.path.dirname(file))
        markerName = os.path.basename(file).split('_')[0]
        logging.info("Checking " + file + " for size")
        # get image size
        im = skimage.io.imread(file)
        height, width = im.shape
        # if image is too small, move to bad images directory
        if (args.minImageWidth and width <= args.minImageWidth) or (args.minImageHeight and height <= args.minImageHeight):
            #move parent directory to bad images directory
            #logging.warn(args.indir + '/' + parentDir, args.outdir + '/' + parentDir)
            logging.info("Images too small: " + parentDir + " moving to " + args.outdir)
            os.rename(args.indir + '/' + parentDir, args.outdir + '/' + parentDir)

# get rid of non biological markers, these are the ones that have the same marker name as antibody name
if (args.filterNonMarkers):
    for file in files:
        logging.info("Checking " + file + " for non biological markers")
        parentDir = os.path.basename(os.path.dirname(file))
        filename_only = os.path.splitext(os.path.basename(file))[0]      
        markerName = filename_only.split('_')[0]
        metalName=filename_only.split('_')[1]
        logging.debug(metalName  + " " + markerName)
        moveToDirectory = args.outdir + '/' + parentDir
        #if (markerName == metalName):
        if CompareStrings(markerName, metalName):
            # make directory if it doesn't exist
            if not os.path.exists(moveToDirectory):
                os.makedirs(moveToDirectory)
            logging.info("Images are not biological markers ("+markerName+" == "+metalName+")" + parentDir + " moving to " + args.outdir)    
            if not os.path.exists(moveToDirectory):
                os.makedirs(moveToDirectory)
            os.rename(file, moveToDirectory + '/' + os.path.basename(file))


