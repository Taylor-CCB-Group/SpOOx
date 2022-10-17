# read in file names from directory

import os
import glob
import pandas as pd
import numpy as np
import re
import sys
import argparse
import logging
import skimage.io  



logger = logging.getLogger(__name__)

# set logging verbosity to 5
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Do a QC check on histocat directoty and remove poor quality images or unnecessary images')
parser.add_argument('-i', '--indir', help='input directory', required=True)
parser.add_argument('-o', '--outdir', help='output directory', required=True)
parser.add_argument('--minImageHeight',
                       type=int,
                       default=25,
                       help='overrides the default fo the minimum image height (in pixels) to be considered for analysis')
parser.add_argument('--minImageWidth',
                       type=int,
                       default=25,
                       help='overrides the default fo the minimum image width (in pixels) to be considered for analysis')
# filter ou images whoe marker begines witha number argument
parser.add_argument('--filterNonMarkers', 
                        action='store_true', 
                        default= True, 
                        help='Remove images that are just used to calibrate the machine' )
args = parser.parse_args()

# find all tiff files in directory that have smaller than 25x25 pixels
files = glob.glob(args.indir + '/*/*.tiff')
if len(files) == 0:
    logger.error("No tiff files found in input directory")
    sys.exit()



# create directories based on file names
for file in files:
    # if file doesn't exist skip
    if not os.path.exists(file):
        continue

    # get parent dir
    parentDir = os.path.basename(os.path.dirname(file))
    markerName = os.path.basename(file).split('_')[0]

    logging.info("Checking " + file + " and " + markerName)
    # get image size
    im = skimage.io.imread(file)
    height, width = im.shape
    # check if image is too small
    # if image is too small, move to bad images directory
    if (width <= args.minImageWidth) or (height <= args.minImageHeight):
        #move parent directory to bad images directory
        #logging.warn(args.indir + '/' + parentDir, args.outdir + '/' + parentDir)
        logging.info("Images too small: " + parentDir + " moving to " + args.outdir)
        os.rename(args.indir + '/' + parentDir, args.outdir + '/' + parentDir)

for file in files:
    # if basename image name starts with a number, move to bad images directory
    if re.match(r'^\d', markerName):
        markerFile = file
        moveToDirectory = args.outdir + '/' + parentDir

        logging.info("Image name " + markerFile + " starts with a number; moving to " + moveToDirectory + " in " + args.outdir)
        if not os.path.exists(moveToDirectory):
            os.makedirs(moveToDirectory)
        os.rename(file, moveToDirectory + '/' + os.path.basename(file))


