#!/usr/bin/env python
# coding: utf-8

DESC = ("Copies the black and white mask images to zegami directory so it can be visualised. Should be run after deepcell but before Zegami collection upload.")
import os
import sys
import argparse
import warnings
import glob
import logging

logger=logging.getLogger(__name__)
# set logger to debug
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


parser = argparse.ArgumentParser(description='''Copy deeper cell images to png directory
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', 
                    help='input directory that will contain images from segmentation')
parser.add_argument('--outdir',
                    help='output image directory for pngs')
    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   


args = parser.parse_args()

inDir = args.indir
outDir = args.outdir
bwImageOut = "bandwmask.png"

# find all the files recursively in directory called bandwmask.png
#bwImageList = glob.glob(inDir + '/**/' + bwImageOut, recursive=True)
cellSegmentationList = glob.glob(inDir + '/**/' + bwImageOut, recursive=True)
# rename the files to the directory name and append "_Cell_Mask.png" and copy to outDir
for cellSegmentation in cellSegmentationList:
    cellSegmentationDir = os.path.dirname(cellSegmentation)
    cellSegmentationDirName = os.path.basename(cellSegmentationDir)
    cellSegmentationOut = cellSegmentationDirName + "_Cell_Mask.png"
    cellSegmentationOut = os.path.join(outDir,cellSegmentationOut)
    cmd = "cp " + cellSegmentation + " " + cellSegmentationOut
    logger.info("copying " + cellSegmentation + " to " + cellSegmentationOut + "...")
    os.system(cmd)

