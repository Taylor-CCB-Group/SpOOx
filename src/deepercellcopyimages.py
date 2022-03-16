#!/usr/bin/env python
# coding: utf-8

DESC = ("Copies the black and white mask images to zegami directory so it can be visualised. Should be run after deepcell but before Zegami collection upload.")

import cv2
from os import listdir, mkdir, read
from os.path import isfile, join
import os
import re
import fnmatch
import sys
from skimage.io import imread
from skimage.io import imsave
from skimage.segmentation import find_boundaries
from skimage import img_as_ubyte
from matplotlib import pyplot as plt
import multipagetiff as tiff
import glob
import tifffile
import numpy as np
import argparse
import warnings

# parses SpOOx config files
import readconfig
warnings.filterwarnings("ignore")

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
    return lines

def TiffNorm(src_image,dest_image):
    img = cv2.imread(src_image, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    normed = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    contrast = 20
    brightness = 0
    out = cv2.addWeighted(normed, contrast, normed, 0, brightness)
    cv2.imwrite(dest_image,out)


def Tiff2Stack(fileName, imageList, dirName):
    with tifffile.TiffWriter(fileName) as stack:
        for imgFilename in imageList:
            #fullPath = dirName + "/" + imgFilename
            # check if the file contains the antibody name
            expression = dirName+"/"+imgFilename+"_"+"*";
            fullPath = glob.glob(dirName+"/"+imgFilename+"_"+"*.tiff")
            #print(dirName+"/"+imgFilename+"_"+"*.tiff")
            #quit()
            try:
                fullPath = glob.glob(dirName+"/"+imgFilename+"_"+"*.tiff")[0]
                print("Ab fullpath = ",fullPath)
                stack.save(tifffile.imread(fullPath), photometric='minisblack',contiguous=True)
            except:
                print("Can't find " + imgFilename)


def BlackAndWhite(tiffFileName,blackAndWhiteOut):
# convert label matrix tif file to b/w png 
    label = imread(tiffFileName)
    boundaries = find_boundaries(label)
    toPlot = label>0
    toPlot = toPlot & ~boundaries
    toPlot = img_as_ubyte(toPlot)
    imsave(blackAndWhiteOut,toPlot)


parser = argparse.ArgumentParser(description='''Generate label matrix mask from tiffs of nuclear and cytoplasm markers.
Input is both a nuc and cyt file of histocat tiff filenames of the markers names.
Output is label matrix file in tiff format.
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--indir', dest='imageDir',
                    help='output directory that will contain images from segmentation')
parser.add_argument('--zegamidir', dest='zegimageDir',
                    help='output image directory for zegami')
    
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   


args = parser.parse_args()
    
# includes miscellaneous output files which can be useful for debugging
imageNuc = os.path.join(args.imageDir,"nuc.tif")
imageCyt = os.path.join(args.imageDir,"cyt.tif")
imageNucFlat = os.path.join(args.imageDir,'nuc_flat.tif')
imageCytFlat = os.path.join(args.imageDir,'cyt_flat.tif')
imageNucFlatNorm = os.path.join(args.imageDir,'nuc_flat_norm.tif')
imageCytFlatNorm = os.path.join(args.imageDir,'cyto_flat_norm.tif')
imageOut = os.path.join(args.imageDir,"deepcell.tif")
bwImageOut = os.path.join(args.imageDir,"bandwmask.png")
zegamiDir = args.zegimageDir

#split pathname
elements = []
elements = args.imageDir.split("/")
condition=elements[-4]
sample=elements[-3]
roi=elements[-2]
newName=condition+"_"+sample+"_"+roi+".png"

if (os.path.isfile(bwImageOut)):
    zegami_deepcell_mask_dir = os.path.join(zegamiDir,args.imageDir.split("/")[-1])
    zegami_deepcell_mask_file = os.path.join(zegami_deepcell_mask_dir,newName) 
    os.makedirs(zegami_deepcell_mask_dir,exist_ok=True)
    print("copying",bwImageOut,"to",zegami_deepcell_mask_file)
    os.system("cp "+bwImageOut+" "+zegami_deepcell_mask_file)



