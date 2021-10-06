#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:43:52 2021

@author: staylor
"""
import sys
import skimage.io
import numpy as np
import argparse
import re
from skimage.segmentation import find_boundaries
from skimage.io import imread
from skimage.io import imsave
from skimage import img_as_ubyte

parser = argparse.ArgumentParser(description='Generate black and white image with boundaries image from label matrix file.')
parser.add_argument('--infile', dest='imageIn',
                    help='input labelled matrix data (usually 16 bit tiff) from segmentation')
parser.add_argument('--outfile', dest='imageOut',
                    help='output file in black and white with boundaries')
parser.add_argument('--writepng', action='store_true',
                    help='output pngfile based on the name of the input file')



if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		
args = parser.parse_args()

file = args.imageIn;

if (args.writepng):
	prefix = re.sub('\.tif$','',file)
	print(prefix)

label = skimage.io.imread(file)
boundaries = find_boundaries(label)
toPlot = label>0
toPlot = toPlot & ~boundaries
toPlot = img_as_ubyte(toPlot)
imsave(args.imageOut,toPlot)
