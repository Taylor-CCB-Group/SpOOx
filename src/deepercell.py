#!/usr/bin/env python
# coding: utf-8

from email.policy import default
from PIL.Image import composite
import cv2
from os import listdir, read
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

def TiffNorm(src_image,dest_image,contrast,brightness):
    img = cv2.imread(src_image, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    normed = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    print("contrast=",contrast)
    print("brightness=",brightness)
    out = cv2.addWeighted(normed, contrast, normed, 0, brightness)
    cv2.imwrite(dest_image,out)


def Tiff2Stack(fileName, imageList, dirName):
    with tifffile.TiffWriter(fileName) as stack:
        for imgFilename in imageList:
            #fullPath = dirName + "/" + imgFilename
            # check if the file contains the antibody name
            expression = dirName+"/"+imgFilename+"_"+"*";
            fullPath = glob.glob(dirName+"/"+imgFilename+"_"+"*.tiff")
            print(dirName+"/"+imgFilename+"_"+"*.tiff")
            #quit()
            try:
                fullPath = glob.glob(dirName+"/"+imgFilename+"_"+"*.tiff")[0]
                print("Ab fullpath = ",fullPath)
                stack.save(tifffile.imread(fullPath), photometric='minisblack',contiguous=True)
            except:
                print("Can't find " + imgFilename)


def BlackAndWhite(tiffFileName,blackAndWhiteOut):
# convert label matrix tif file to b/w png 
    print("Converting ",tiffFileName,"to ",blackAndWhiteOut)
    # test if tiffFileName is a label matrix tif exists
    if os.path.isfile(tiffFileName):  
        label = imread(tiffFileName)
        boundaries = find_boundaries(label)
        toPlot = label>0
        toPlot = toPlot & ~boundaries
        toPlot = img_as_ubyte(toPlot)
        imsave(blackAndWhiteOut,toPlot)
    else:
        print("Can't find ",tiffFileName) 
   

def ZProject(tiffFileName,zProjectOut,channelList):
    # convert tiff stack to z-projected tiff stack, unless only one slice
    # then just copy the file
    # we need the channelList to know how many channels so can fix a weird feature
    # where axes get swapped
    # for example if there are <3 or >4 channels img.shape= (27, 1000, 1000) which is incorrect
    # since surely it should be (1000,1000,27) but if there are 3 or 4 channels it is correct
    # img.shape= (3, 1000, 1000) so we must make it consistent for mp.max
    # I do this by swapping the axes (see below)
    numChannels = len(channelList)
    print("Converting ",tiffFileName,"to flattened ",zProjectOut)
    # read in a tiff stack
    if os.path.isfile(tiffFileName):
        img = imread(tiffFileName,cv2.IMREAD_GRAYSCALE)
        h,w,c = img.shape
        print("h=",h,"w=",w,"c=",c)
        # if only one slice, just copy the file
        if img.shape[0] == 1:
            print("Only one slice, copying file")
            os.system("cp "+tiffFileName+" "+zProjectOut)
        else:
            # z-project the tiff stack
            # due to a strange quirk the height and channel are swapped when there is 3 or 4 channels
            # so need to swap them back
            if numChannels == 3 or numChannels == 4:
                img = np.moveaxis(img,-1,0)
            print("New shape=",img.shape)
            imgNew = np.max(img,axis=0).copy()
            imsave(zProjectOut,imgNew)

parser = argparse.ArgumentParser(description='''Generate label matrix mask from tiffs of nuclear and cytoplasm markers.
Input is both a nuc and cyt file of histocat tiff filenames of the markers names.
Output is label matrix file in tiff format.
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--markerfile', dest='markerFile',
                    help='config file containing list of markers with columns for nucleus and cytoplasm')
parser.add_argument('--cyt', dest='cytolist', metavar='FILE', type=lambda x: IsValidFile(parser, x),
                    help='file containing list of tiff files defining cytoplasm')
parser.add_argument('--nuc', dest='nuclist', metavar='FILE', type=lambda x: IsValidFile(parser, x),
                    help='file containing list of tiff files defining nucleus')
parser.add_argument('--outdir', dest='imageOutDir',
                    help='output directory that will contain images from segmentation')
parser.add_argument('--indir', dest='dirName',
                    help='directory of images (imctools histocat directory)')
parser.add_argument('--brightness',type=int, default =0,
                    help='brightness of the image')
parser.add_argument('--contrast',type=int, default=5,
                    help='contrast of the image')
# Note default is 1 which is hyperion specific
parser.add_argument('--imageMpp', dest='imageMpp', default = 1,
                    help='microns per pixel')		
# overwrite the image directory if true or false
parser.add_argument('--overwrite', dest='overwrite', default = False, action='store_true',
                    help='If a deepcell file exists overwrite it if true.')	


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   


args = parser.parse_args()
os.makedirs(args.imageOutDir,exist_ok=True)

#try:
#    os.mkdir(args.imageOutDir)
#except:
#    print("Can't make " + args.imageOutDir + e)
    
# includes miscellaneous output files which can be useful for debugging
imageNuc = os.path.join(args.imageOutDir,"nuc.tif")
imageCyt = os.path.join(args.imageOutDir,"cyt.tif")
imageNucFlat = os.path.join(args.imageOutDir,'nuc_flat.tif')
imageCytFlat = os.path.join(args.imageOutDir,'cyt_flat.tif')
imageNucFlatNorm = os.path.join(args.imageOutDir,'nuc_flat_norm.tif')
imageCytFlatNorm = os.path.join(args.imageOutDir,'cyto_flat_norm.tif')
imageNucFlatNormPNG = os.path.join(args.imageOutDir,'nuc_flat_norm.png')
imageCytFlatNormPNG = os.path.join(args.imageOutDir,'cyto_flat_norm.png')
imageOut = os.path.join(args.imageOutDir,"deepcell.tif")
mergedImage = os.path.join(args.imageOutDir,"mergenucandcyt.png")
bwImageOut = os.path.join(args.imageOutDir,"bandwmask.png")
contrast = args.contrast
brightness = args.brightness


print("*** Analysing",args.dirName," ***")


#check if deepcell.tif exists
if(os.path.isfile(imageOut) and args.overwrite == False):
    print(imageOut," already exists ")
    quit()


# get the channels from the config file and convert to a list
cytoList = readconfig.GetMarkerList(args.markerFile, "cytoplasm")
nucList = readconfig.GetMarkerList(args.markerFile,"nucleus")


# convert those images to a tiff stack
print("Assembling nuclear channel containing ", nucList)
Tiff2Stack(imageNuc,nucList,args.dirName)
print("Assembling cytoplasmic channel containing ", cytoList)
Tiff2Stack(imageCyt,cytoList,args.dirName)

print(len(nucList)," nuclear images")  
print(len(cytoList)," cytoplasmic images")

# Z project and flatten nuc to 1 image
ZProject(imageNuc,imageNucFlat,nucList)
ZProject(imageCyt,imageCytFlat,cytoList)

#enhance contrast on image
TiffNorm(imageNucFlat, imageNucFlatNorm, contrast,brightness)
TiffNorm(imageCytFlat, imageCytFlatNorm, contrast,brightness)

imNuc = imread(imageNucFlatNorm)
imCyt = imread(imageCytFlatNorm)
height=imNuc.shape[0]
width=imNuc.shape[1]
blankImage = np.zeros((height,width), np.uint16)

# write out png version of the image
cv2.imwrite(imageNucFlatNormPNG,imNuc)
cv2.imwrite(imageCytFlatNormPNG,imCyt)
#overlay the images
m = np.dstack((imNuc,imCyt,imNuc))
#m = np.dstack((imNuc,imCyt,mask))
cv2.imwrite(mergedImage,m)


# Combined together and expand to 4D
im = np.stack((imNuc, imCyt), axis=-1)
im = np.expand_dims(im,0)

# now import the deepcell stuff
from deepcell.applications import Mesmer

# Create the application
app = Mesmer()

# create the label matrix
labeled_image = app.predict(im, image_mpp=args.imageMpp)
imsave(imageOut, labeled_image[0,:,:,0])

# create b/w image if specified for visualisation
BlackAndWhite(imageOut,bwImageOut)


