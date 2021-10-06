#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jbull, staylor
"""
# SignalExtraction_CommandLine.py
import argparse
import os
from pathlib import Path
import sys
import skimage.io
import skimage.color
import skimage.measure
from skimage.transform import rescale
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from cv2 import bitwise_and
from skimage.segmentation import find_boundaries
import re

# ST library
import readconfig

# Command line script to:
# INPUTS:
#       Path to Label matrix
#       Path to folder of multiplex images - as numChannels images, size (h x w)
#       Path to file with channel names and stains

# Currently running as a script for testing purposes
runAsScript = False
statisticsTable = "cellData.tab"

def AnalysisNameFromPath(path):
    # Returns a string that uniquely identifies the roi based on the directory name
    # input:
    # /t1-data/project/covidhyperion/staylor/covid2/tree/COVID/SAMPLE_1/ROI_1/deepcell
    # output:
    # COVID_SAMPLE_1_ROI_1
    
    # normpath collapses redundant separators and uplevel references 
    (p0,analysis)=os.path.split(os.path.normpath(os.path.dirname(path)))
    (p1,roi)=os.path.split(p0)
    (p2,sample)=os.path.split(p1)
    (p3,source)=os.path.split(p2)
    name = '_'.join([source,sample,roi])
    return name
    

if not runAsScript:
    # Create the parser
    parser = argparse.ArgumentParser(description='Extract marker intensities on a per-cell basis, based on existing cell segmentation data.')

    # Add the arguments
    parser.add_argument('pathToLabelMatrix',
                       type=str,
                       help='the path to the label matrix (as .tiff)')
    parser.add_argument('pathToIntensityDataFolder',
                       type=str,
                       help='the path to the folder of intensity data (currently histocat)')
    parser.add_argument('pathToRgbImage',
                       type=str,
                       help='the path to the rgb image to use for Zegami collections. To make a new rgb image based on histocat, specify "MAKE_NEW" and use the --set_rgb flag')
    parser.add_argument('pathToSaveAnalysis',
                       type=str,
                       help='file path to folder to save the analysis under')
    parser.add_argument('--verbose',
                        help='display progress messages',
                        action="store_true")
    parser.add_argument('--saveSegmentationImage',
                        help='save a png image of the label matrix used',
                        action="store_true")
    parser.add_argument('--makeCellImages',
                        help='generate png images of each individual cell for use in Zegami',
                        action="store_true")
    parser.add_argument('--cellSegmentationImagesExtension',
                        help='file extension used to save cell images (default .png)',
                        type=str,
                        choices=['.png','.jpg','.tif'],
                        default='.png')
    parser.add_argument('--doArcsinh',
                        help='apply arcsinh transform to intensity values before outputting cell data. If true, both transformed and non-transformed values will be returned',
                        action="store_true",
                        default=False)
    parser.add_argument('--arcsinhCofactor',
                        help='cofactor in the arcsinh transform (result = arcsinh(input / cofactor) ). Default 5.',
                        type=float,
                        default=5)
    parser.add_argument('--set_rgb',
                        help='used when pathToRgbImage equals MAKE_NEW. Specify as a single string of 3 histocat names separated by commas, e.g. "DNA1_Ir191,PanK_Nd148,Ki67_Er168". These will be assigned red, green, and blue respectively. ',
                        type=str)
    parser.add_argument('--backgroundColour',
                    help='specify background colour',
                    type=str,
                    choices=['white','black','grey'],
                    default='white')
    parser.add_argument('--distribution_crop_percentile',
                        help='After transforming the mean intensity data, this argument allows any data above the \'distribution_crop_percentile\' to be capped at that percentile, effectively trimming the data above that threshold to that level',
                        type=float,
                        default=0.99)
    parser.add_argument('--normalise_intensities',
                        help='Rescale the intensity distributions to the range [0,1]. If applicable, rescaling happens after taking an arcsinh transform or capping intensity values to a fixed percentile',
                        action="store_true",
                        default=True)
    parser.add_argument('--analysisName',
                       type=str,
                       help='overrides the default the name of this analysis which is derived from the directory name')

    # Execute the parse_args() method
    args = parser.parse_args()
    # print(args)
    #st if analysis name not specified autogenerate from the file path of the label matrix path
    #if (args.analysisName):
    #    analysisName = args.analysisName
    #else:
    
    analysisName = AnalysisNameFromPath(args.pathToLabelMatrix)
    print(args.pathToLabelMatrix)
    print("Name of analysis  : " + analysisName)
    
    pathToLabelMatrix = args.pathToLabelMatrix
    pathToIntensityDataFolder = args.pathToIntensityDataFolder
    pathToRgbImage = args.pathToRgbImage
    pathToSaveAnalysis = args.pathToSaveAnalysis + '/' # to ensure its a directory
    verbose = args.verbose
    saveSegmentationImage = args.saveSegmentationImage
    makeCellImages = args.makeCellImages
    cellSegmentationImagesExtension = args.cellSegmentationImagesExtension
    backgroundColourStr = args.backgroundColour
    rgb_string = args.set_rgb
 
    doArcsinh = args.doArcsinh
    arcsinhCofactor = args.arcsinhCofactor
    distribution_crop_percentile = args.distribution_crop_percentile
    normalise_intensities = args.normalise_intensities
    
    
    if backgroundColourStr == 'white':
        backgroundColour = [1,1,1]
    elif backgroundColourStr == 'black':
        backgroundColour = [0,0,0]
    elif backgroundColourStr == 'grey':
        backgroundColour = [0.5,0.5,0.5]         
        
    if not os.path.isfile(pathToLabelMatrix):
        print('The label matrix specified does not exist')
        sys.exit()

    if not os.path.isdir(pathToIntensityDataFolder):
        print('The intensity data folder specified does not exist')
        sys.exit()
        
    if pathToRgbImage != 'MAKE_NEW':
        if not os.path.isfile(pathToRgbImage):
            print('The rgb image specified does not exist')
            sys.exit()
        
else:
    # For testing purposes
    analysisName = 'COVIDpanel2_SAMPLE_1_ROI_1'
    pathToLabelMatrix = '/stopgap/hyperion/lho/covid2/deepcell/'+analysisName+'.tif'
    pathToIntensityDataFolder = '/stopgap/hyperion/lho/covid2/histocat/'+analysisName+'/'
    pathToRgbImage = 'MAKE_NEW'#'/stopgap/hyperion/lho/covid/covid/ilastikprobabilities/COVID_Sample_1_ROI_2_Probabilities.tiff'
    pathToSaveAnalysis = '/stopgap/hyperion/lho/covid2/signal_extraction/'+analysisName+'_TEMP_TEST/'
    verbose = True
    saveSegmentationImage = True
    makeCellImages = True
    cellSegmentationImagesExtension = '.png'
    backgroundColour = [1,1,1]
    doArcsinh = True
    arcsinhCofactor = 5
    rgb_string = 'DNA1_Ir191,Collagen1_Tm169,CD68_Tb159'
    
    distribution_crop_percentile = 0.99
    normalise_intensities = True

    
path = Path(pathToSaveAnalysis)
path.mkdir(parents=True, exist_ok=True)

# if not os.path.exists(pathToSaveAnalysis):
#     os.mkdirs(pathToSaveAnalysis)
if makeCellImages:
    path = Path(pathToSaveAnalysis + 'SegmentationImages/')
    path.mkdir(parents=True, exist_ok=True)
    # if not os.path.exists(pathToSaveAnalysis + 'SegmentationImages/'):
    #     os.mkdirs(pathToSaveAnalysis + 'SegmentationImages/')
    
    

# Set up colormap for the label matrix
color_bins = np.linspace(0, 1, 1000)
np.random.shuffle(color_bins)
norm_bins = mpl.colors.Normalize(vmin=0, vmax=1)
m = plt.cm.ScalarMappable(norm=norm_bins, cmap=plt.get_cmap('jet')) # There aren't many occasions when jet is an appropriate colormap to use, so let's make the most of it
colormap = m.to_rgba(color_bins)[:, :3]
        
# Read in the label matrix
if verbose:
    print('Reading label matrix')
label = skimage.io.imread(pathToLabelMatrix)
result = skimage.color.label2rgb(label, colors=colormap, alpha=0.3, bg_label=0, bg_color=(0, 0, 0))
[h,w] = np.shape(label)
if saveSegmentationImage:
    if verbose:
        print('Saving segmentation image to: ' + pathToSaveAnalysis + 'SegmentationImage.png')
    plt.figure(figsize=(9,9))
    plt.imshow(result)
    plt.savefig(pathToSaveAnalysis + 'SegmentationImage.png',dpi=600)
    plt.close()
    
    if verbose:
        print('Saving bw segmentation image to: ' + pathToSaveAnalysis + 'SegmentationImage_bw.png')
    boundaries = find_boundaries(label)
    toPlot = label>0
    toPlot = toPlot & ~boundaries
    
    plt.figure(figsize=(9,9))
    plt.imshow(toPlot,cmap='binary_r')
    plt.savefig(pathToSaveAnalysis + 'SegmentationImage_bw.png',dpi=600)
    plt.close()

if verbose:
    print('Reading intensity images')

intensityImagesList = [pathToIntensityDataFolder + v for v in os.listdir(pathToIntensityDataFolder) if v.endswith('.tiff')]
stainList = [v[:-5] for v in os.listdir(pathToIntensityDataFolder) if v.endswith('.tiff')]

#ST ignore the metal ending 
#todo convert the source names to ones without the marker
print("List of markers: ",stainList)
stainList = [re.sub(r'_\w+', '', stain) for stain in stainList]
print("List of markers after trimming: ",stainList)

images = np.zeros(shape=(h,w,len(intensityImagesList)))
for i, image in enumerate(intensityImagesList):
    print(image)
    images[:,:,i] = skimage.io.imread(image)

# Can do this more efficiently later using the intensity_image keyword; at the moment only supported in skimage 0.18 onwards (cluster has 0.17)
properties = skimage.measure.regionprops(label) 

if makeCellImages:
    # Check whether we're making a new image or not
    if pathToRgbImage == 'MAKE_NEW':
        if verbose:
            print('Collecting histocat images for rgb output')
        channelNames = rgb_string.split(',')
        if len(channelNames) != 3:
            print('When using MAKE_NEW, ensure --set_rgb is specified as a single string of comma separated histocat channel names (e.g., "DNA1_Ir191,PanK_Nd148,Ki67_Er168")')
            sys.exit()
        try:
            rgb_indices = [stainList.index(v) for v in channelNames]
        except Exception as e:
            print(rgb_indices)
            print('Cannot generate rgb images because:')
            print(e)
            print('Ensure --set_rgb is correctly formatted, and that each of the three given stains appears in the histocat folder.')
            sys.exit()
        # Great, now create the rgb image
        rgb = images[:,:,rgb_indices]
        rgb = np.arcsinh(rgb/arcsinhCofactor)
        maxes = np.max(np.max(rgb,axis=0),axis=0)
        rgb = rgb = rgb / maxes
        # For debugging: plt.imshow(rgb)

    else:
        # Read specified rgb image
        if verbose:
            print('Opening rgb image: ' + pathToRgbImage)
        rgb = skimage.io.imread(pathToRgbImage)
    
    if rgb.dtype == 'uint16':
        # Convert to uint8
        from skimage import exposure, img_as_ubyte
        rgb = img_as_ubyte(exposure.rescale_intensity(rgb))
    
    [h,w] = np.shape(label)
    [h2,w2] = np.shape(rgb)[0:2]
    # Rescale image if necessary
    if h != h2 or w != w2:
        if h/h2 == w/w2:
            print('Rescaling rgb image')
            rgb = rescale(rgb, h/h2, anti_aliasing=False, multichannel=True)
        else:
            # TODO - proper error detection here!
            print('rgb image and label matrix have different aspect ratios')
            sys.exit()
        
        
# We'll put the outputs into a dataframe in a minute - from a list of dicts
# One dict per row
outputList = []
if verbose:
    print('Collating cell statistics')
    
[h,w] = np.shape(label)
for i in range(len(properties)):
    # make the cell id uniquely identifyable
    cellIdImage = analysisName +'_CELL_' + str(i) + cellSegmentationImagesExtension
    cellId = analysisName +'_CELL_' + str(i)
    if makeCellImages:
        im = rgb[properties[i].bbox[0]:properties[i].bbox[2],properties[i].bbox[1]:properties[i].bbox[3],:]
        im[np.where(~properties[i].image)[0],np.where(~properties[i].image)[1],:] = backgroundColour
        # plt.imshow(im)
        plt.imsave(pathToSaveAnalysis + 'SegmentationImages/' + cellIdImage, im)
        # plt.savefig(pathToSaveAnalysis + 'SegmentationImages/' + analysisName + '_cell_' + str(i) + cellSegmentationImagesExtension, bbox_inches='tight')
        # plt.close()
        
    stats = {}
    stats['Image Name'] = cellIdImage
    stats['x'] = properties[i].centroid[1]
    stats['y'] = h-properties[i].centroid[0]
    stats['cellID'] = cellId
    stats['area'] = properties[i].area
    stats['bbox'] = properties[i].bbox
    #stats['coords'] = properties[i].coords
    stats['eccentricity'] = properties[i].eccentricity
    stats['equivalent_diameter'] = properties[i].equivalent_diameter
    stats['euler_number'] = properties[i].euler_number
    #stats['image'] = properties[i].image
    #stats['intensity_image'] = properties[i].intensity_image
    stats['label'] = properties[i].label
    stats['local_centroid'] = properties[i].local_centroid
    stats['major_axis_length'] = properties[i].major_axis_length
    stats['minor_axis_length'] = properties[i].minor_axis_length
    stats['orientation'] = properties[i].orientation
    stats['perimeter'] = properties[i].perimeter
    stats['label'] = properties[i].label
    stats['Full filepath'] = pathToSaveAnalysis + 'SegmentationImages/' + cellId

    pixel_intensities = images[properties[i].coords[:,0],properties[i].coords[:,1],:] # A nPixels by nStains array with intensities

    if doArcsinh:
        arcsinh_pixel_intensities = np.arcsinh(pixel_intensities/arcsinhCofactor)
    for j, stain in enumerate(stainList):
        # We only output one value per stain, or everyone complains that the spreadsheet is too big
        if ~doArcsinh:
            stats[stain] = pixel_intensities[:,j].mean()
        else:
            stats[stain] = arcsinh_pixel_intensities[:,j].mean()
        
            
    outputList.append(stats)
    
    
if verbose:
    print('Cell statistics complete')    
df = pd.DataFrame(outputList)
# We'll do the data trimming and normalising once we've made a dataframe, so we can use pandas and don't have to worry about it on a cell by cell basis

for j, stain in enumerate(stainList):
    # Crop the distribution
    limit = df[stain].quantile(distribution_crop_percentile)
    ind = df[stain] > limit
    df.loc[ind, stain] = limit
    
    # Now normalise
    if normalise_intensities:
        columnMax = df[stain].max()
        df.loc[:,(stain)] = df[stain]/df[stain].max()

if verbose:
    print('Writing cell statistics file to: ' + pathToSaveAnalysis + statisticsTable)
    #st
#df.to_csv(pathToSaveAnalysis + analysisName + '_cellData.csv',index=False)
#df.to_csv(pathToSaveAnalysis + analysisName + '_cellData.txt', index=False, sep='\t')
df.to_csv(pathToSaveAnalysis + statisticsTable, index=False, sep='\t')



