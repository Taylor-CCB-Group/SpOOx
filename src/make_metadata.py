#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import pandas as pd
from tifffile import TiffFile
from sklearn.decomposition import PCA
import shutil


parser = argparse.ArgumentParser(description='''Generate a metadata file describing the data analysed. Output will be tab-delimited: sample_id, sample_name, condition, ROI, path.
''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir', dest='indir', default = "deepcell",
                    help='the location of files used to generate this metadata')		
parser.add_argument('--outfile', dest='outfile', default = "metadata.tsv",
                    help='output file')
parser.add_argument('--mergefile', dest='mergefile', default = "signalextraction/mergecellData.tab",
                    help='merged cell data')
parser.add_argument('--panel_file', dest='panel_file',
                    help='the file to write out')		

args = parser.parse_args()

indir=args.indir
outfile=args.outfile

cwd = os.getcwd()

#get number of cells per sample
df =pd.read_csv(args.mergefile,sep="\t")
cellnums= dict(df.groupby(["sample_id"]).size().items())

pcas=None
if args.panel_file:
    mdf = pd.read_csv(args.panel_file,sep="\t")
    #get the markers for PCAs
    pca_cols= [x for x,y in zip(mdf["marker_name"],mdf["clustering"]) if y==1]
    #collapse dataframe to sample_id and take median of markers
    meds = df.groupby("sample_id").median()[pca_cols]
    #number of PCs is 4 unless markers/samples are less
    mm= min(4,min(meds.shape[0],meds.shape[1]))
    #get the PCs
    pca = PCA(n_components=mm)
    pcas = pca.fit_transform(meds)
    #dict of sample_id to pcs plus number of components
    pcas ={"pcs":{x:y for x,y in zip(meds.index,pcas)},"ncomp":mm}

    #write the marker file for easy acess
    shutil.copyfile(args.panel_file,os.path.join("signalextraction","markers.tsv"))
       
   



filedata = "sample_id\tsample_name\tcondition\tROI\tpath\tregion_width\tregion_height\tcell_number"
#add the PC headers
if pcas:
    for n in range(1,pcas["ncomp"]+1):
        filedata+=f'\tPCA_{n}'
filedata+="\n"
for d in os.listdir(indir):
    pattern = '_ROI_'
    if (re.search(pattern, d)):
        (sample_name, ROI) = re.split(pattern, d)
        #get the condition  - relies on sample names being in a specific format
        arr = d.split("_SAMPLE_")
        #sometimes they are lower case
        if len(arr)==1:
            arr=d.split("_sample_")
        condition = arr[0]
        full_path = cwd + '/' + 'signalextraction' + '/' + d
        metadata = f'{d}\t {sample_name}\t{condition}\tROI_{ROI}\t{full_path}'
        #try and get the deep cell tiff
        tiff = os.path.join(indir,d,"deepcell.tif")
        if os.path.exists(tiff):
            #extract width and height
            t=TiffFile(tiff)
            t_dim =t.pages[0].shape
            metadata +=f'\t{t_dim[0]}\t{t_dim[1]}'
        else:
            metadata+="\tNA\tNA"
        cn = cellnums.get(d,"ND")
        metadata+=f'\t{cn}'
        #add pcas if any
        if pcas:
            for n in range(0,pcas["ncomp"]):
                metadata+=f'\t{pcas["pcs"][d][n]}'

        metadata+="\n"


        filedata += metadata


        
metadata_file = cwd + '/' + outfile
f = open(metadata_file,'w')  
f.write(filedata)
f.close()
