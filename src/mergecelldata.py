#!/usr/bin/env python
# coding: utf-8
import argparse
import sys
import os
import os.path
import pandas as pd
from pandas_profiling import ProfileReport


def FindAll(name, indir):
    # name of file
    # starting path
    result = []
    for root, dirs, files in os.walk(indir):
        if name in files:
            path = os.path.join(root, name)
            if "ROI" in path:
                result.append(os.path.join(root, name))
    return result


def AppendFiles(fileList, outputFile, minAreaVal, maxAreaVal):
    # take each directory and append the celldata file therein to a dataframe
    df = pd.concat((pd.read_csv(file,sep="\t") for file in fileList))
    # filter using values in 'area' column of dataframe
    df = df[(df.area >= minAreaVal) & (df.area <= maxAreaVal)]
   
    outputDir = os.path.dirname(outputFile)
    beforeClean = os.path.join(outputDir,"beforecleaning.tab")
    df.to_csv(beforeClean,sep="\t", index=False)
    # gets rid of any columns that don't line up with the header
    afterClean = df.dropna(axis='columns')
    delim ="_SAMPLE_"
    if not "_SAMPLE_" in afterClean.cellID.iloc[0]:
        delim = "_sample_"
    afterClean["sample_id"]= afterClean.apply(lambda x:  x.cellID.split("_CELL_")[0],axis=1)
    afterClean["condition"]= afterClean.apply(lambda x:  x.cellID.split(delim)[0],axis=1)
    afterClean["sample_name"]= afterClean.apply(lambda x:  x.cellID.split("_ROI_")[0],axis=1)
    afterClean["ROI"]= afterClean.apply(lambda x: "ROI_"+x.cellID.split("_CELL_")[0].split("_ROI_")[1],axis=1)


    #calculate PCAs



    #create sample metadat from cell id
   

    afterClean.to_csv(outputFile,sep="\t", index=False)
    # commented out the following line as it is nor working at the moment
    #pandasReport = os.path.join(outputDir,"pandasprofile.html")
    #profile = ProfileReport(df, title="Pandas Profiling Report", minimal=True, progress_bar=False)
    #profile.to_file(pandasReport)

def main():
    parser = argparse.ArgumentParser(description='''Generate merged cellData.tab files so can run analyses at the source/sample level based on the files one level below
    ''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--indir', dest='indir', default = "signalextraction",
                        help='initial directory level')
    parser.add_argument('--excludeList', dest='excludeList', nargs='+',
                        help='list of dirs to exclude as a pattern if you have bad samples for example e.g. --excludeList AA_SAMPLE_11_ROI_1 BB_SAMPLE_4_ROI_12')
    parser.add_argument('--infile', dest='infile', default = "cellData.tab",
                        help='the file it will look to merge')		
    parser.add_argument('--outfile', dest='outfile', default = "mergecellData.tab",
                        help='the file to write out')			
    parser.add_argument('--minarea', dest='minarea', default = "50", nargs='?', const=1, type=int,
                        help='minimum value for area of each cell')		
    parser.add_argument('--maxarea', dest='maxarea', default = "300", nargs='?', const=1, type=int,
                        help='maximum value for area of each cell')		

    args = parser.parse_args()


    infile=args.infile
    indir=args.indir
    name=args.infile
    outfile=args.outfile
    excludeList = args.excludeList
    minArea = args.minarea
    maxArea = args.maxarea

    allCellDataFiles = FindAll(name, indir)
    print("Initial data files to process",allCellDataFiles)

    if (args.excludeList):
        allCellDataFiles = [sent for sent in allCellDataFiles 
            if not any(word in sent for word in excludeList)]
        print("Remaining files after exclusion:",allCellDataFiles)

    outfile = os.path.join(indir,outfile)
    print("Writing files from:\n","\n",allCellDataFiles,"\n*** to ***\n",outfile)

    AppendFiles(allCellDataFiles, outfile, minArea, maxArea)


if __name__ == "__main__":
    main()