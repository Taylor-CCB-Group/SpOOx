#!/usr/bin/env python
# coding: utf-8

DESC = ("Apply spatial statistics pipeline with command line inputs.")

# Imports
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils_alt import dataset, EfficientPCF_AtoB, plotPCFWithBootstrappedConfidenceInterval
import matplotlib.path as mpltPath
import skimage
from pathlib import Path

sns.set_style('white')
sns.set(font_scale=2)


#name of current program
prog = os.path.basename(__file__)
description = DESC

def main():     
        parser = argparse.ArgumentParser(
        )
        parser.add_argument(
                '-i', '--pathToData', 
                help = 'File path to data table (must be a tab-delimited file).',
                required = True
        ) 
        parser.add_argument(
                '-o', '--output', 
		help = 'Path to write all the outputs to.',
                required = True
        ) 
        parser.add_argument(
                '-c', '--cluster_annotations', 
                help = 'Input file with cluster annotations.',
                required = True
        )
        parser.add_argument(
                '-f', '--functions',
                nargs = '+', 
                help = 'Input the functions that are to be carried out (celllocationmap, contourplots, quadratcounts, quadratcelldistributions, paircorrelationfunction). If no argument is given, all functions will be carried out.',
                required = False
        ) 
        parser.add_argument(
                '-v', '--verbose', 
                action = 'store_true', 
                required = False
        )



        args = parser.parse_args()

        pathToData = args.pathToData 
        pathToWriteOutput = args.output
        cluster_annotations = args.cluster_annotations
        functions = args.functions
     
        if bool(functions) == False:
                print("Running all functions:")
                functions = ['celllocationmap', 'contourplots', 'quadratcounts', 'quadratcelldistributions','paircorrelationfunction','networkdiagram']

        df_annotations, ds = preprocessing(pathToData, cluster_annotations)
        #todo should not need to add a / 
        ds.pathToWriteOutput = pathToWriteOutput+"/"
        clusteringToUse = 'phenoGraph_cluster' ###ADD to cmd line?
        colors = [plt.cm.tab20(v) for v in range(len(df_annotations.ClusterNumber))]
        clusterNames = {df_annotations.ClusterNumber.iloc[v] : df_annotations.Annotation.iloc[v].strip() for v in range(len(df_annotations))}

        #make an output directory
        Path(pathToWriteOutput).mkdir(parents=True, exist_ok=True)
        print("Output file directory created.")


        #Calling the functions.
        functions = [x.casefold() for x in functions]
        if 'quadratcounts' in functions or 'quadratcelldistributions' in functions:
                r, VMRs, counts = quadratMethods(ds, df_annotations, clusteringToUse)

        for i in functions:
                if i == 'celllocationmap':
                        cellLocationMap(ds, df_annotations, clusteringToUse, clusterNames, colors)
                if i == 'contourplots':
                        contourPlots(ds, df_annotations, clusteringToUse, clusterNames)
                if i == 'quadratcounts':
                        quadratCounts(ds, df_annotations, colors, clusterNames, r, VMRs)
                if i == 'quadratcelldistributions':
                        quadratCellDistributions(ds, df_annotations, colors, clusterNames, counts)
                if i == 'paircorrelationfunction':
                        pairCorrelationFunction(ds, df_annotations, clusteringToUse, clusterNames)



def preprocessing(pathToData, cluster_annotations):
        #PART 1 - Preprocess main data table & make path to write output.
        ds = dataset(pathToData)
        print("Main data table processed.")
        
        #PART 2 - Preprocess annotation table.
        df_annotations = pd.read_csv(cluster_annotations, delimiter='\t')
        # Find the correct ROI.
        #df_annotations = df_clusterAnnotations[(df_clusterAnnotations.Indication.str.casefold()==ds.indication.casefold()) & (df_clusterAnnotations.Sample.str.casefold()==ds.sample.casefold()) & (df_clusterAnnotations.ROI.str.casefold()==ds.roi.casefold())]
        print("Annotation table processed.")
        
        return df_annotations, ds
        

def cellLocationMap(ds, df_annotations, clusteringToUse, clusterNames, colors):
        #PART 4 - Generate full map of cell locations
        plt.figure(figsize=(12,12))
        plt.title(ds.name)
        for i, cl in enumerate(df_annotations.ClusterNumber):
                p = ds.points[ds.df[clusteringToUse] == int(cl)]
                plt.scatter(p[:,0],p[:,1],color=colors[i],s=10,label=clusterNames[cl] + ' (cluster ' + str(cl) + ')')
        plt.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)
        plt.axis('square')
        plt.xlim([0,ds.domainX])
        plt.ylim([0,ds.domainY])
        plt.savefig(ds.pathToWriteOutput + ds.name + '__Overview.png',bbox_inches='tight')
        print("Cell location map completed.")


def contourPlots(ds, df_annotations, clusteringToUse, clusterNames):
        #PART 5
        for i, cl in enumerate(df_annotations.ClusterNumber):
                plt.figure(figsize=(12, 9))
                p = ds.points[ds.df[clusteringToUse] == int(cl)]
                plotting_points = pd.DataFrame(np.array(p), columns=['x', 'y'])
                ax = sns.kdeplot(data=plotting_points, shade=True, cmap="PuRd", cut=0, levels=15)#200
                ax.axis('square')
                plt.title(clusterNames[cl] + ' (cluster ' + str(cl) + ')')
                plt.xlim([0, ds.domainX])
                plt.ylim([0, ds.domainY])
                plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[cl] + ' (cluster ' + str(cl) + ')' + '_KDEplot.png',bbox_inches='tight')
        print("Contour plots completed.")


def quadratMethods(ds, df_annotations, clusteringToUse):
        #PART 6 - Time to implement some quadrat methods
        quadratEdgeLength = 100 # microns
        x = np.arange(0,ds.domainX,quadratEdgeLength)
        y = np.arange(0,ds.domainY,quadratEdgeLength)
        x = np.append(x,ds.domainX)
        y = np.append(y,ds.domainY)

        xv, yv = np.meshgrid(x,y)
        nx = len(x)
        ny = len(y)
        nQuadrats = (nx-1)*(ny-1)

        quadrats = []
        for i in range(nx-1):
                for j in range(ny-1):
                        vertices = [(x[i],y[j]), (x[i],y[j+1]), (x[i+1],y[j+1]), (x[i+1],y[j])]
                        quadrats.append(vertices)
        assert(nQuadrats == len(quadrats))

        counts = np.zeros(shape=(nQuadrats,len(df_annotations.ClusterNumber)))
        totals = np.zeros(shape=(len(df_annotations.ClusterNumber)))

        VMRs = []
        means = []
        variances = []
        for ind, cl in enumerate(df_annotations.ClusterNumber):
                p = ds.points[ds.df[clusteringToUse] == int(cl)]
                totals[ind] = len(p)
                for i in range(nQuadrats):
                        path = mpltPath.Path(quadrats[i])
                        inside2 = path.contains_points(p)
                        count = sum(inside2)
                        counts[i,ind] = count

                mean = np.nanmean(counts[:,ind])
                variance = np.nanstd(counts[:,ind])**2
                VMR = variance / mean
                VMRs.append(VMR)
                means.append(mean)
                variances.append(variance)

        r = np.corrcoef(counts.transpose())

        print("Quadrat methods completed.")

        return r, VMRs, counts


def quadratCounts(ds, df_annotations, colors, clusterNames, r, VMRs):
        #PART 7
        plt.figure(figsize=(12,9))
        g = sns.heatmap(r, cmap='RdBu_r', vmin=-1, vmax=1, annot=False, xticklabels=clusterNames.values(), yticklabels=clusterNames.values())
        plt.title('Quadrat counts: correlation coefficient')
        plt.savefig(ds.pathToWriteOutput + ds.name + '__QuadratCountsCorrelationCoefficient.png',bbox_inches='tight')

        plt.figure(figsize=(12,12))
        plt.scatter(VMRs,df_annotations.ClusterNumber,c=colors,s=400)
        # for i, cl in enumerate(clustersToCompare):
        #     plt.plot([cl,cl],[lowVMR99[i],highVMR99[i]],color=colors[i])
        plt.ylabel('Cluster ID')
        plt.xlabel('VMR')
        plt.gca().set_yticks(list(clusterNames.keys()))
        plt.gca().set_yticklabels(clusterNames.values())
        plt.plot([1,1],[np.sort(list(clusterNames.keys()))[0],np.sort(list(clusterNames.keys()))[-1]],'k:')
        plt.savefig(ds.pathToWriteOutput + ds.name + '__VarianceToMeanRatio.png',bbox_inches='tight')
        print("Quadrat counts completed.")


def quadratCellDistributions(ds, df_annotations, colors, clusterNames, counts):
        #PART 8 - Visualise distributions of cells in each quadrat
        import scipy.stats as ss

        sns.set(font_scale=1)

        for i, cl in enumerate(df_annotations.ClusterNumber):
                plt.figure(figsize=(6,6))
                b = np.arange(0,30,1)

                mu = sum(counts[:,i])/len(counts[:,i])
                #print(sum(counts[:,i]), mu, len(counts[:,i]))
                #print(ss.poisson.pmf(b, mu))
                #print(sum(counts[:,i])*ss.poisson.pmf(b, mu))
                #plt.figure()
                #plt.plot(b,ss.poisson.pmf(b, mu))
                f_exp = ss.poisson.pmf(b, mu)

                plt.plot(b, f_exp, 'bo', ms=2, label='poisson pmf')
                plt.hist(counts[:,i],color=colors[i], bins=b, density=True)
                plt.title(clusterNames[cl])

                plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[cl] + ' (cluster ' + str(cl) + ')' + '_HistogramOfQuadratCounts.png',bbox_inches='tight')

        #     f_obs = out[0]
        #     #print(out[1])
        #     chi_sq, p = ss.chisquare(f_obs, f_exp[0:-1])
        #     print(chi_sq, p)

        print("Quadrat cell distributions completed.")

def pairCorrelationFunction(ds, df_annotations, clusteringToUse, clusterNames):
        #PART 9 - Cell-cell methods
        # Comute Pair-correlation function
        maxR_mum = 100
        dr_mum = 10

        gs = np.zeros(shape=(len(df_annotations.ClusterNumber), len(df_annotations.ClusterNumber), len(np.arange(0,maxR_mum,dr_mum))))
        # lengths = []
        for a, clusterA in enumerate(df_annotations.ClusterNumber):
                for b, clusterB in enumerate(df_annotations.ClusterNumber):
                        plt.figure(figsize=(12, 9))
                        print(a,b)

                        pair = [clusterA, clusterB]
                        p0 = ds.points[ds.df[clusteringToUse] == int(pair[0])]
                        p1 = ds.points[ds.df[clusteringToUse] == int(pair[1])]
                #         lengths.append(len(p0))
                        g, radii, contributions = EfficientPCF_AtoB(p0, p1, ds.domainX, ds.domainY, dr_mum, maxR_mum)
                        gs[a, b, :] = g
                        N0 = len(p0)
                        N1 = len(p1)
                #         bootstrapSamples = CalculateBootstrapAroundCSRForPValues(N0, N1, ds.domainX, ds.domainY, dr, maxR, numBootstrapSamples=100)
                        # plt.figure(figsize=(12,9))
                        # plotPCF(plt.gca(), radii, g, label=indication)
                        plotPCFWithBootstrappedConfidenceInterval(plt.gca(), radii, g, contributions, p0, ds.domainX, ds.domainY, label=ds.indication, includeZero=True)
                        plt.title(clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]])
                #       assert(1==2)
                #       plt.show()
                        plt.savefig(ds.pathToWriteOutput + ds.name + '_' + clusterNames[pair[0]] + ' to ' + clusterNames[pair[1]] + '_PCF.png',bbox_inches='tight')
                #        plt.close()
        # plt.legend()
        # plt.close()
        # plt.show()
        # plt.close('all')
        print("Pair correlation function completed.")


        #PART 10 - Make a heatmap
        #plt.imshow(gs[:,:,1],cmap='RdBu_r',vmin=0,vmax=2)
        fig = plt.figure(figsize=(11, 9))
        radius = 1
        labs = clusterNames.values()
        g = sns.heatmap(gs[:, :, radius], cmap='RdBu_r', vmin=0, vmax=2, annot=True, xticklabels=labs, yticklabels=labs, cbar_kws={'label': '$g(r=' + str(radii[radius]) + '\mu m)$'})
        # plt.title(indications[i] + ' - ' + str(radii[radius]) + '$\mu$m')
        plt.title(ds.name + ' - ' + str(radii[radius]) + '$\mu$m')
        # fig.autofmt_xdate()
        plt.savefig(ds.pathToWriteOutput + ds.name + '_PCF_Heatmap.png',bbox_inches='tight')
        # Make graph image
        # Save everything
        print("Heatmap completed.")


main()
