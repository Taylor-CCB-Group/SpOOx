#!/usr/bin/env python
# coding: utf-8

DESC = ("Clustering using PhenoGraph - for use in Zegami processing to define cell types from Hyperion data"
        "Saves scatterplot displaying a colour-coded 2 dimensional representation of clustering /and/ saves output table")
import argparse
import logging
import os
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
#st added
import readconfig
from pathlib import Path
import re
from pathlib import Path

from pandas_profiling import ProfileReport
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


logging_format = '%(asctime)s\t%(levelname)s\t%(message)s'
mainOutputName = "cellData.tab"


class dataset:

    def __init__(self, pathToData):
        head_tail = os.path.split(pathToData)
        file_string = head_tail[0].split('/')
        self.roi = file_string[-1]
        self.sample = file_string[-2]
        self.indication = file_string[-3]
        dirname = os.path.dirname
        abspath = os.path.abspath(pathToData)
        # parent of the specified cellData.tab 
        self.baseFolderPath = dirname(dirname(abspath))
        self.name = self.indication + "_" + self.sample + "_" + self.roi
        self.pathToCellData = pathToData
        self.pathToWriteOutput = self.baseFolderPath + '/'+args.outputDir+'/'
        self.df = pd.read_csv(self.pathToCellData,delimiter='\t', header = 0)


def main():
    #Read in celldata table and process marker input.
    # marker input is a list so take first element for path and then next is the name of the column
    #todo test for if its a file
    # first param is the name of the marker panel and second is column from that file
    marker_list = readconfig.GetMarkerList(marker_input[0], marker_input[1])

    #ST ignore the metal ending 
    #todo convert the source names to ones without the metal ending
    marker_list = [re.sub(r'_\w+', '', marker) for marker in marker_list]
    print("Marker", marker_list)
    data_X , marker_types = process_marker_input(ds.df, marker_list)
    #print("Datax_60",data_X)
    #quit()
    print("Analysing the following markers :", marker_types)
    data_X = pd.DataFrame(ds.df[marker_types].values)
    data_X.columns = marker_types
    #distribution_crop_percentile = 0.99
    #normalise_intensities = True

    # ST now do cropping / normalisation
    for j,stain in enumerate(marker_list):
        # Crop the distribution
        limit = ds.df[stain].quantile(distribution_crop_percentile)
        ind = ds.df[stain] > limit
        ds.df.loc[ind, stain] = limit
    
        # Now normalise
        if normalise_intensities:
            columnMax = ds.df[stain].max()
            ds.df.loc[:,(stain)] = ds.df[stain]/ds.df[stain].max()

    #profile = ProfileReport(data_X, title="Pandas Profiling Report", minimal=True, progress_bar=False)
    #profile.to_file("/home/s/staylor/hyperion/SpOOx/src/output.html")

    #Make anndata object of relevant data, where adata.X is the marker data, adata.obs is the image name, and adata.var is the marker labels.
    obs = pd.DataFrame(ds.df["Image Name"].values)
    adata = ad.AnnData(X = data_X, obs=obs)

    #Run phenograph, and assign output to anndata object.
    communities, graph, Q = sc.external.tl.phenograph(adata.X, k=k, directed=False, prune=False, min_cluster_size=2, jaccard=True, primary_metric='euclidean', n_jobs=- 1, q_tol=0.001, louvain_time_limit=200)
    adata.obs['phenoGraph_cluster'] = pd.Categorical(communities)
    adata.uns['phenoGraph_Q'] = Q
    adata.uns['phenoGraph_k'] = k

    #Dimensionality reduction and scatterplot.
    sc.pp.neighbors(adata)
    reduce_dimensions(adata, analysis_type)
    sc.pl.embedding(adata, basis = analysis_type, color=['phenoGraph_cluster'], projection=plot, palette='tab20', title=["Clusters"], save="_"+ds.name+"_scatterplot.png", show = False)

    #Producing output file with 3D coordinates and cluster name.
    columns_to_add = pd.DataFrame(adata.obsm["X_"+analysis_type])
    columns_to_add.columns = ["x_"+analysis_type, "y_"+analysis_type, "z_"+analysis_type]
    phenoGraph_cluster = pd.DataFrame(adata.obs['phenoGraph_cluster'])
    clusterIds=phenoGraph_cluster.values
    # convert the cluster number to a string
    newClusterIds = ["cl"+(str(clusterId).replace('[','')).replace(']','') for clusterId in clusterIds]
    columns_to_add['phenoGraph_cluster'] = newClusterIds
    #zt_out = pd.merge(ds.df, columns_to_add, on = adata.obs[0], how = 'outer')
    #st added to get rid of key_0 which keeps getting added to the output
    zt_out = pd.merge(ds.df, columns_to_add, left_index=True, right_index=True, how = 'outer')

    #zt_out.to_csv(ds.pathToWriteOutput+ds.name+"_phenograph_output.tab", na_rep = "NaN", index = False, sep ='\t')
    zt_out.to_csv(ds.pathToWriteOutput+mainOutputName, na_rep = "NaN", index = False, sep ='\t')
    #logging.info('Coordinates and phenograph clusters saved to: %s', ds.pathToWriteOutput+ds.name+"_phenograph_output.tab")
    logging.info('Coordinates and phenograph clusters saved to: %s', ds.pathToWriteOutput+mainOutputName)


    #Produce heatmap mean marker intensity for each cluster.
    adata2, n_groups = grouped_obs_mean(adata)
    logging.info("There are %s clusters present.", n_groups)
    
    #heatmap(adata, marker_types, "_big.png")
    heatmap(adata2, marker_types, "_"+ds.name+".png")
    #group_file(n_groups)


def process_marker_input(table, marker_input):
    print(">>>>>>>>>>>>>>>>>",table,">>>>>>>>>>>>",marker_input)
    #Takes input table and marker_input path, outputs (a) numpy.ndarray of marker values, (b) list of markers.
    marker_types = []
    marker_input_string = ''.join(marker_input)
    isfile = os.path.isfile(marker_input_string)
    if isfile == True:
        with open(marker_input_string) as file:
            for line in file:
                line_temp = []
                line = line.strip()
                line_temp = line.split("\t")
                if line in table.columns:
                    marker_types.append(line)
                elif line_temp[0] in table.columns:
                    for element in line_temp:
                        if element == 'type':
                            marker_types.append(str(line_temp[0]))
                        elif element == 'state':
                            logging.warning('%s is state not type. Skipping.', line_temp[0])
                else:
                    logging.warning('Row %s is invalid, or does not match the marker table. Skipping.', line)
    else:
        marker_types = marker_input
    #Check whether any markers are found.
    if not marker_types:
        logging.warning("Invalid marker input file (incompatible formatting /or/ no valid column heading in file)")
    #Marker types is a list of the markers that should be used in clustering.
    data_X = pd.DataFrame(table[marker_types].values)
    print("dataX_151",data_X)
    data_X.columns = marker_types

    #ST worth a try
    pd.DataFrame(data_X.fillna(0))
    return data_X , marker_types


def reduce_dimensions(adata, analysis_type="umap"):
    #update adata.obsm with TSNE or UMAP with (3 principle components).
    if analysis_type == "tsne":
        ###Need to update scanpy.tl.tsne with the following: https://github.com/andrea-tango/scanpy/blob/master/scanpy/tools/_tsne.py
        ###This will allow the addition of the following argument: n_components = 3
        sc.tl.tsne(adata, n_pcs = 2)
        logging.info("TSNE complete.")
    elif analysis_type == "umap":
        sc.tl.umap(adata, n_components=3)
        logging.info("UMAP complete.")
    else:
        logging.warning("Invalid analysis type.")


def heatmap(adata, marker_types, heatmap_name=".png"):
    sc.pl.heatmap(adata=adata, var_names=marker_types, groupby='phenoGraph_cluster', dendrogram = True, save=heatmap_name, cmap='RdBu_r',show = False)


def grouped_obs_mean(adata, layer=None, gene_symbols=None):
    #Function reference: https://github.com/theislab/scanpy/issues/181#issuecomment-534867254
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby("phenoGraph_cluster")
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float32),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    clusters = []
    n_groups = 0
    for group, idx in grouped.indices.items():
        n_groups += 1
        clusters.append(group)
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float32))
    
    #Make new anndata object of the group means of each cluster.
    out = pd.DataFrame(out).transpose()
    adata2 = ad.AnnData(X=out)
    adata2.obs['phenoGraph_cluster'] = pd.Categorical(clusters)

    return adata2, n_groups


parser = argparse.ArgumentParser(
        prog = "phenographclustering.py",
        description = DESC
)
parser.add_argument(
        '-i', '--input', 
        help = 'Input file name (must be a tab-delimited file)',
        required = True
) 

parser.add_argument(
        '-o', '--outputDir', 
        help = 'output directory name',
        required = True
) 

parser.add_argument(
        '-m', '--marker_input', 
        nargs = '+',
        help = 'Input a file containing column names /or/ column names separated by spaces', 
        required = False,
)    

parser.add_argument(
        '-a', '--analysis_type', 
        help = 'TSNE or UMAP',
        required = False, 
        default = "UMAP"
)

parser.add_argument(
        '-k', '--k_value', 
        help = 'Value for k',
        required = False,
        default = 30
)
parser.add_argument(
        '-p', '--plot', 
        help = '2d or 3d',
        required = False, 
        default = "2d"
)

parser.add_argument(
        '-n','--normalise_intensities', 
        dest='normalise_intensities', 
        default=False, 
        action='store_true',
        help='Rescale the intensity distributions to the range [0,1]. If applicable, rescaling happens after taking an arcsinh transform or capping intensity values to a fixed percentile')

parser.add_argument(
        '-d','--distribution_crop_percentile',
        help='After transforming the mean intensity data, this argument allows any data above the \'distribution_crop_percentile\' to be capped at that percentile, effectively trimming the data above that threshold to that level',
        type=float,
        default=0.99)

parser.add_argument(
        '-v', '--verbose', 
        action = 'store_true', 
        required = False
)

args = parser.parse_args()

#Assign command line inputs to variables.
table_input = args.input 

ds = dataset(table_input)
path = Path(ds.pathToWriteOutput)
path.mkdir(parents=True, exist_ok=True)


# logfile to go in every output directory
logFile = os.path.join(path,"README.txt")

if args.verbose:
    logging.basicConfig(
            filename = logFile,
            level = logging.DEBUG,
            format = logging_format
    )
else:
    logging.basicConfig(
            filename = logFile,
            level = logging.INFO,
            format = logging_format
    )

print("******************************")
print(args)
print("******************************")
logging.info(args)

#print(table_input)
marker_input =  args.marker_input
analysis_type = args.analysis_type
normalise_intensities = args.normalise_intensities
distribution_crop_percentile = args.distribution_crop_percentile

k = int(args.k_value)
plot = args.plot.lower()
analysis_type = analysis_type.lower()

#for input_file in table_input:
#    ds = dataset(input_file)
#    path = Path(ds.pathToWriteOutput)
#    path.mkdir(parents=True, exist_ok=True)
#    #os.mkdir(ds.pathToWriteOutput)
#    os.chdir(ds.pathToWriteOutput)
#    #Determine the directory in which to save the figures generated.
#    sc._settings.ScanpyConfig(figdir = ds.pathToWriteOutput + 'figures/')
#    main()
#    logging.info("Completed input file %s.", input_file)



#os.mkdir(ds.pathToWriteOutput)
os.chdir(ds.pathToWriteOutput)
#Determine the directory in which to save the figures generated.
sc._settings.ScanpyConfig(figdir = ds.pathToWriteOutput + 'figures/')
main()
logging.info("Completed input file %s.", table_input)


### logging thing ###
#logging.basicConfig(filename = ds.pathToWriteOutput+"/README.txt", encoding='utf-8', level=logging.DEBUG)
#logging.info(params)

