#!/usr/bin/env python
# coding: utf-8

DESC = ("Run clustering pipeline with command line inputs.")

# Imports
import os
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import harmonypy as hm
import argparse
import sys
from sklearn.preprocessing import MinMaxScaler

# define function for extracting the metadata from the cellID column of the input table
def extractmeta(df):
	if not 'cellID' in df.columns:
		print('cellID column is missing from the input file.')
	md_values = df['cellID'].str.split("_", expand=True)
	condition = md_values[0]
	sample_col = md_values.iloc[0,].str.contains('sample', case=False)
	roi_col = md_values.iloc[0,].str.contains('roi', case=False)
	sample = 'SAMPLE_' + md_values[np.where(sample_col)[0][0]+1]
	roi = 'ROI_' + md_values[np.where(roi_col)[0][0]+1]
	sample_name = condition + '_' + sample + '_' + roi
	md = pd.DataFrame({'cellID':df['cellID'], 'sample_name':sample_name, 'condition':condition, 'sample':sample, 'roi':roi})

	return md

# define function for preprocessing the data - run minmax scaling and PCA, if necessary
def preprocess(df, markersforcluster, md, minmaxtr, runPCA, quant):
	df = df[markersforcluster]
	if minmaxtr:
		# first trim
		for i in markersforcluster:
			df[i] = np.where(df[i]>df[i].quantile(quant), df[i].quantile(quant),df[i])
		# then scale
		scaler = MinMaxScaler()
		df_scaled = scaler.fit_transform(df.to_numpy())
		adata = ad.AnnData(df_scaled, obs=md, dtype=np.float32)
	else:
		adata = ad.AnnData(df, obs=md, dtype=np.float32)

	if runPCA:
		sc.tl.pca(adata, svd_solver='arpack') # this is just doing zero centering of all markers by default

	return adata


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-i', '--pathToData', 
		help = 'File path to data table (must be a tab-delimited file). Usually output from the signalextraction.',
		required = True
		) 
	parser.add_argument(
		'-o', '--output', 
	help = 'Path to write all the outputs to.',
		required = True
		) 
	parser.add_argument(
		'-p', '--panelinfo', 
		help = 'Input file with the panel information.',
		required = True
		) 
	parser.add_argument(
		'-n', '--analysisname', 
		help = 'Name of the output table with the clustered results.',
		default = 'clustereddata'
		) 
	parser.add_argument(
		'--seed', 
		default = 1,
		help = 'Provide a seed for the run.',
		required = False
		)
	parser.add_argument(
		'--minmaxtr', 
		default = False, 
		help = 'Do a min-max scaling on the data. Default=False'
		) 
	parser.add_argument(
		'--quant', 
		default = 0.999, 
		help = 'Quantile to use for trimming, set as 1 for no trimming. Default=0.999'
		) 
	parser.add_argument(
		'--runPCA', 
		default = True, 
		help = 'Run harmony on the PCA reduction. Default=True'
		) 
	parser.add_argument(
		'--npcs', 
		default = 20, 
		help = 'npcs if running PCA before Harmony. Default=20'
		) 
	parser.add_argument(
		'--runHarmony', 
		default = True, 
		help = 'Run harmony to remove unwanted variation (to be used with the variable to remove --varremove argument.) Default=False'
		) 
	parser.add_argument(
		'--varremove', 
		default = 'sample', 
		help = 'The var from the data that will be removed with Harmony. It needs to be one of condition/sample.'
		) 
	parser.add_argument(
		'--kpheno', 
		default = 30, 
		help = 'k parameter for phenograph. Default=30'
		) 
	parser.add_argument(
		'--exportAnnData', 
		default = False, 
		help = 'Export AnnData object. Default=False'
		) 

	args = parser.parse_args()

	if (os.path.exists(args.panelinfo) == False | os.path.exists(args.pathToData) == False):
			print("Panel information file (and/)or data file does not exist.")
			sys.exit()

	pathToData = args.pathToData #'/Users/erepapi/Documents/My_Documents/Fellowship/Projects/Hyperion/Hyp_COVID19/python_scripts/example_files/cellData_concat.tab'#
	pathToWriteOutput = args.output #'/Users/erepapi/Documents/My_Documents/Fellowship/Projects/Hyperion/Hyp_COVID19/python_scripts/example_files'#
	pathTomarkers = args.panelinfo #'/Users/erepapi/Documents/My_Documents/Fellowship/Projects/Hyperion/Hyp_COVID19/python_scripts/example_files/markers_panel2.tsv'
	analysisname = args.analysisname
	kpheno = args.kpheno
	runharmony = args.runHarmony
	varremove = args.varremove
	runPCA = args.runPCA
	minmaxtr = args.minmaxtr
	quant = args.quant
	seed = args.seed
	npcs = args.npcs
	exportAnnData = args.exportAnnData
	
	df = pd.read_csv(pathToData, delimiter='\t')
	markers = pd.read_csv(pathTomarkers, delimiter='\t')

	markersforcluster = markers.loc[markers['clustering']==1,'marker_name']

	if not all(markersforcluster.isin(df.columns)):
		print('Clustering markers not in dataframe. Check namings of protein markers!')
		sys.exit()

	md = extractmeta(df)
	adata_obj = preprocess(df, markersforcluster, md, minmaxtr, runPCA, quant)

	# Run harmony (optional)
	if runharmony:
		if not (varremove in adata_obj.obs.columns):
			print("Error can not remove unwanted variation with harmony, variable not valid.")
			sys.exit()
		if runPCA:
			my_harmony_embeddings = hm.run_harmony(adata_obj.obsm['X_pca'], md, varremove)
			adata_obj.obsm['X_harmony'] = my_harmony_embeddings.Z_corr.T
		else:
			my_harmony_embeddings = hm.run_harmony(adata_obj.X, md, varremove)
			adata_obj.obsm['X_harmony'] = my_harmony_embeddings.Z_corr.T			
		# run phenograph with harmony
		communities, graph, Q = sce.tl.phenograph(adata_obj.obsm['X_harmony'], k=kpheno, clustering_algo='louvain', seed=seed)
		sc.pp.neighbors(adata_obj, use_rep='X_harmony')
	else:
	# run phenograph without harmony
		if runPCA:
			communities, graph, Q = sce.tl.phenograph(adata_obj.obsm['X_pca'], k=kpheno, clustering_algo="louvain", seed=seed)
			sc.pp.neighbors(adata_obj, use_rep='X_pca')
		else:
			communities, graph, Q = sce.tl.phenograph(adata_obj.X, k=kpheno, clustering_algo="louvain", seed=seed)
			sc.pp.neighbors(adata_obj)

	adata_obj.obs['phenoGraph_cluster'] = pd.Categorical(communities)
	adata_obj.uns['phenoGraph_Q'] = Q
	adata_obj.uns['phenoGraph_k'] = kpheno

	# Do the dim reduction with UMAP and TSNE
	sc.tl.umap(adata_obj, n_components=3)
	sc.tl.tsne(adata_obj)
	
	df_out = pd.concat([md, df[markersforcluster]], axis=1)
	df_out = pd.merge(df_out, adata_obj.obs[['cellID','phenoGraph_cluster']], left_on='cellID', right_on='cellID')
	umap_dims = pd.DataFrame(adata_obj.obsm['X_umap'], index=adata_obj.obs_names)
	umap_dims.columns = ['UMAP1', 'UMAP2', 'UMAP3']
	tsne_dims = pd.DataFrame(adata_obj.obsm['X_tsne'], index=adata_obj.obs_names) #not tested
	tsne_dims.columns = ['TSNE1', 'TSNE2']
	## merge the dim red to the obj
	df_out.to_csv(os.path.join(pathToWriteOutput, (analysisname + '.txt')), header=True, sep='\t')

	if exportAnnData:
		print('Exporting data to: ', os.path.join(pathToWriteOutput, 'phenographobj.h5ad'))
		adata_obj.write(os.path.join(pathToWriteOutput, 'phenographobj.h5ad'))

if __name__ == "__main__":
	main()

