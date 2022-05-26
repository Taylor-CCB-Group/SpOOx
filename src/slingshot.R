#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the pseudotime analysis with slingshot. It will load an RData file and run slingshot on the phenograph clusters specified with the relevant flag.
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
  	require(SingleCellExperiment, quietly=T)
  	require(slingshot, quietly=T)
  	require(harmony, quietly=T)
  	require(cowplot, quietly=T)
  	require(dplyr, quietly=T)
  	require(ggplot2, quietly=T)
  	require(viridis, quietly=T)
  	require(readxl, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the integration and clustering", metavar="character"),
    make_option(c("--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--analysisName"), type="character", default=NULL, help="the name of this analysis", metavar="character"),

    make_option(c("--var"), type="character", default=NULL, help="The vars (comma separated) from the data that will be removed with Harmony. They need to be columns in colData."),
    make_option(c("--pca"), type="logical", default="TRUE", help="Run harmony on the PCA reduction. (default=T)"),
    make_option(c("--npcs"), type="integer", default=20, help="If doing PCA on input matrix, number of PCs to compute. (default=20)"),

    make_option(c("--annot"), type="character", default=NULL,  help="Excel file with annotations to be loaded for the analysis. The clusters need to be ordered (from 1 to n) and it needs to have a column named ClusterID and one named Annotations. (optional file)", metavar="character"),
    make_option(c("--start"), type="character", default=NULL,  help="Cluster to use as starting point in the plotting of slingshot. (optional)", metavar="character"),
    make_option(c("--sht"), type="integer", default=1,  help="Sheet number to be read from the annotations excel file (optional)", metavar="character"),
    make_option(c("--stage_info"), type="character", default=NULL,  help="Excel file with stages to be loaded for the density plots. (optional file)", metavar="character"),
    make_option(c("--sht_stage"), type="integer", default=1,  help="Sheet number to be read from the stage info excel file. Needs columns Stage and sample_id. (optional)", metavar="character"),
    make_option(c("--clusters"), type="character", default=NULL,  help="Phenograph cluster ids on which the analysis will be run (comma separated). If NULL, running the pseudotime analysis on all clusters.", metavar="character"),
    make_option(c("--markers"), type="character", default=NULL,  help="Markers which will be used for the analysis (comma separated). Need to be an exact match with the names of the markers in the column marker_name from the rowData of the obj. If NULL, using the same markers that phenograph was run.", metavar="character")

); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$outdir),is.null(opt$analysisName))) {
  print_help(opt_parser)
  stop("Arguments missing.", call.=FALSE)
}

source("/stopgap/hyperion/lho/scripts/v2/hyperion/scripts/plotting_functions.R")

###########################
c38 <- c("dodgerblue2", "#E31A1C", # red
                "green4", 
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black", "gold1",
                "skyblue2", "#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "darkgrey", "khaki2",
                "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                "darkturquoise", "green1", "yellow4", "yellow3",
                "darkorange4", "grey","darkorchid4","brown",
                "azure3","burlywood1","aquamarine","cadetblue","coral2",
                "darkolivegreen1","darkgoldenrod2","brown1","darkgreen","blueviolet","darkred"
              )

###########################
# loading the files and setting the parameters

cat("Loading ",opt$infile, " \n")	
load(opt$infile)

clusters <- unlist(strsplit(opt$clusters, split=","))
markers <- opt$markers

if(!is.null(opt$annot)) {
	annotations <- read_excel(opt$annot, sheet=opt$sht)
	if (any(annotations$ClusterID!=1:length(annotations$ClusterID))) stop("Clusters not ordered!")
	colData(sce)$annotations <- annotations$Annotations[colData(sce)$harmony_phenograph_exprs]
}

if(!is.null(opt$stage_info)) {
	othervars <- read_excel(opt$stage_info, sheet=opt$sht_stage)
}

if(!is.null(markers)) {
	rowData(sce_pseudo)$marker_class <- "state"
	rowData(sce_pseudo)$marker_class[rowData(sce_pseudo)$marker_name %in% markers] <- "type"
}

vars_harmony <- unlist(strsplit(opt$var, split=","))
if (length(vars_harmony>1)) {
	vars_use <- vars_harmony
} else {
	vars_use <- NULL
}

# running the analysis - first redoing harmony on the subset of clusters and then slingshot

if(!is.null(clusters)) {
	sce_pseudo <- sce[,(colData(sce)$harmony_phenograph_exprs %in% clusters)]
	exp <- assay(sce_pseudo, "exprs")[rowData(sce_pseudo)$channel_name[rowData(sce_pseudo)$marker_class=="type"],]
	my_harmony_embeddings_subset <- HarmonyMatrix(
		data_mat = exp,
		meta_data = colData(sce_pseudo)[,vars_harmony],
		vars_use = vars_use,
		do_pca = opt$pca,
		npcs = opt$npcs)
	sce_umap <- uwot::umap(my_harmony_embeddings_subset, n_components = 2)
	reducedDims(sce_pseudo)[["HarmIntegr_UMAP_subset"]] <- sce_umap
	sce_pheno <- data.frame(reducedDims(sce_pseudo)[["HarmIntegr_UMAP_subset"]])
	slingshot_redDim <- "HarmIntegr_UMAP_subset"
	colnames(sce_pheno) <- c("HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2")
	sce_pheno <- cbind(sce_pheno, sce_pseudo@colData)
	p1 <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_UMAP_1, HarmIntegr_UMAP_2, col=annotations)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))+ scale_color_manual(values=c38)
	ggsave(p1, filename=paste0(opt$outdir,"/", opt$analysisName, "_UMAPsubset.pdf"), width = 8, height = 7)
	p2 <- p1 + facet_wrap(~annotations)
	ggsave(p2, filename=paste0(opt$outdir,"/", opt$analysisName, "_UMAPsubset_facet.pdf"), width = 10, height = 7)

	if(!is.null(opt$stage_info)) {
		colData(sce_pseudo)$stage <- factor(othervars$Stage[match(colData(sce_pseudo)$sample_id, othervars$sample_id)])
		pdf(file=paste0(opt$outdir,"/", opt$analysisName, "_density.pdf"), width = 10, height = 8)
			print(plot_density(sce_pseudo, nrow=3, umap_reduction="HarmIntegr_UMAP_subset", var="stage"))
		dev.off()
	}
} else {
	sce_pseudo <- sce
	slingshot_redDim <- "HarmIntegr_UMAP"
}

sce_pseudo <- slingshot(sce_pseudo, clusterLabels = 'annotations', reducedDim = slingshot_redDim)
if (!is.null(opt$start)) lin1 <- getLineages(sce_pseudo, colData(sce_pseudo)$annotations, start.clus = opt$start, reducedDim = slingshot_redDim)

pdf(file=paste0(opt$outdir,"/", opt$analysisName, "_slingshot.pdf"), width = 8, height = 7)
	plot(reducedDims(sce_pseudo)[[slingshot_redDim]], col = c38[as.factor(sce_pseudo$annotations)], pch=16, asp = 1, cex=0.5)
	if (!is.null(opt$start)) {
		lines(SlingshotDataSet(lin1), lwd=2, col='black')
	} else {	
		lines(SlingshotDataSet(sce_pseudo), lwd=2, col='black')
	}
dev.off()

save(sce_pseudo, file=paste0(opt$outdir,"/", opt$analysisName, "_slingshot.RData"))

