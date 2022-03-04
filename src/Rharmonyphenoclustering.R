#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the integration with harmony and then clustering with Rphenograph. It will load an RData file and run harmony on the various datatransformations specified with the relevant flag and then
## run Phenograph on the results. 
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
#  	require(CATALYST, quietly=T) # this is causing harmony to quit with weird exit message, but sourcing the plotting functions script that loads CATALYST too seems to be fine (!?!?!?)
  	require(SingleCellExperiment, quietly=T)
  	require(Rphenograph, quietly=T)
  	require(harmony, quietly=T)
  	require(cowplot, quietly=T)
  	require(gridExtra, quietly=T)
  	require(ggplot2, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the integration and clustering", metavar="character"),
    make_option(c("--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--analysisName"), type="character", default=NULL, help="the name of this analysis", metavar="character"),

    make_option(c("--datatransf"), type="character", default="exprs",  help="Assay of the SCE on which to be running Harmony. (default=exprs)"),
    make_option(c("--var"), type="character", default=NULL, help="The var from the data that will be removed with Harmony. It needs to be a column in colData."),
    make_option(c("--pca"), type="logical", default="TRUE", help="Run harmony on the PCA reduction. (default=T)"),
    make_option(c("--npcs"), type="integer", default=20, help="RIf doing PCA on input matrix, number of PCs to compute. (default=20)"),
    make_option(c("--annot"), type="character", default=NULL,  help="File with annotations to be loaded for the plots.Can include more than one option with comma separated files (will be added as annotation2, etc). (optional)", metavar="character"),
    make_option(c("--k"), type="integer", default=30,  help="k parameter for Rphenograph."),

    make_option(c("--save_sceobj"), type="logical", action="store_true", help="Include flag to save the SingleCellExperiment object as well."),
    make_option(c("--draw_charts"), type="logical", default=FALSE ,help="draw charts - default is true")

); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$var),is.null(opt$outdir),is.null(opt$analysisName))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}



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
# loading the files and running the analysis

cat("Loading ",opt$infile, " \n")	
load(opt$infile)

datatransf <- opt$datatransf
k <- opt$k

if(!is.null(opt$annot)) {
	annotations <- unlist(strsplit(opt$annot, split=","))
	for (i in annotations) {
		annot <- read.table(i, header=T, stringsAsFactors=F, sep="\t")
		annot <- annot[match(colData(sce)$cellID_name,annot$cellID_name),]
		if (i==annotations[1]) {
			name <- "annotation"
		} else {
			name <- paste0("annotation", which(i==annotations))
		}
		colData(sce)[,name] <- annot$annotation
	}
}

pca <- opt$pca
cat("Running Harmony and Rphenograph for ", datatransf, ", PCA=",pca," \n")	

exp <- assay(sce, datatransf)[rowData(sce)$channel_name[rowData(sce)$marker_class=="type"],]
my_harmony_embeddings <- HarmonyMatrix(
	data_mat  = exp,
	meta_data = colData(sce)[,opt$var],
	do_pca    = pca,
	npcs=opt$npcs)

R_pheno_out <- Rphenograph(my_harmony_embeddings, k=k)
cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
## temp save:
#	save(sce, R_pheno_out, file=file.path(opt$outdir, paste0(opt$analysisName,"_Harmonysceobj_k", k,"temp.RData")))

sce[[paste0("harmony_phenograph_",datatransf)]] <- factor(membership(R_pheno_out[[2]]))
sce_umap <- uwot::umap(my_harmony_embeddings,  n_components = 3)
reducedDims(sce)[[paste0("HarmIntegr_UMAP_",datatransf)]] <- sce_umap

if  (opt$draw_charts){
  source("/stopgap/hyperion/lho/scripts/v2/hyperion/scripts/plotting_functions.R")
  sce_pheno <- data.frame(reducedDims(sce)[["UMAP"]],reducedDims(sce)[[paste0("HarmIntegr_UMAP_",datatransf)]])
  colnames(sce_pheno ) <- c("UMAP_1","UMAP_2","HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2", "HarmIntegr_UMAP_3")
  sce_pheno <- cbind(sce_pheno, sce@colData)
  
  p1a <- ggplot(as.data.frame(sce_pheno), aes(UMAP_1, UMAP_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
  p2a <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_UMAP_1, HarmIntegr_UMAP_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
  
  if(is.null(opt$annot)) {
  	previouspheno <- grep("phenograph_cluster", names(colData(sce)), value=T)
  	if (length(unique(sce_pheno[[previouspheno]]))<=30) {
  		p1b <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_1", "UMAP_2", col=previouspheno))+ scale_color_manual(values=c38) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
  		p2b <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2", col=previouspheno)) + scale_color_manual(values=c38)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
  	} else {
  		p1b <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_1", "UMAP_2", col=previouspheno)) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
  		p2b <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2", col=previouspheno)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
  	}
  } else {
  	p1b <- ggplot(as.data.frame(sce_pheno), aes(UMAP_1, UMAP_2, col=annotation))+ scale_color_manual(values=c38) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
  	p2b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_UMAP_1, HarmIntegr_UMAP_2, col=annotation)) + scale_color_manual(values=c38)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
  }
  
  if (length(unique(sce_pheno[,paste0("harmony_phenograph_",datatransf)]))<=30) {
  	p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_1", "UMAP_2", col=paste0("harmony_phenograph_",datatransf))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c38)+ guides(colour = guide_legend(override.aes = list(size=4)))
  	p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2", col=paste0("harmony_phenograph_",datatransf))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c38)+ guides(colour = guide_legend(override.aes = list(size=4)))
  } else {
  	p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_1", "UMAP_2", col=paste0("harmony_phenograph_",datatransf))) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
  	p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_UMAP_1", "HarmIntegr_UMAP_2", col=paste0("harmony_phenograph_",datatransf))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
  }
  
  pall <- plot_grid(p1a,p1b,p1c, p2a,p2b,p2c, ncol=3)
  ggsave(pall, filename=file.path(opt$outdir, paste0(opt$analysisName,"_HarmonyInt_UMAP_",datatransf,"_k", k,".pdf")), width = 30, height = 15)
  
  one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id=paste0("harmony_phenograph_",datatransf),
  	features=rowData(sce)$channel_name[rowData(sce)$marker_class=="type"], row_anno = TRUE, bars=T) 
  cat("Plotting the heatmap plots of the  phenograph clusters... \n")
  pdf(file=file.path(opt$outdir, paste0(opt$analysisName,"_HarmonyInt_heatmap_",datatransf,"_k", k,".pdf")), height=10, width=15)
  print(one_plot_heatmap)
  dev.off()
  
  one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id=paste0("harmony_phenograph_",datatransf), features = rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]) 
  cat("Plotting the expression density plots of the phenograph clusters... \n")
  pdf(file=file.path(opt$outdir, paste0(opt$analysisName,"_HarmonyInt_exprsdens_",datatransf,"_k", k,".pdf")),  height=15, width=25)
  print(one_plot_exprs)
  dev.off()
  
  one_plot_exprs_scaled <- plotClusterExprs_updated(sce, cluster_id=paste0("harmony_phenograph_",datatransf), features = rowData(sce)$channel_name[rowData(sce)$marker_class=="type"], assay="scaledtrim") 
  pdf(file=file.path(opt$outdir, paste0(opt$analysisName,"_HarmonyInt_exprsdens_",datatransf,"_k", k,"_scaledtrim.pdf")),  height=15, width=25)
  print(one_plot_exprs_scaled)
  dev.off()
  
  one_plot_abund_box <- plotAbundances_updated(sce, cluster_id=paste0("harmony_phenograph_",datatransf), group_by="sample_name", by="cluster_id" ) 
  one_plot_abund_stripe <- plotAbundances_updated(sce, cluster_id=paste0("harmony_phenograph_",datatransf), group_by="sample_name", by="sample_id") 
  cat("Plotting the cluster abundances of the phenograph clusters per sample... \n")
  pdf(file=file.path(opt$outdir, paste0(opt$analysisName,"_HarmonyInt_abund_",datatransf,"_k", k,".pdf")), height=10, width=15)
  print(one_plot_abund_box)
  print(one_plot_abund_stripe)
  dev.off()
}

##### saving 
cat("Saving everything in ", opt$outdir, "... \n")
if (opt$save_sceobj) {
	save(sce, my_harmony_embeddings, file=file.path(opt$outdir, paste0(opt$analysisName,"_Harmonysceobj_k", k,".RData")))
}

my_mat <- t(assay(sce, datatransf))

other_vars <- sce@colData
pheno_cols <- other_vars[,grepl("phenograph", names(other_vars))] 
pheno_cols <- sapply(pheno_cols, function(x) sprintf("cl%02d", x))
other_vars[,grepl("phenograph", names(other_vars))] <- pheno_cols

my_dims <- data.frame(reducedDims(sce)[[1]])
colnames(my_dims) <- paste(names(reducedDims(sce))[1], c(1,2), sep="_")
for (i in 2:length(names(reducedDims(sce)))) {
	my_dims_temp <- data.frame(reducedDims(sce)[[i]])
	colnames(my_dims_temp) <- paste(names(reducedDims(sce))[i], 1:ncol(reducedDims(sce)[[i]]), sep="_")
	my_dims <- cbind(my_dims, my_dims_temp)
}

output_table <- cbind(other_vars, my_dims, my_mat)

write.table(output_table, file=paste0(opt$outdir,"/", opt$analysisName, "_cellData_Harmonyclustered_k", k,".txt"), quote=F, row.names=F, sep="\t")
