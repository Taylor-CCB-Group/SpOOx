#!/usr/bin/env Rscript
##################
## Description 
##
## This script runs the clustering with Rphenograph.
##
###########################
## libraries 

suppressPackageStartupMessages({
  require(optparse, quietly=T)
  require(CATALYST, quietly=T)
  require(flowCore, quietly=T)
  require(Rphenograph, quietly=T)
  require(SingleCellExperiment, quietly=T)
  require(ggplot2, quietly=T)
  require(cowplot, quietly=T)
  require(gridExtra, quietly=T)
  require(purrr, quietly=T)
  require(plyr, quietly=T)
  require(tidyr, quietly=T)
  require(Hmisc, quietly=T)
  require(dplyr, quietly=T)
})

# parsing the arguments
option_list = list(

    make_option(c("--analysisName"), type="character", default="test", help="the name of this analysis", metavar="character"),
    make_option(c("--panel_file"), type="character", default = "/project/taylorlab/sergeant/rtest/markers.tsv", help="Config panel with columns marker_name and clustering defining markers to use for clustering.", metavar="character"),
    make_option(c("--metadata_file"), type="character", default= "/project/taylorlab/sergeant/rtest/metadata.txt", help="Text file with metadata information (sample_id, sample_name, condition, ROI)", metavar="character"),
    make_option(c("--only_include"), type="character", default=NULL, help="Name of specific ROI/sample/condition (sample_id/sample_name/condition) to process, if NULL it will process everything in the metadata file.", metavar="character"),
    make_option(c("--samples_to_exclude"), type="character", default=NULL, help="Sample_id to exclude from the analysis(comma separated).", metavar="character"),

    make_option(c("--filters"), type="character", default=NULL, help="Comma separated filters to be evaluated. The variables need to be column names (attention to the case). (e.g. Area>30,Eccentricity<=0.4)", metavar="character"),
    make_option(c("--datatransf"), type="character", default="scaledtrim",  help="Data transformation for running the clustering on. Options: nothing=no transformation, arcsinh= arcsinh 
		transformation, scaled=scaled transformation of the input, arcsinhscaled=scaled transformation of the arcsinh, scaledtrim=the scaled and trimmed values of the input (q=0.01 if not specified),
		arcsinhscaledtrim=the scaled and trimmed values of the arcsinh (q=0.01 if not specified). (default=exprs)"),
    make_option(c("--q"), type="numeric", default=0.01,  help="Quantile for the trimming of the intensity values (if using the scaledtrim option for the transformation). default=0.01"),

    make_option(c("--k"), type="integer", default=30,  help="k parameter for Rphenograph."),
    make_option(c("--run_dimRed"), type="logical", default=TRUE, action="store_true", help="Include flag to also run dimensionality reduction -default is True."),
    make_option(c("--save_sceobj"), type="logical", default= TRUE ,action="store_true", help="Include flag to save the SingleCellExperiment object as well."),
    make_option(c("--out_dir"), type="character",default="/project/taylorlab/sergeant/rtest", help="the output directory - default is to calculate it from the input.", metavar="character"),
    make_option(c("--draw_charts"), type="logical", default=FALSE ,help="draw charts - default is true"),
    make_option(c("--use_subdirectory"), type="logical", default=FALSE ,help="If true  celldata.tab should e in signal extraction subdirectory in the specified path. default is FASLE")

); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$panel_file),is.null(opt$metadata_file),is.null(opt$analysisName))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

if (is.null(opt$run_dimRed)) opt$rerun_dimRed <- F
if (is.null(opt$save_sceobj)) opt$save_sceobj <- F
k <- opt$k
opt$analysisName <- gsub(",", replacement="_", opt$analysisName)
if (!is.null(opt$filters)) {
	filters <- gsub(" ", replacement="", opt$filters)
	filters <- unlist(strsplit(filters, split=","))
}
if (opt$draw_charts){
  source("/stopgap/hyperion/lho/scripts/v2/hyperion/scripts/plotting_functions.R")
}


### add the colour scheme that I like better:
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
### import and process the data

cat("Starting the preprocessing part... \n")
prepanel <- read.table(opt$panel_file, header=T, stringsAsFactors=F,sep="\t")
#prepanel$marker_name <- gsub("-", replacement=".", prepanel$marker_name)
#if cluster column has 1 use (type), anyother values will not be used
panel <- data.frame(channel_name=prepanel$marker_name, marker_name=prepanel$marker_name, marker_class=ifelse(prepanel$clustering=="1", "type", "state"))

md_file <- read.table(opt$metadata_file, header=T, stringsAsFactors=F, sep="\t")
md_file$n_cells <- NA
if (!is.null(opt$samples_to_exclude)) {
	cat("Excluding samples: ",opt$samples_to_exclude,"\n")
	samples_to_exclude <- unlist(strsplit(opt$samples_to_exclude, split=","))
	md_file <- md_file[!(md_file$sample_id %in% samples_to_exclude),]
}

if (!is.null(opt$only_include)) {
	cat("Only running samples: ",opt$only_include,"\n")
	samples_to_include <- unlist(strsplit(opt$only_include, split=","))
	if(grepl("ROI", samples_to_include[1])) {
		md_file <- md_file[md_file$sample_id %in% samples_to_include,]
		if (nrow(md_file)==0) stop("Sample not valid. Check the metadata file to make sure the naming scheme is correct.")
		pathToSaveAnalysis <- paste0(md_file$path,"clustering")
		#if (!dir.exists(pathToSaveAnalysis)) dir.create(pathToSaveAnalysis, showWarnings = FALSE)
	} else if (grepl("SAMPLE", toupper(samples_to_include[1]))) {
		md_file <- md_file[md_file$sample_name %in% samples_to_include,]
		if (nrow(md_file)==0) stop("Sample name not valid. Check the metadata file to make sure the naming scheme is correct.")
		if (length(samples_to_include)>1 & length(unique(md_file$condition))>1) {
			temp_dir <- strsplit(md_file$path[1], split="/")[[1]]
			pathToSaveAnalysis <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering",paste(samples_to_include,collapse="_")), collapse="/")
			temp_path <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering"), collapse="/")
			if (!dir.exists(temp_path)) dir.create(temp_path, showWarnings = FALSE)
		} else if (length(samples_to_include)>1 & length(unique(md_file$condition))==1) {
			pathToSaveAnalysis <-  paste0(strsplit(md_file$path[1], split="SAMPLE")[[1]][1],"clustering")
		} else {
			pathToSaveAnalysis <-  paste0(strsplit(md_file$path[1], split="ROI")[[1]][1],"clustering")
		}
		#if (!dir.exists(pathToSaveAnalysis)) dir.create(pathToSaveAnalysis, showWarnings = FALSE)
	} else {
		md_file <- md_file[md_file$condition %in% samples_to_include,]
		if (nrow(md_file)==0) stop("Condition not valid. Check the metadata file to make sure the naming scheme is correct.")
		if (length(samples_to_include)==1) {
			pathToSaveAnalysis <-  paste0(strsplit(md_file$path[1], split="SAMPLE")[[1]][1],"clustering")
		} else {
			temp_dir <- strsplit(md_file$path[1], split="/")[[1]]
			pathToSaveAnalysis <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering",paste(samples_to_include,collapse="_")), collapse="/")
			temp_path <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering"), collapse="/")
			if (!dir.exists(temp_path)) dir.create(temp_path, showWarnings = FALSE)
		}
		#if (!dir.exists(pathToSaveAnalysis)) dir.create(pathToSaveAnalysis, showWarnings = FALSE)
	}
} else {
	temp_dir <- strsplit(md_file$path[1], split="/")[[1]]
	samples_to_include <- opt$analysisName
	pathToSaveAnalysis <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering",paste(samples_to_include,collapse="_")), collapse="/")
	temp_path <-  paste(c(temp_dir[1:(length(temp_dir)-3)],"clustering"), collapse="/")
	if (!dir.exists(temp_path)) dir.create(temp_path, showWarnings = FALSE)
	#if (!dir.exists(pathToSaveAnalysis)) dir.create(pathToSaveAnalysis, showWarnings = FALSE)
}

if (! is.null(opt$out_dir)){
  pathToSaveAnalysis=opt$out_dir
  
}
if (!dir.exists(pathToSaveAnalysis)) dir.create(pathToSaveAnalysis, showWarnings = FALSE)

antib_data <- data.frame()
celldata <- data.frame()

### loading in the data
for (i in md_file$sample_id) { 

	cat("Uploading sample ", i, "\n")		

	filename <- file.path(md_file$path[md_file$sample_id==i],ifelse(opt$use_subdirectory,"signalextraction",""),"cellData.tab")
	#filename <- file.path(md_file$path[md_file$sample_id==i],"SignalExtraction",paste(i,"cellData.csv",sep="_"))

	exp1 <- read.table(filename, header=T, check.names=FALSE, stringsAsFactors=F, sep="\t")
	exp1$cellID_name <- gsub(".png", replacement="", exp1$"Image Name")

	### applying filters
	if (!is.null(opt$filters)) {
		for (i in filters) {
			var <- strsplit(i, split="<|>|<=|>=")[[1]][1]
			num <- gsub("[A-Za-z]", "", i)

			filter <- paste0("exp1[['",var,"']]",num)
			index <- eval(parse(text=filter))
			cat("keeping only cells for which ", filter, " : removing ", length(which(!index)), " cells. \n")
			exp1 <- exp1[index,]
		}
	}

	expr_only <- exp1
	names(expr_only) <- sapply(names(expr_only), function(x) strsplit(x,"_")[[1]][1])

	antib_data_temp <- select(expr_only, panel$channel_name)
	index <- grep("Full filepath", names(exp1))
	celldata_temp <- cbind(exp1[,1:index],select(expr_only, setdiff(names(expr_only),panel$channel_name)))

	antib_data <- rbind(antib_data, antib_data_temp)
	celldata <- rbind.fill(celldata, celldata_temp)

	md_file[md_file$sample_id==i,"n_cells"] <-  nrow(exp1)

}

row.names(antib_data) <- celldata$cellID_name
row.names(celldata) <- celldata$cellID_name
row.names(panel) <- panel$channel_name
celldata$sample_id <- sapply(celldata$"Image Name", function(x) strsplit(x, split="_CELL")[[1]][1])

### creating the sce object

if (opt$datatransf=="arcsinh") {
	sce <- SingleCellExperiment(list(counts=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	sce <- CATALYST:::.transform(sce, 5) # asinh(value/5) - stored in the "exprs" slot of the assay()
	datatransf <- "exprs"
} else if (opt$datatransf=="nothing") {
	sce <- SingleCellExperiment(list(exprs=t(as.matrix(antib_data))),
    		colData=celldata,
    		rowData=panel,
    		metadata=list(experiment_info=md_file)
	)
	datatransf <- "exprs"
} else if (opt$datatransf=="scaled") {
	sce <- SingleCellExperiment(list(exprs=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	assay(sce, "scaled", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, 0)
	datatransf <- "scaled"
} else if (opt$datatransf=="arcsinhscaled") {
	sce <- SingleCellExperiment(list(counts=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	sce <- CATALYST:::.transform(sce, 5) # asinh(value/5) - stored in the "exprs" slot of the assay()
	assay(sce, "scaled", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, 0)
	datatransf <- "scaled"
} else if (opt$datatransf=="scaledtrim") {
	sce <- SingleCellExperiment(list(exprs=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	assay(sce, "scaledtrim", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, opt$q)
	datatransf <- "scaledtrim"
} else if (opt$datatransf=="arcsinhscaledtrim") {
	sce <- SingleCellExperiment(list(counts=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	sce <- CATALYST:::.transform(sce, 5) # asinh(value/5) - stored in the "exprs" slot of the assay()
	assay(sce, "scaledtrim", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, opt$q)
	datatransf <- "scaledtrim"
} else {
  stop(paste("The argument datatransf needs to be one of arcsinh/arcsinhscaled/scaled/arcsinhscaledtrim/scaledtrim/nothing."), call.=FALSE)
}

colData(sce)$sample_name <- sapply(colData(sce)$sample_id, function(x) sce@metadata$experiment_info$sample_name[sce@metadata$experiment_info$sample_id==x])
colData(sce)$ROI <- sapply(colData(sce)$sample_id, function(x) sce@metadata$experiment_info$ROI[sce@metadata$experiment_info$sample_id==x])
colData(sce)$condition <- sapply(colData(sce)$sample_id, function(x) sce@metadata$experiment_info$condition[sce@metadata$experiment_info$sample_id==x])

### Clustering

cat("Running the Rphenograph analysis \n")	

R_pheno_out <- Rphenograph(t(assay(sce, datatransf))[,rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]], k=k)

# Modularity is one measure of the structure of networks or graphs. It was designed to measure the strength of division of a network into modules (clusters/communities).
# Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules. Modularity is often used in optimization methods
# for detecting community structure in networks.
# It has been shown that modularity suffers a resolution limit and, therefore, it is unable to detect small communities.
cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")

sce[[paste0("phenograph_cluster_", datatransf,"_k", k)]] <- factor(membership(R_pheno_out[[2]]))

### Plots 
if  (opt$draw_charts){
  one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id=paste0("phenograph_cluster_", datatransf,"_k", k),
  	features=rowData(sce)$channel_name[rowData(sce)$marker_class=="type"], row_anno = TRUE, bars=T) 
  cat("Plotting the heatmap plots of the  phenograph clusters... \n")
  pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_Rpheno_heatmap_k", k,".pdf")), height=10, width=15)
  print(one_plot_heatmap)
  dev.off()

  one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id=paste0("phenograph_cluster_",datatransf,"_k", k), features = rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]) 
  cat("Plotting the expression density plots of the phenograph clusters... \n")
  pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_Rpheno_exprsdens_k", k,".pdf")), height=15, width=25)
  print(one_plot_exprs)
  dev.off()

  one_plot_exprs_scaled <- plotClusterExprs_updated(sce, cluster_id=paste0("phenograph_cluster_",datatransf,"_k", k), features = rowData(sce)$channel_name[rowData(sce)$marker_class=="type"], assay="scaledtrim") 
  pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_Rpheno_exprsdens_k", k,"_scaledtrim.pdf")), height=15, width=25)
  print(one_plot_exprs_scaled)
  dev.off()
  
  if (nrow(sce@metadata$experiment_info)>1) {
  	one_plot_abund_box <- plotAbundances_updated(sce, cluster_id=paste0("phenograph_cluster_",datatransf,"_k", k), group_by="sample_name", by="cluster_id" ) 
  	one_plot_abund_stripe <- plotAbundances_updated(sce, cluster_id=paste0("phenograph_cluster_",datatransf,"_k", k), group_by="sample_name", by="sample_id") 
  	cat("Plotting the cluster abundances of the phenograph clusters per sample... \n")
  	pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_Rpheno_clusterabund_k", k,".pdf")), height=10, width=15)
  	print(one_plot_abund_box)
  	print(one_plot_abund_stripe)
	dev.off()
}

#  require(ggcorrplot, quietly=T) missing
#	if (length(unique(colData(sce)$sample_name))>1) {
#		for (j in unique(colData(sce)$sample_name)) {
#			expr_mat <- t(assay(sce, "exprs")[,colData(sce)$sample_name==j])[,rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]]
#			mydata.rcorr = rcorr(expr_mat)
#			pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_corrplot_", j,".pdf")), height=10, width=10)
#			ggcorrplot(mydata.rcorr$r, hc.order = TRUE, outline.col = "white") + ggtitle(i)
#			dev.off()
#		}
#	}

}

#### Do UMAP and TSNE

if (opt$run_dimRed) {

	cat("Running dimentionality reduction for ", datatransf, " ... \n")
	cat("starting with TSNE... \n")
	sce <- runDR(sce, dr = "TSNE",  features = "type", assay=datatransf)
	cat("UMAP... \n")
	sce <- runDR(sce, dr = "UMAP", features = "type", assay=datatransf)

	sce_pheno <- data.frame(reducedDims(sce)[["TSNE"]], reducedDims(sce)[["UMAP"]])
	colnames(sce_pheno) <- c("TSNE1","TSNE2", "UMAP1","UMAP2")

	pheno_ks <- paste0("phenograph_cluster_", datatransf,"_k", k)
	sce_pheno <- cbind(sce_pheno, sce@colData)
	sce_pheno[[pheno_ks]] <- as.factor(sce_pheno[[pheno_ks]])
	if (opt$draw_charts){
  	if (length(levels(sce_pheno[[pheno_ks]]))<=27) {
  		one_plot_tsne <- ggplot(sce_pheno, aes_string(x="TSNE1", y="TSNE2", col=pheno_ks)) + geom_point(size = 0.5) + theme_classic() + scale_color_manual(values=c38)
  		one_plot_umap <- ggplot(sce_pheno, aes_string(x="UMAP1", y="UMAP2", col=pheno_ks)) + geom_point(size = 0.5) + theme_classic() + scale_color_manual(values=c38)
  	} else {
  		one_plot_tsne <- ggplot(sce_pheno, aes_string(x="TSNE1", y="TSNE2", col=pheno_ks)) + geom_point(size = 0.5) + theme_classic() 
  		one_plot_umap <- ggplot(sce_pheno, aes_string(x="UMAP1", y="UMAP2", col=pheno_ks)) + geom_point(size = 0.5) + theme_classic() 
  	}
  
  	tsne_sample <- ggplot(sce_pheno, aes_string(x="TSNE1", y="TSNE2", col="sample_id")) + geom_point(size = 0.5) + theme_classic() + scale_color_manual(values=c38)
  	umap_sample <- ggplot(sce_pheno, aes_string(x="UMAP1", y="UMAP2", col="sample_id")) + geom_point(size = 0.5) + theme_classic() + scale_color_manual(values=c38)
  
  	pdf(file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_Rpheno_TSNEUMAP_k", k,".pdf")), height=8, width=10)
  	    print(one_plot_tsne)
  	    print(one_plot_umap)
  	    print(tsne_sample)
  	    print(umap_sample)
  	dev.off()
	}

	my_dims <- data.frame(reducedDims(sce)[[1]])
	colnames(my_dims) <- paste(names(reducedDims(sce))[1], c(1,2), sep="_")
	for (i in 2:length(names(reducedDims(sce)))) {
		my_dims_temp <- data.frame(reducedDims(sce)[[i]])
		colnames(my_dims_temp) <- paste(names(reducedDims(sce))[i], c(1,2), sep="_")
		my_dims <- cbind(my_dims, my_dims_temp)
	}

}

### saving 
cat("Saving everything in ", pathToSaveAnalysis, "... \n")
if (opt$save_sceobj) {
	save(sce, file=file.path(pathToSaveAnalysis, paste0(opt$analysisName,"_sceobj_k", k,".RData")))
}

my_mat <- t(assay(sce, "exprs"))

other_vars <- sce@colData
pheno_cols <- other_vars[,grepl("phenograph", names(other_vars))] 
pheno_cols <- sapply(pheno_cols, function(x) sprintf("cl%02d", x))
other_vars[,grepl("phenograph", names(other_vars))] <- pheno_cols

if (opt$run_dimRed) {
	output_table <- cbind(other_vars, my_dims, my_mat)
} else {
	output_table <- cbind(other_vars, my_mat)
}

write.table(output_table, file=paste0(pathToSaveAnalysis,"/", opt$analysisName, "_cellData_clustered.txt"), quote=F, row.names=F, sep="\t")
