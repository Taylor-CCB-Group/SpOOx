#!/usr/bin/env Rscript
##################
## Description
##
## This script runs the clustering with Rphenograph.
##
###########################
## libraries

suppressPackageStartupMessages({
  require(optparse, quietly = T)
  require(CATALYST, quietly = T)
  require(Rphenograph, quietly = T)
  require(SingleCellExperiment, quietly = T)
  require(harmony, quietly = T)
  require(Rtsne, quietly = T)
})



# parsing the arguments
option_list <- list(

    make_option(
       c("--input_file"), type = "character",
        default = "/t1-data/project/covidhyperion/sergeant/scripts/HYPETMFPV/cellData.tab",#/t1-data/project/covidhyperion/sergeant/scripts/R_scripts/test_data/cellData.tab",
      help = "the input file",
        metavar = "character"
    ),
    make_option(
       c("--output_dir"), type = "character",
        default = "/t1-data/project/covidhyperion/sergeant/scripts/HYPETMFPV/markers_panel1.tsv", help = "the output directory",
        metavar = "character"
    ),
   make_option(
      c("--panel_file"), type = "character",
      default =  "/t1-data/project/BM_hyperion/shared/data/panel1/config/markers_panel1.tsv", 
      help = "Config panel with columns marker_name and clustering
         defining markers to use for clustering.",
      metavar = "character"
   ),
   make_option(
       c("--set_seed"), type = "integer",
       help = "set the same seed to get the same results
               if it is not set, a default random on will be used",
       default = NULL
    ),
    make_option(
       c("--datatransf"), type = "character", 
        default = "scaledtrim",
        help = "Data transformation for running the clustering on
        Options: nothing=no transformation, arcsinh= arcsinh transformation,
        scaled=scaled transformation of the input,
        arcsinhscaled=scaled transformation of the arcsinh,
        scaledtrim=the scaled and trimmed values of the input (q=0.01 if not specified),
        arcsinhscaledtrim=the scaled and trimmed values of the arcsinh
        (q=0.01 if not specified). (default=exprs)"
    ),
    make_option(
       c("--q"), type = "numeric", default = 0.01,
       help = "Quantile for the trimming of the intensity values 
       (if using the scaledtrim option for the transformation). default=0.01"
    ),
    make_option(
       c("--metadata_cols"), type = "numeric", default = 2,
       help = "the number of left hand metdata columns
       only used if panel_file is NULL"
    ),

    make_option(
       c("--k"), type = "integer", default = 30,
       help = "k parameter for Rphenograph."
    ),
    make_option(
        c("--run_dimRed"), type = "character",
        default = "UMAP", action = "store_true",
        help = "UMAP,TSNE,BOTH,NONE"
    ),
   make_option(
        c("--npcs"), type = "integer", 
        default = 20, action = "store_true", 
        help = "Include flag to also run dimensionality reduction -default is True."
    ),
    make_option(
        c("--num_dimensions"), type = "integer", 
        default = 3, action = "store_true", 
        help = "Number of TNSE/UMAP dimensions"
    ),
	 make_option(
        c("--var"), type = "character", 
        default = "sample_name", action = "store_true", 
        help = "Include flag to also run dimensionality reduction -default is True."
    ),
   make_option(
      c("--save_sceobj"), type = "logical",
      default = TRUE, action = "store_true",
      help = "Include flag to save the SingleCellExperiment object as well."
   ),

   make_option(
      c("--save_trimmed_markers"), type = "logical",
      default = FALSE, action = "store_true",
      help = "save the scaled markers (not the raw ones)"
   ),
   make_option(
      c("--output_clusters_only"), type = "logical",
      default = FALSE,
      help = "output only cellID, clusters and  dimension reduction"
   ),

    make_option(
      c("--group_pcas"), type = "character",
      default = "sample_id",
      help = "Create pcas based on the media of markers for each group
             specified e.g. --group_pcas sample_id "
   ),

   make_option(c("--create_sample_columns"), type = "logical", 
   default = TRUE,
   help = "columns condition,sample_id,sample_name and ROI are
           generated from the cellID column, which must be preent,
           and in the format <cond>_<sample>_<x>_<roi>_<y>"
   )
);

 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$input_file), is.null(opt$output_dir))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call. = FALSE)
}


k <- opt$k



if (!is.null(opt$set_seed)) {
   set.seed(opt$set_seed)
}


#read in the table (allowing spaces and funny characters)
exp <- read.table(opt$input_file, header = T, check.names = FALSE, 
                   stringsAsFactors = F, sep = "\t")


#There are two options either a file showing markers (plus which ones used in clustering)
# Or just x columns 


if (! is.null(opt$panel_file)) {
    #work out markers and which to use from the panel file
    prepanel <- read.table(opt$panel_file, header = T,
                           stringsAsFactors = F, sep = "\t")

    #empty clustering values cause a problem -change to 0
    prepanel$clustering= ifelse(is.na(prepanel$clustering),"0",prepanel$clustering)

   panel <- data.frame(channel_name = prepanel$marker_name,
            marker_name = prepanel$marker_name,
            marker_class = ifelse(prepanel$clustering == "1", "type", "state")
   )
   all_markers <- panel$marker_name
   #check all specified markers are in the table 
   markers_present <- panel$marker_name %in% colnames(exp)
   absent_indexes <- which(FALSE == markers_present)
   if (length(absent_indexes) > 0){
      missing  <- panel$marker_name[absent_indexes]
      li <- paste(missing, collapse = ",")
      stop(paste(li, "specified markers are missing in the input table",
               sep = " "),
            call. = FALSE
      )
   }
   count_data <- subset (exp, select = all_markers)
   exp_data <- subset(exp, select = !(colnames(exp) %in% all_markers))
}else{
   #use all markers in the input table 
   # (all columns after the index specified by metada_cols)
   n <- opt$metadata_cols
   count_data <- subset (exp, select = c((n+1):ncol(exp)))
   exp_data <- subset (exp, select =  c(1:n))
   panel <- data.frame (channel_name = colnames(count_data),
                     marker_name = colnames(count_data),
                     marker_class = "type")
}



#create ROI,sample_name and condition (all derived from cellID)
if (opt$create_sample_columns) {
   if (!"cellID" %in% colnames(exp_data)){
        stop("No cellID column present in input table
            to extract sample information",
         call. = FALSE
      )
   }
    gcond <- function(x) {
       return(strsplit(x, split = "_")[[1]][1])
    }
    exp_data$condition <- sapply(exp_data$cellID,gcond)[]

    groi <- function(x){
        return(paste(strsplit(x, split = "_")[[1]][4],
        strsplit(x, split = "_")[[1]][5],sep = "_"))

    }
    exp_data$ROI <- sapply(exp_data$cellID,groi)[]


    gsn <- function(x) {
        return (paste(strsplit(x, split = "_")[[1]][1],
                strsplit(x, split = "_")[[1]][2],
                strsplit(x, split = "_")[[1]][3], sep = "_"))
    }

    exp_data$sample_name <- sapply(exp_data$cellID,gsn)[]
    exp_data$sample_id <- paste(exp_data$sample_name,
                              exp_data$ROI, sep = "_")

}
#create the single cell object
sce <- SingleCellExperiment(list(exprs = t(as.matrix(count_data))),
                            colData = exp_data,
                            rowData = panel)

datatransf <- "exprs"
#scale the data as specified in args
if (opt$datatransf == "scaledtrim"){
   assay(sce, "scaledtrim", FALSE) <-
      CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, opt$q)
   datatransf <- "scaledtrim"
}else if (opt$datatransf == "trim"){
   assay(sce, "scaled", FALSE) <- 
      CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, 0)
   datatransf <- "scaled"
}


cat("Running the Rphenograph analysis \n")

R_pheno_out <- Rphenograph(t(assay(sce, datatransf))[,rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]], k=k)

# Modularity is one measure of the structure of networks or graphs. It was designed to measure the strength of division of a network into modules (clusters/communities).
# Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules. Modularity is often used in optimization methods
# for detecting community structure in networks.
# It has been shown that modularity suffers a resolution limit and, therefore, it is unable to detect small communities.
cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")

sce[["pheno_cluster"]] <- factor(membership(R_pheno_out[[2]]))



#### Do Harmony
#number of pcs must be less than number of markers 
n_markers_used <-  length(which(panel$marker_class == "type"))
n_pcs <-  ifelse(opt$npcs >= n_markers_used, n_markers_used - 1, opt$npcs)
#get the raw marker values
exp <- assay(sce, "exprs")[rowData(sce)$channel_name[rowData(sce)$marker_class=="type"],]

outdir <- opt$output_dir
if (!dir.exists(outdir)) {
   dir.create(outdir, showWarnings = FALSE)
}

#calculate sample based pcas
gpcs <- opt$group_pcas
if (gpcs != "none") {
    cs_by_s <- split(seq_len(ncol(sce)), sce[[gpcs]])
    es <- as.matrix(exp)
    ms <- vapply(cs_by_s, function(cs) rowMedians(es[, cs, drop = FALSE]),
            numeric(length(which(rowData(sce)$marker_class == "type"))))
    rownames(ms) <-
        rowData(sce)$channel_name[rowData(sce)$marker_class == "type"]
    gp_pcas <- prcomp(t(ms), scale = T)
    write.table(gp_pcas$x[, 1:4],
        file = paste(outdir,"pcas.txt", sep="/"),
        col.names = NA, quote = F, sep="\t")

}




#use harmony to create PCAs that will be fed into phenograph
my_harmony_embeddings <- HarmonyMatrix(
   data_mat  = exp,
   meta_data = colData(sce)[, opt$var],
   do_pca    = T,
   npcs = n_pcs
)

#run phenograph on the PCAs and store in sce
R_pheno_out <- Rphenograph(my_harmony_embeddings, k = k)
cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
sce[["harmony_pheno_cluster"]] <- factor(membership(R_pheno_out[[2]]))


#run any dimension reduction required and store in sce
if (opt$run_dimRed == "UMAP" | opt$run_dimRed == "BOTH") {
   reducedDims(sce)[["UMAP"]] <-
      uwot::umap(my_harmony_embeddings,  n_components = opt$num_dimensions)
}

if (opt$run_dimRed == "TSNE" |  opt$run_dimRed == "BOTH") {
   tsne <- Rtsne(my_harmony_embeddings)
   reducedDims(sce)[["TSNE"]] <- tsne$Y
}

cat("Saving everything in ", opt$output_dir, "... \n")


  

if (opt$save_sceobj) {
   save(sce, my_harmony_embeddings, file = file.path(outdir, "clusters.RData"))
}

#get the sclaled or non-scaled marker values (depnds on the args)
my_mat <- t(assay(sce, ifelse(opt$save_trimmed_markers, datatransf, "exprs")))


#get the other columns and add cl to the cluster values
other_vars <- sce@colData
pheno_cols <- other_vars[,grepl("pheno_", names(other_vars))] 
pheno_cols <- sapply(pheno_cols, function(x) sprintf("cl%02d", x))
other_vars[,grepl("pheno_", names(other_vars))] <- pheno_cols
if (opt$output_clusters_only){
   #only require clusters and cellID
   cell_id <-  data.frame(other_vars[,"cellID"])
   colnames(cell_id) <- "cellID"
   output_table <- cbind(cell_id,data.frame(pheno_cols))
}else{
   output_table <- cbind(other_vars, my_mat)
}


#extract the dimension reduction (if there is any)
dim_list  <- reducedDims(sce)
my_dims <-  NULL
for (dimr in c("UMAP", "TSNE")) {
   d <- dim_list[[dimr]]
   if (is.null(d)) {
      next
   }
   mdf <-  data.frame(d)
   colnames(mdf) <-  paste(dimr, 1:opt$num_dimensions, sep = "_")
   if (is.null(my_dims)) {
      my_dims <- mdf
   }
   else{
      my_dims <- cbind(my_dims, mdf)
   }
}
if (! is.null(my_dims)){
   output_table <- cbind(output_table, my_dims)
}

write.table(output_table, file = paste0(outdir,
            "/", "clusters.txt"), quote = F, 
            row.names = F, sep = "\t"
)



