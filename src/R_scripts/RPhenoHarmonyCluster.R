#!/usr/bin/env Rscript
##################
## Description
##
## This script runs the clustering with Rphenograph.
##

#insstall Rphenograph as its not in conda
if (!requireNamespace("Rphenograph", quietly = TRUE)) remotes::install_github("JinmiaoChenLab/Rphenograph")



###########################
## libraries



require(optparse, quietly = T)
require(Rphenograph, quietly = T)
require(SingleCellExperiment, quietly = T)
require(harmony, quietly = T)
require(Rtsne, quietly = T)
require(scales)
require(ggplot2)
require(cowplot, quietly=T)
require(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


cluster_colors <- c(
"#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
"#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
"#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
"#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
"#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
"#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


#scales values in each row of datamatrix between 0 and 1
#optionally if q is supplied will scale between q and 1-q percentile
#which removes the skewing due to outliers
scale_and_trim <- function(scaled,q=NULL){
    for (row in 1:nrow(scaled)) {
    if (!is.null(q)){
        qv <- quantile(scaled[row,], c(opt$q,1-opt$q))
        d3 <- rescale(scaled[row,], to = c(0, 1), from = qv)
        d3 <- pmax(0, d3)
        d3 <- pmin(1, d3)
        scaled[row,]=d3
    }
    else{
            scaled[row,] <- rescale(scaled[row,],to=c(0,1))  
        }
    }
    return(scaled)
}


####################PLOTTING FUNCTIONS####################


# Heat map of marker expression (x axis) for each cluster (y axis)
# sce - the data
# features -  a vector containing the markers to use in the heat map
# cluster id - the column to use for the x category i.e. the cluster/annotation column
# q the percentile to use for scaling and trimming of the data
plotExprHeatmap_updated <- function(sce,
     features = NULL, 
     cluster_id = "cluster_id",
     assayv="exprs")

{
  

    mat <- assay(sce, assayv)[features,]
    #transpose
    mat <- t(mat)
    # add the cluster column (converting to factor if not a factor already)
    cl_col <- factor(sce[[cluster_id]])
    mat <- cbind(mat,cl_col)
    colnames(mat) <- c(features,cluster_id)
    #aggregate to median
    mat <- aggregate(mat,list(mat[,cluster_id]),median)
    #remove the group columns added by aggregate
    mat <- as.matrix(mat[3:length(mat)-1])
    #create the legends
    qs <- round(quantile(mat, c(0.01, 0.99)) * 5)/5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))

    #left annotation of clusters (tree plus colored bars)
    kids <- levels(cl_col)
    nk <- length(kids)
    myColors <- cluster_colors
    #add extra colors i no. clusters > no. colors
    if (nk > length(myColors))
        myColors <- colorRampPalette(myColors)(nk)
    myColors <- myColors[seq_len(nk)]
    names(myColors) <- kids
    df <- data.frame(cluster_id = kids)
    col <- list(cluster_id = myColors)
    left_anno <- rowAnnotation(df = df, col = col)

    #right annontaion (bars showing number cells in each cluster)
    ns  <- table(cl_col)
    fq <- round(ns/sum(ns)*100, 2)
    txt <- sprintf("%s%%(%s)", fq, names(fq))
    foo <- row_anno_text(txt, 
            just = "center", 
            gp = gpar(fontsize = 8),  
            location = unit(0.5, "npc"))
    right_anno <-  rowAnnotation(
        "n_cells" = row_anno_barplot(
            x = as.matrix(ns), width = unit(2, "cm"),
            gp = gpar(fill = "grey", col = "white"),
            border = FALSE, axis = TRUE, bar_width = 0.8),
        "foo" = foo)

    Heatmap(matrix = mat,show_column_dend = F, name ="median scaled expression",
            col = colorRamp2(seq(min(mat), 
            max(mat), l = n <- 100), 
            colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(n)),
            heatmap_legend_param = lgd_aes,
            left_annotation = left_anno, 
            right_annotation = right_anno, 
            rect_gp = gpar(col = "white"))
}


plotAbundances_updated <- function (x, cluster_id, by = c("sample_id", "cluster_id"), 
    group_by = "condition", shape_by = NULL, col_clust = TRUE, 
    distance = c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski"), linkage = c("average", "ward.D", 
        "single", "complete", "mcquitty", "median", "centroid", 
        "ward.D2"), k_pal =cluster_colors) 
{
    x$cluster_id <- x[[cluster_id]]
    by <- match.arg(by)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)
    stopifnot(is.logical(col_clust), length(col_clust) == 1)
    #shapes <- CATALYST:::.get_shapes(x, shape_by)
    if (is.null(shapes)) 
        shape_by <- NULL
    if (by == "sample_id") {
        nk <- nlevels(x[[cluster_id]])
        if (length(k_pal) < nk) 
            k_pal <- colorRampPalette(k_pal)(nk)
    }
    ns <- table(cluster_id = x[[cluster_id]], sample_id = x[["sample_id"]])
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
    if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 
        1) {
        d <- dist(t(fq), distance)
        h <- hclust(d, linkage)
        o <- colnames(fq)[h$order]
        df$sample_id <- factor(df$sample_id, o)
    }
    p <- ggplot(df, aes_string(y = "Freq")) + labs(x = NULL, 
        y = "Proportion [%]") + theme_bw() + theme(panel.grid = element_blank(), 
        strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, 
            color = NA), axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.key.height = unit(0.8, "lines"))
    switch(by, sample_id = p + (if (!is.null(group_by)) facet_wrap(group_by, 
        scales = "free_x")) + geom_bar(aes_string(x = "sample_id", 
        fill = "cluster_id"), position = "fill", stat = "identity") + 
        scale_fill_manual("cluster_id", values = k_pal) + scale_x_discrete(expand = c(0, 
        0)) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 
        100, 25)) + theme(panel.border = element_blank(), panel.spacing.x = unit(1, 
        "lines")), cluster_id = {
        p <- p + scale_shape_manual(values = shapes) + guides(col = guide_legend(order = 1, 
            override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 3)))
        if (is.null(group_by)) {
            p + geom_boxplot(aes_string(x = "cluster_id"), alpha = 0.2, 
                position = position_dodge(), outlier.color = NA) + 
                geom_point(aes_string("cluster_id", shape = shape_by), 
                  position = position_jitter(width = 0.2))
        } else {
            p + facet_wrap("cluster_id", scales = "free_y", ncol = 4) + 
                geom_boxplot(aes_string(x = group_by, color = group_by, 
                  fill = group_by), position = position_dodge(), 
                  alpha = 0.2, outlier.color = NA, show.legend = FALSE) + 
                geom_point(aes_string(x = group_by, col = group_by, 
                  shape = shape_by), position = position_jitter(width = 0.2))
        }
    })
}


plotClusterExprs_updated <- function (x, cluster_id = "cluster_id", features = "type", assayv="exprs") 
{
    x$cluster_id <- x[[cluster_id]]
    # just to order the clusters - which is not the same as in the heatmap
    #ms <- t(CATALYST:::.agg(x[features, ], "cluster_id", "median"))
    #d <- dist(ms, method = "euclidean")
    #o <- hclust(d, method = "average")$order
    cd <- colData(x)
    o <- 1:length(levels(x$cluster_id))
    es <- assay(x[features, ], assayv)
    df <- data.table::data.table(data.frame(t(es), cd, check.names = FALSE))
    df <- data.table::melt(df, id.vars = names(cd), variable.name = "antigen", 
        value.name = "expression")
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    df$cluster_id <- factor(df$cluster_id, levels = rev(c("avg", 
        levels(x$cluster_id)[o])), labels = rev(c("average", 
        paste0(names(fq), " (", fq, "%)")[o])))
    ggplot(df, aes_string(x = "expression", y = "cluster_id", 
        col = "avg", fill = "avg")) + facet_wrap(~antigen, scales = "free_x", 
        nrow = 2) + ggridges::geom_density_ridges(alpha = 0.2) + ggridges::theme_ridges() + 
        theme(legend.position = "none", strip.background = element_blank(), 
            strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45))
}



# parsing the arguments
option_list <- list(

    make_option(
        c("--input_file"), type = "character",
        help = "the input file",
        metavar = "character"
    ),
    make_option(
        c("--output_dir"), type = "character",
        default = ".",
        help = "the output directory",
        metavar = "character"
    ),
    make_option(
        c("--panel_file"), type = "character",
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
        Options: exprs=no transformation, 
        scaled=scaled transformation of the input (0 to 1)
        scaledtrim=the scaled and trimmed values of the input (q=0.01 if not specified),
        (q=0.01 if not specified). (default=scaledtrim)"
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
        default = "UMAP,TSNE", action = "store_true",
        help = "A comma delimited list of reduction TSNE,UMAP or NONE"
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
   ),
    make_option(c("--draw_charts"), type="logical", default=TRUE ,help="draw charts - default is true")
);


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (any(is.null(opt$input_file),is.null(opt$panel_file))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
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
    #gcond <- function(x) {
    #   return(strsplit(x, split = "_")[[1]][1])
    #}
    #exp_data$condition <- sapply(exp_data$cellID,gcond)[]

 
    
    exp_data$ROI <-  sapply(exp_data$cellID,function(x)substr(x,gregexpr(pattern ='_ROI_',x)[[1]]+1,gregexpr(pattern ='_CELL_',x)[[1]]-1))
    exp_data$sample_name <- sapply(exp_data$cellID,function(x) substr(x,1,gregexpr(pattern ='_ROI_',x)[[1]]-1))
    exp_data$sample_id <- sapply(exp_data$cellID,function(x) substr(x,1,gregexpr(pattern ='_CELL_',x)[[1]]-1))
}



#create the single cell object
sce <- SingleCellExperiment(list(exprs = t(as.matrix(count_data))),
                            colData = exp_data,
                            rowData = panel)

#do any scaling required
if (opt$datatransf != "exprs"){
    scaled = assay(sce, "exprs")
    for (row in 1:nrow(scaled)) {
        q= quantile(scaled[row,],c(opt$q,1-opt$q))
        if (opt$datatransf == "scaledtrim"){
            d3=rescale(scaled[row,],to=c(0,1),from=q)
            d3 = pmax(0,d3)
            d3 = pmin(1,d3)
            scaled[row,]=d3

        }
        else{
             scaled[row,]=rescale(scaled[row,],to=c(0,1))  
         }
    }
     if (opt$datatransf == "scaledtrim"){
          assays(sce)[["scaledtrim"]] = scaled
     }
     else{
         assays(sce)[["scaled"]]  = scaled
     }
}

datatransf <- opt$datatransf
#scale the data as specified in args


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
    #get indexes of each group
    cs_by_s <- split(seq_len(ncol(sce)), sce[[gpcs]])
    #if  4 groups or less PCA cannot be run
    if (length(cs_by_s)>4){
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


# get any dimension reduction specified
dtr <- strsplit(opt$run_dimRed,",")[[1]]
dtr <- dtr[dtr %in% c("UMAP","TSNE")]


for (i in 1:length(dtr)){
    if (dtr[i] =="UMAP"){
        cat("Calculating UMAP.....\n")
         reducedDims(sce)[["UMAP"]] <-
      uwot::umap(my_harmony_embeddings,  n_components = opt$num_dimensions)

    }
    if (dtr[i]=="TSNE"){
        cat("Calculating tSNE.....\n")
        tsne <- Rtsne(my_harmony_embeddings,dims= opt$num_dimensions )
        reducedDims(sce)[["TSNE"]] <- tsne$Y
    }
}


 



features <- rowData(sce)$channel_name[rowData(sce)$marker_class=="type"]

if  (opt$draw_charts){

  #draw the dim reduction charts each with sample_id if present
  #For each DR draws 3 scatter plots colored by sample_id, phenograph and harmony clusters
  if (length(dtr) !=0) {
    sce_pheno = NULL
    cnames = c()
    #create the data structures depending on which dimreds have been used
    for (dr in dtr){
        if (is.null(sce_pheno)){
            sce_pheno <- data.frame(reducedDims(sce)[[dr]])
        }
        else{
            sce_pheno = cbind(sce_pheno,reducedDims(sce)[[dr]])
        }
        cnames = c(cnames,paste0(dr,"_1"),paste0(dr,"_2"))
    }

    colnames(sce_pheno ) <- cnames
    sce_pheno <- cbind(sce_pheno, sce@colData)
    charts <- list()
    #loop though DMs and the columns to color the points
    for (dr in dtr){
        for (color_by in c("sample_id","pheno_cluster","harmony_pheno_cluster")){
             p2a <- ggplot(as.data.frame(sce_pheno), 
                aes_string(paste0(dr,"_1"), paste0(dr,"_2"), col=color_by)) +
                 geom_point(size = 0.5)+theme_bw()+ 
                 scale_color_manual(values=cluster_colors)+
                 guides(colour = guide_legend(override.aes = list(size=4)))

             charts[[length(charts)+1]] <- p2a
        }
    }
    pall <- plot_grid(plotlist=charts,ncol=3)
    ggsave(pall, filename=file.path(outdir, "DR_scatterplots.pdf"), width = 30, height = 15)
  }
  

  # Draw the heat maps
  one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id="harmony_pheno_cluster",
                    features = features, assayv=datatransf) 
  pdf(file=file.path(outdir, "Harmony_heatmap.pdf"), height=10, width=15)
  print(one_plot_heatmap)
  dev.off()

  one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id="pheno_cluster",features = features) 
  pdf(file=file.path(outdir, "Phenograph_heatmap.pdf"), height=10, width=15)
  print(one_plot_heatmap)
  dev.off()

  one_plot_abund_box <- plotAbundances_updated(sce, cluster_id="harmony_pheno_cluster", group_by="sample_name", by="cluster_id" ) 
  one_plot_abund_stripe <- plotAbundances_updated(sce, cluster_id="harmony_pheno_cluster", group_by="sample_name", by="sample_id") 
  pdf(file=file.path(outdir,"Harmony_cluster_abundance.pdf"), height=10, width=15)
  print(one_plot_abund_box)
  print(one_plot_abund_stripe)
  dev.off()

  one_plot_abund_box <- plotAbundances_updated(sce, cluster_id="pheno_cluster", group_by="sample_name", by="cluster_id" ) 
  one_plot_abund_stripe <- plotAbundances_updated(sce, cluster_id="pheno_cluster", group_by="sample_name", by="sample_id") 
  pdf(file=file.path(outdir,"Phenograph_cluster_abundance.pdf"), height=10, width=15)
  print(one_plot_abund_box)
  print(one_plot_abund_stripe)
  dev.off()

  one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id="harmony_pheno_cluster", 
  assayv=datatransf, features = features) 
  pdf(file=file.path(outdir,"Harmony_expression_density.pdf"),  height=15, width=25)
  print(one_plot_exprs)
  dev.off()

  one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id="pheno_cluster", 
  assayv=datatransf, features = features) 
  pdf(file=file.path(outdir,"Phenograph_expression_density.pdf"),  height=15, width=25)
  print(one_plot_exprs)
  dev.off()


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
for (dimr in dtr) {
   d <- dim_list[[dimr]]
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

write.table(output_table, file = file.path(outdir, "clusters.txt"), quote = F, 
            row.names = F, sep = "\t"
)

params <- DataFrame(param= names(opt),values=unlist(opt))
params[nrow(params) + 1,] = list("markers",paste(features, collapse = ","))

write.table(params, file = file.path(outdir,"params.txt"), quote = F, 
            row.names = F, sep = ":"
)



