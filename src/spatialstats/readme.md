# Spatial Stats

These scripts analyze each ROI and produce various images and stats showing which type of cells are statistically co-related

## spatialstats.py

### Inputs

* A tab delimited text file containing at the bare minimum the following columns
    * cellID
    * cluster - contains cl01,cl02 etc
    * sample_id
    * x - the x position of the cell 
    * y - the y position of the cell

* a file that names the clusters
    * ClusterNumber
    * Annotation
```
ClusterNumber Annotation
1             T Cells
2             NK cells
.             .....
```

### Basics 
To run on the  clusters file and images produced by the pipeline :-  

```
python spatialstats/spatialstats.py \
       -i clusters.txt \
       -o <output_folder> \
       -cl harmony_pheno_cluster \
       -c annotations.txt \
       -d deepcell
```


### Arguments

* **--i --pathToData** The path to the input data. This should be a tab delimited text file with the following columns (other columns can be present, but will be ignored)
    * cellID - a unique identifier
    * cluster - the column can have any name, which is specified in the -cl argument. Cluster ids should be in the format cl01,cl02 etc. 
    * sample_id - the name of the sample that each cell belongs to e.g. HC_sample1_ROI_1
    * x - the x position of the cell 
    * y - the y position of the cell

* **--cl --clusteringToUse** The name of the column in the input column that contains the cluster ids

* **-c --cluster_annotations** A tab delimited text file which matches the cluster id with names e.g.
```
    ClusterNumber Annotation
    1             T Cells
    2             NK cells
```

* **-d --deepcellPath** The path to the deepcell folder created by the pipeline. Only required for the network analysis

* **-r --rois** The ROIs(sample ids) to be used in the analysis. By default all ROIs will be processed, but if only a selection are required, they can be specified here

```
   -r HC_Sample_1_ROI_1 HC_Sample_2_ROI_2
```

* **-f --functions** The types of analysis to carry out, which can include 
    * paircorrelationfunction
    * localclusteringheatmaps
    * celllocationmap
    * contourplots
    * quadratcounts
    * quadratcelldistributions
    * morueta-holme
    * networkstatistics

    If not specified then all the functions will be run

* **-q --quadrat_size** The size of the quadrats to use for the quadrat and morueta-holme analysis. The default value is 100 



### Outputs