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
    * ClusterNumber - the number (name) of the cluster in the clusters column (specified by the -cl argument)
    * Annotation - the name of the cluster
```
e.g.
ClusterNumber Annotation
cl01             T Cells
Cl02             NK cells
.             .....
```
If an annotation file is not supplied than the names in the cluster column will be used.

```

```


### Basics 
To run on the  clusters file and images produced by the pipeline :-  
```
 python spatialstats.py -i clusters/clusters.txt \
                        -o output_folder \
                        -cl harmony_pheno_cluster \
                        -d deepcell
```


### Arguments

* **--i --pathToData** The path to the input data (required). This should be a tab delimited text file with the following columns (other columns can be present, but will be ignored)
    * cellID - a unique identifier
    * cluster - the column can have any name, which is specified in the -cl argument. Cluster ids should be in the format cl01,cl02 etc. 
    * sample_id - the name of the sample that each cell belongs to e.g. HC_sample1_ROI_1
    * x - the x position of the cell 
    * y - the y position of the cell

* **-o --output** A folder where all the output is written to (required)

* **-cl --clusteringIdColumn** The name of the column in the input  file that contains the cluster ids (required)

* **-c --cluster_annotations** A tab delimited text file which matches the cluster id with names e.g.
```
    ClusterNumber Annotation
    cl01             T Cells
    cl02             NK cells
```
if not supplied the cluster names/ids in the specified column will be used

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



## average_by_condition.py

### Basics 
Produces spatial stats based on ROIs grouped by condition (once the main spatial stats script has been run)
```
python spatialstats/average_by_condition.py \
        -i clusters/clusters.txt \
        -p <spatial_stats_outputdir>
        -o <outputdir> \
        -cl harmony_pheno_cluster \
        -j conditions.json
```


### Arguments


* **-p** is the output folder of the main spatial stats script -the -o argument of spatialstats.py (required)
* **-j** specifies a json file which contains the conditions(groupings) of the sample ids (required)
```
{
	"conditions":{
        "Healthy":[
	        "HC_sample_1_ROI_6",
            "HC_sample_3_ROI_3"
        ],
        "Diseased":[
            "DS_sample_1_ROI_1",
            "DS_sample1_ROI_2"
        ]
    }
}

```
* **--i --pathToData** The path to the input data used to create the original statistics -see above (required)
* **-cl --clusteringIdColumn** The name of the column in the input file that contains the cluster ids (required)
* **-o --output** A folder where all the output is written to (required)




## summary.py

to summarise the main spatial stats:-
```
python spatialstats/summary.py -p <spatial_stats_ouputdir>
```


