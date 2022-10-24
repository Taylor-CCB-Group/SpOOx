## SpOOx - Spatial Omics Oxford Analysis Pipeline
This pipeline has automated the analysis of mcd image files (generated by the [Fluidigm Hyperion Imaging System](https://www.fluidigm.com/products-services/instruments/hyperion)) using [imctools](https://github.com/BodenmillerGroup/imctools) and [DeepCell](https://simtk.org/projects/deepcell) for image processing and segmentation, and custom python and R scripts to extract image intensities for a variety of cell markers. The image intensities can be further analysed using [MLV - Multi-locus View](https://www.nature.com/articles/s42003-021-02097-y), and the associated images can be analysed using [Zegami](https://zegami.com/). The pipeline itself has been written using the python pipelining tool, [Ruffus](http://www.ruffus.org.uk).

### Set-up the pipeline ##################
The following notes assume that conda has been installed and in the system path - guidance on setting up conda in a Linux environment can be found [here](https://github.com/OBDS-Training/OBDS_Open_Workshop_Materials/blob/master/1_Conda/3_Conda_intro_CCB_Rcourse.md).

A zipped version of the SpOOx code can be [downloaded](https://github.com/Taylor-CCB-Group/SpOOx/archive/refs/heads/main.zip), or the repository can be cloned using:
```
git clone https://github.com/bioinfbloke/SpOOx.git
```
The conda environment is created using:
```
conda activate base
conda env create -n hyperion -f </path/to/>hyperion.yml
conda activate hyperion
```




### Pipeline Parameters ###################
Parameters are specified in the pipeline.yaml file in the working directory

For example you will see:
`deepcell_options: --contrast 5`

If you change the contrast option this will change it for every job that is run the pipeline.

#### Specifying the queue manager

If the pipeline is not run with the --local option, then it will attempt to submit jobs to the systems cluster mechanism, which needs to be specified in the  `batch_queue` and `cluster` parameters, the following example is for slurm and a queue named 'batch'
```
batch_queue: batch
cluster:
  queue_manager: slurm
```
#### Path to the scripts

The pipeline needs to know the path to the src directory of the pipeline, which is specified in `scripts_dir`. Paths can either be absolute or relative to where the pipeline is being run. Hence if the pipeline was installed in /home/user1/SpOOx and being run from the user1 directory ,either of the following would work:-
```
scripts_dir: /home/user1/SpOOx/src
```

```
scripts_dir: SpOOx/src
```

#### Path to the markers file

`marker_file` specifies the path to the marker file (see below). Again the path can be relative to the current working directory or  absolute. Hence to use the file generated by `hyperion_pipeline.py make make_config` the following would be added
```
marker_file: markers_panel.tsv
```


### Specifying Markers ################

The pipeline needs to know which nucleus and cytoplasm markers are to be used in the segmentation step . This is a simple tab delimited text  file with 4 columns:-
* **marker_name**  the name of the marker
* **nucleus** whether this is a nucleus marker
* **cytoplasm**  whether this is a cytoplasm marker (1 or 0)
* **clustering** whether to use this marker in the clustering step
e.g.
```
marker_name	nucleus	cytoplasm	clustering
DNA1	    1	    0	        0
RAGE        0       1           1
CD56        0       0           1
.....       ..      ..          ..
```
This file can be generated by the `hyperion_pipeline.py make make_config` step of the pipeline, which extracts the name of all markers from the mcd filed. This then has to be filled in manually with the appropriate values to show which markers are to be used for segmentation (nucleus/cytoplasm) and clustering. 

Whether the file is generated by the pipeline or an existing file  with the correct markers is used, the path of the file must be specified in `pipeline.yml`
```
marker_file: markers_panel.tsv
```

### Run the pipeline ##################

Config:
```
cd <your_working_dir>

cp pipeline.yml .
# edit pipeline.yml params (eg, cluster queue name;  zegami options)
```

Input data - create a 'mcd' dir to contain all the .mcd data files to be analysed:
```
mkdir mcd 
```
Copy (or symlink) .mcd files to mcd/ dir. 

Each .mcd file should be in its own, named dir within the mcd dir. SpOOx and MDV relies on a specific naming scheme to process and parse files during the analysis. The file and dir names are constructed using three elements (\<condition\>, "sample", \<sampleId\>) separated by underscores. For example:
```
<your working dir>/mcd/COVID_sample_1/COVID_sample_1.mcd
<your working dir>/mcd/COVID_sample_2/COVID_sample_2.mcd
<your working dir>/mcd/INFLUENZA_sample_1/INFLUENZA_sample_1.mcd
<your working dir>/mcd/INFLUENZA_sample_2/INFLUENZA_sample_2.mcd
<your working dir>/mcd/INFLUENZA_sample_3/INFLUENZA_sample_3.mcd
<your working dir>/mcd/CONTROL_sample_1/CONTROL_sample_1.mcd
```
Subsequently, each region of interest (ROI) will have files named using this scheme:

\<condition\>, "sample", \<sampleId\>, "ROI", \<roi_id\>

For example:
```
COVID_sample_1_ROI_3.tiff
```
Similarly, in each row in the output data table will get a unique id:

\<condition\>, "sample", \<sampleId\>, "ROI", \<roi_id\>, "CELL", \<cell_id\>

For example:
```
COVID_sample_1_ROI_3_CELL_1.tiff
```

### Introduction to running the pipeline ###
**It is recommended you run the pipeline step by step (see below), when running for the first time to understand the process or for trouble shooting.** The SpOOx pipeline can be run at any step and all steps upstream will be attempted. For example, you may wish to exclude certain ROIs because the staining was not optimal or the region captured was too small. In this case you could run up to the generation of ROIs for each marker staining (see `tiff_to_histocat`). There are other things to bear in mind when running your samples, such as naming of markers, which if inconsistent will cause an error in the running of the pipeline. 

### Running the pipeline using a single command ###

If you already have a file that specifies which markers are to be used for segmentation and clustering, the whole pipeline can be run with
```
python SpOOx/hyperion_pipeline.py make 
```
You  will need to specify the markers file in the pipeline.yml file

Otherwise you can run up to the step that produces an empty markers file
```
python SpOOx/hyperion_pipeline.py make make_config
```
Then fill in the file and continue the rest of the pipeline

```
python SpOOx/hyperion_pipeline.py make
```


### Running the pipeline step by step ###
A basic Ruffus function 'show' can be used to describe the steps that make up the pipeline:
```
python SpOOx/hyperion_pipeline.py show
```

[imctools](https://github.com/BodenmillerGroup/imctools) works on folders not files, but folders are poor for tracking within a pipeline as their modicfication time is updated on access. To mitigate this we use a hidden dummy file (.ruffus) which is added to each directory of interest in the mcd dir for tracking purposes.

```
python SpOOx/hyperion_pipeline.py make mark_input_folders
```

The initial inputs to the pipeline are MCD image and metadata files from the hyperion and we want to convert these to OME_TIFF format files and extract each ROI for ease of processing. Each SAMPLE is put in a directory that contain the sub ROIs. There is other metadata extracted about running the machine but that is not used. It is important to make sure marker names are consistent across all of your MCD files as this will cause errors in the pipeline. Also make sure you are only using one marker panel each instance of the pipeline as this will also cause errors.
```
python SpOOx/hyperion_pipeline.py make mcd_to_tiff
```

The OME TIFF format files for each ROI are  processed to export a TIFF file based on each marker used as part of the panel (we call this the 'histocat' output because this produces tiff files which are compatible with histocat, but it is really just a convenient way for storing these files to enable processing by SpOOx).
```
python SpOOx/hyperion_pipeline.py make tiff_to_histocat
```
When using the Hyperion, some ROIs or images are not useful for further analysis. 

To remove these so they don't get analysed any further you can use:
```
python SpOOx/hyperion_pipeline.py make removebadimages
```
and this will move these into a directory called 'badimages'. Here are some examples (the command line options are shown in brackets which can be edited in the ```pipeline.yml``` file). Images/ROIs that you may want to exclude are ones: 
* created to test ablation strength for a ROI. These are often small test images e.g. 25 x 25 pixels. (--minImageWidth <size> or --minImageHeight <size>)
* used to calibrate a the machine. The name of these usually begins with a number e.g. 80ArAR,134XE or where the marker name is the same as the metal ion (--filterNonMarkers)
* that may be poor quality or abnormal in someway (--roi <list of ROIs>)

Next we create a 'markers' configuration file which is used to specify nuclear and cytoplasm channels (used in the following segmentation step). The name of the configuration file is specified using the 'marker_file' param in the 'pipeline.yml' file. 
```
python SpOOx/hyperion_pipeline.py make make_config
```

Note - the configuration SHOULD be edited at this point before proceeding to the next step. You need to specify the markers as 'Nuclear' and 'Cytoplasm' by setting the respective column / row to 1 or 0. 

Nuclear and Cytoplasm single channel images are constructed and passed to Deep Cell to create a label matrix file for each ROI.
```
python SpOOx/hyperion_pipeline.py make deepcell
```
Please note DeepCell is a great piece of software but not infallible given the wide variety of cell types that can be potentially analysed. There will almost certainly some cells that will be missed, but this is the case with any segmentation package currently. This is an area of active research!
You may need to increase the --contrast option in ```pipeline.yml``` file if you images have low signal.
Note - DeepCell requires tensorflow and defaults to looking for GPUs. Depending on your hardware configuration you may see several errors in this step. You may see errors such as "Could not load dynamic library" which can be ignored.

We next make PNG images that will be uploaded to MDV to help analyse the data generated and visualise it with the ROIs.
```
python SpOOx/hyperion_pipeline.py make roi2pngs
````

The following step uses the segmented cell masks generated by Deep Cell to extract average intensity information for each marker for each cell. This is then written to a tab separated file. Each cell has a row contains the cell id and the metadata associated with that cell including shape, marker name, etc.
```
python SpOOx/hyperion_pipeline.py make signal_extraction
```

The next step  concatenates all the data into a single file (```signalextraction/mergecellData.tab```) and also produces a table with metadata about the samples used (```metadata.tsv```)
```
python SpOOx/hyperion_pipeline.py make mergecelldata
```
The final step clusters the cells based on expression of markers (specified in the markers file) and produces various charts (see src/R_scripts). A folder *clustering* in the working directory is produced which contains various metrics, including *clusters.txt*, a tab delimited text containing cluster designations and UMAP/tSNE values  
```
python SpOOx/hyperion_pipeline.py make phenoharmonycluster
```
The recommended next step is to upload these data in Multi Dimensional Viewer (MDV). This is a web based tool that organises the data into a series of views that can be queried and visualised more easily than looking at static outputs. Details on how to upload the data into MDV are in https://github.com/Taylor-CCB-Group/SpOOx/tree/main/src/MLVUpload.

When you have investigated your data and annotated the phenotypes for each cell using the tools in MDV you can export these as annotations (see below, annotations.txt) so they can be used when you interrogate cell to cell interactions.

###  1. Spatial Stats:

To run the main spatial stats 
```
python spatialstats/spatialstats.py \
       -i clusters.txt \
       -o <output_folder> \
       -cl harmony_pheno_cluster \
       -c annotations.txt \
       -d deepcell
```
* *-i* specifies the input file, which is the clusters.txt file produced by the clustering script
* *-cl* is the column name with the cluster ids , in this case it can either be harmony_pheno_cluster of pheno_cluster
* *-c* is file containing the cluster names. The cluster id does should not contain a cl prefix as it does in the clusters.txt file e.g.
```
ClusterNumber Annotation
1             T Cells
2             NK cellls
```
* *-d* is the path to the deepcell directory produced by the pipeline, as images are required for certain parts of the analysis 


To run the PCFs averaged on groups of sample id
```
python spatialstats/AveragingPCFsByDisease.py \
-i clusters.txt \
-o <outputdir> \
-cl harmony_pheno_cluster \
-c annotations.txt \
-j conditions.json
```

The *-j* parameter specified a json file which contains the conditions(groupings) of the sample ids

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


To run the script which averages the MoruetaHolme stats by group, you first need to run the spatialstats/py script and then :-
```
python spatialstats/AveragingMoruetaHolmeByDisease.py \
-i morueta-holme \
-o <outputfolder> \
-c annotations.txt \
-j /t1-data/project/covidhyperion/sergeant/test_sard_ss/conditions1.json \

```
* *-i* is the path to the morueta-holme folder created from the spatialstats.py script
* *-j* the path to a json file describing the groupings (see above)


====================================

Optional step for visualisation using Zegami (account required: https://zegami.com/). To upload to Zegami you have to install the Zegami command line and enter your username and password to save a token that can be used by SpOOx. More details on this <a href="https://github.com/zegami/zegami-cli">here</a>.
```
python SpOOx/hyperion_pipeline.py make zegami_roi
```

====================================
