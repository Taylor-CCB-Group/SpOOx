# Uploading To MDV

The mdvupload.py script in src/MDVupload takes the various outputs of the SpOOx pipeline, compresses them and uploads it to MDV, where it is processed to produce a project where it can be interactively viewed in a web browser. 

## Quick Start

1.  Goto https://mdv.molbiol.ox.ac.uk/projects/hyperion/home (you will need to have registered and logged in) and press 'Get Token'. Then copy the token (long string of characters) that appears. Tokens last a couple of days

2. cd to the base to directory where you ran the pipeline and run the following, with 'token' being the large string you have just copied and 'my_project' a suitable name for the project.
    ```
    /path/to/SpOOx/src/MDVupload/mdvupload.py -t token -n my_project
    ```    
3. You will be sent an email with link to your project once it has been created


## What will be Uploaded?

The script will detect which parts of the pipeline have been run and upload data accordingly

### Basic data

The minimum that is required is to run the pipeline to the clustering stage
```
python SpOOx/hyperion_pipeline.py make phenoharmonycluster
```
When mdvupload,py is run it will then upload data for every cell including size, marker values, cluster designations UMAP data etc. and information about each sample. It will produce various charts and you will have the option to re-cluster choosing a subset of cells and/or markers.

### Images

MDV allows background iages to be displayed in cell centroid plots. By default, only the cell mask image is uploaded. If you run the following pipeline step then black and white images showing each channel will also be uploaded.
```
python SpOOx/hyperion_pipeline.py make roi2pngs
````
In fact you can put any image in the pngs/inmg folder with the format `<sampleid>_<imagename>.png` e.g. `COVID_SAMPLE_1_ROI_2_my_he_image.png` and it will be uploaded. The images need not be the same size or offset as the ROI e.g HE stained images of a slightly different section, as they can be resized and aligned with the cell centroids in MDV. 

### OME-TIFFS

OME-TIFFs contain information on every channel in the ROI and you can view many channels at once in different colors to produce false color images of the ROI. If you want to view the  ome-tiffs, then you need to have created pyramidal ome-tiffs
```
python SpOOx/hyperion_pipeline.py make make_pyramidal_ometiff
```
### Spatial Stats Data

If you run the spatial stats summary script (which will run the main spatial stats script if not already run) then this data will also be uploaded and can be interactively viewed
```
python SpOOx/hyperion_pipeline.py make spatialstats_summary
```
In addition, you can average the spatialstats for each condition (or another parameter)
```
python SpOOx/hyperion_pipeline.py make spatialstats_average
```

## Arguments

* **-n/--name** The name of the project - required

* **-t/--token** The token obtained from the MDV website (https://mdv.molbiol.ox.ac.uk/projects/hyperion/home) - required

* **-u/--url** The url of the site hosting MDV - default is https://mdv.molbiol.ox.ac.uk

* **-d/--directory** The base directory where the pipeline was run  - defaults to the current working directory

* **-g/--genome** The genome build  - this is not required as no genomic information is actually used (yet) defaults to hg19

* **-td/--temp_dir** The folder to crate the temporary tar ball to upload- default is _temp in the  data directory 

* **-i/--images** By default all images in the png/img folder will be uploaded, but if yu just want a certain subset, you can specify them here e.g. to specify only DNA and GranzymeB images  _-i DNA1_Ir19 GranzymeB_Er167_

* **-d/--dry_run** If true, the final tar file will be created, but not uploaded, default is false









