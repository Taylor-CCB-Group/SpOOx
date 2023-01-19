# Tutorial 1
In this tutorial, two mcd files HC_SAMPLE_1 and COVID_SAMPLE_6, which contain 2 and 3 images respectively are going to be processed by the pipeline. This will involve cell segmentation, followed clustering by clustering. 

##  1.Organizing the file structure

1. Create a directory ('test' in the example below) and then cd into it, this will be the run's root directory. The pipeline will be run from this directory and all results will be written to it. Then then download the sample data (https://zenodo.org/api/files/7ed9992a-b32d-47f5-8c8d-2b15caf6a862/sample_data.tar) and untar it. Next, create a subdirectory called mcd and within this create two folders COVID_SAMPLE_6 and HC_SAMPLE_1.

    ```
    mkdir test
    cd test
    wget  https://zenodo.org/api/files/7ed9992a-b32d-47f5-8c8d-2b15caf6a862/sample_data.tar
    mkdir mcd
    mkdir mcd/HC_SAMPLE_1
    mkdir mcd/COVID_SAMPLE_6 
    ```


2. Untar and move the mcd files to their corresponding folders and remove the compressed files. The markers_panel.tsv file can stay in the root directory.
    ```
    tar -xvf COVID_SAMPLE_6.tar.gz -C  mcd/COVID_SAMPLE_6/
    tar -xvf HC_SAMPLE_1.tar.gz -C  mcd/HC_SAMPLE_1/ 
    rm *.tar *.gz
    ```

3. Copy the pipeline.yml from where you installed SpOOX (/path/to/spoox) into the current directory (the run's root directory)
   ```
    cp /path/to/spoox/pipeline.yml ./
   ```

You should end up with a directory structure shown below:-
```
    test
    |--mcd
    |   |--HC_SAMPLE_1
    |   |  |--HC_SAMPLE_1.mcd
    |   |--COVID_SAMPLE_6
    |      |--COVID_SAMPLE_6.mcd
    |--markers_panel.tsv
    |--pipeline.yml
```

## 2.Changing pipeline parameters 

Edit the copied pipeline.yml file to alter parameters for the pipeline so they are suitable for this tutorial.

* **cluster** - change this if you are using a different queue manager and queue. You can run the pipeline with --local if you are not using any queue manager
* **scripts_dir**  change to /path/to/spoox/src (where /path/to/spoox is where you installed the pipeline)
* **hyperion_dir** keep as *mcd*. This is the location of the the mcd files, in this case the mcd folder that you have just populated.
* **marker_file**  keep as *markers_panel.tsv*. This specifies the markers to be used in segmentation, in this case the markers_panel.tsv that was in the tar file you extracted


## 3.Running the Pipeline To the Clustering Stage 

If you are sure your images are good and you already know which markers you are going to use for segmentation and clustering (as specified in the markers_panel.tsv) then you can run the pipeline up to the clustering stage. This will segment the cells, cluster them and produce a summary table of all cells. The segmentation will be based on the markers_panel.tsv that you downloaded.

Alternatively you can only run the pipeline to the stage of creating tiffs and then create a markers file based on the markers the pipeline has detected and remove any bad images before continuing with the pipeline (see Running the pipeline in two stages)

1. Ensure you have activated your environment

    ```
     conda activate hyperion
    ```
    or if conda is not in your path
    ```
    source/path/to/conda/install/bin/activate hyperion
    ```

2. Make sure you are in the run's root directory and run the pipeline up to the step of clustering the cells:-

   If you have you are using a queueing mechanism (which you have specified in pipeline.yml)
   ```
    python /path/to/spoox/hyperion_pipeline.py make phenoharmonycluster
   ```
   or if you are running locally 

   ```
    python /path/to/spoox/hyperion_pipeline.py make --local phenoharmonycluster
   ```
   If you get any errors on the console, look in the newly created log folder, where a log of each stage is generated and may be give a better idea of the problem

3. Once the pipeline is run there will be various folders in the run's root directory. The main summary will be clusters directory, where there will be a summary table, clusters.txt and various pdfs containing charts. 



## 4.Running the pipeline in two stages

Run the pipeline to the stage of creating tiffs. Then manual intervention to choose which markers are to be used in cell segmentation and also the opportunity

1. Activate the hyperion environment `conda activate hyperion`

2. Run the pipeline up to mark_histocat step

    ```
    python /path/to/spoox//hyperion_pipeline.py make mark_histocat
    ```

3. Create a new markers file. First,so we don't overwrite markers_panel.tsv, change the marker_file parameter in pipeline.yml to *mymarkers.tsv*. Then run

    ```
    python /path/to/spoox/hyperion_pipeline.py make make_config
    ```

    A file *mymarkers.tsv* should have been created. Fill this in appropriately - put a 1 in the appropriate columns to show which markers you want used in the segmentation (nucleus and cytoplasm) and clustering.

4. Remove any bad images which may interfere in downstream processes (in this case no images will be removed as they all pass the criteria set)
    ```
    python /path/to/spoox//hyperion_pipeline.py make removebadimages
    ```
5. Run the rest of the pipeline up to the clustering stage
    ```
    python /path/to/spoox/hyperion_pipeline.py make phenoharmonycluster
    ```