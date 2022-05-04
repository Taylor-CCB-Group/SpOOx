## SpOOx - Spatial Omics Oxford Analysis Pipeline

### Set-up the pipeline ##################
Note: assuming conda has been installed and in the system path
```
git clone https://github.com/bioinfbloke/SpOOx.git
```

```
conda activate base
conda env create -n hyperion -f </path/to/>hyperion.yml
conda activate hyperion
```

### Run the pipeline ##################

\# config:
```
cd <your_working_dir>

cp pipeline.yml .
# edit pipeline.yml params (eg, cluster queue name;  zegami options)
```

\# input data:
```
mkdir mcd 
```
\# copy (or symlink) mcd files to mcd/ dir. 

\# each mcd file should be in is own named dir within the mcd dir. For example:
```
<your working dir>/mcd/AB_sample_1/AB_sample_1.mcd
<your working dir>/mcd/AB_sample_2/AB_sample_2.mcd
<your working dir>/mcd/CD_sample_1/CD_sample_1.mcd
<your working dir>/mcd/CD_sample_2/CD_sample_2.mcd
<your working dir>/mcd/CD_sample_3/CD_sample_3.mcd
<your working dir>/mcd/EF_sample_1/EF_sample_1.mcd
```

### Pipeline commands:
```
python hyperion_pipeline.py show
```
```
python hyperion_pipeline.py mark_input_folders
```
```
python hyperion_pipeline.py make mcd_to_tiff
```
```
python hyperion_pipeline.py make tiff_to_histocat
```
```
python hyperion_pipeline.py make make_config
```
\# note - the config SHOULD be edited at this point before proceding to the next step
```
python hyperion_pipeline.py make deepcell
```
```
python hyperion_pipeline.py make signal_extraction
```
Optional step for visualisation using Zegami (account required: https://zegami.com/)
```
python hyperion_pipeline.py make zegami_roi
```

### Final steps after pipeline:
Note: outside conda environment:

###  1. clustering:
```
python make_metadata.py
```
```
Rscript Rphenoclustering.R \
--panel_file <markers.tsv> \
--metadata_file <metadata.tsv> \
--analysisName <analysis_name> \
--datatransf scaledtrim
--k 15 \
--q 0.001 \
--run_dimRed TRUE \
--out_dir <path_to_out_dir> \
--save_sceobj
```

### 2. spatial stats

\# [ not tested ] something like.....
```
python spatialstats.py -i StructuralIteration2.tab -o output -cl harmony_phenograph_exprs -c structural_iteration2.tab --roi HEALTHY_SAMPLE_1_ROI_1
```




