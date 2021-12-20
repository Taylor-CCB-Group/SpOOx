module load python-cbrg
module load imctools/2.1.0
module load zegami-cli/1.5.1 

test=$1

echo Pipeline start

# set up
incomingDir=/stopgap/hyperion/lho/incoming/
projectName=covid2
baseDir=/project/covidhyperion/shared/data/panel2/
treeDir=$baseDir"tree"
histocatDir=$baseDir"/histocat/"
omeDir=$baseDir"/ometiff/"
confDir=$baseDir"/config/"
zegamiDir=$baseDir"/zegami/"
mcdDir=$baseDir"/mcd/"
zegamiStackDir=$baseDir"/zegamistacks/"
#markerFile=$confDir"/markers.tsv"
markerFile=$confDir"/markers_panel2.tsv"

#names of dirs for subdirectries in the tree
deepCellDir="deepcell/"
signalExtractionDir="signalextraction/"
spatialStatsDir="spatialstats/"

scriptsDir=/home/s/staylor/hyperion/SpOOx/src

# the names of the markers to show in the stacked view in a file
markersToShow=$baseDir"/config/markers_zegami_stacks.txt"

function CreateDirs() {
# creates inital directories for uncompressing data etc
	mkdir $baseDir
	mkdir $omeDir $deepcellDir $histocatDir $deepcellDir $zegamiDir $zegamiStackDir
}

function MCDToTiff() {
	#input
	# find all the directories in incoming and make tiffs
	find $mcdDir -maxdepth 1 -type d -exec imctools mcdfolder-to-imcfolder {} $omeDir \;
}

function Tiff2Histocat() {
	find $omeDir -maxdepth 1 -type d -exec imctools omefolder-to-histocatfolder {} $histocatDir \;
}

function Renamer() {
	#the metal suffixes e.g. Sm147 in CCR2_Sm147.tiff often get incorrectly entered and are generally not needed, so remove them at this stage
	#using built in BASH cli tool prename NOT USED but noted here
	#find $histocatDir -maxdepth 1 -type d -exec prename 's/_\w+?\.tiff/\.tiff/' {}/* \;
	#uses lookup table to rename the histocat directories into something sensible
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern COVID1panel2 --disease_sample COVID_SAMPLE_1
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern CUN10_panel2 --disease_sample COVID_SAMPLE_10
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern CUN11_panel2 --disease_sample COVID_SAMPLE_11
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern CUN16_Panel2 --disease_sample COVID_SAMPLE_16
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern CUN17_panel2 --disease_sample COVID_SAMPLE_17
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern Healthy_0028000B_Panel2 --disease_sample HEALTHY_SAMPLE_1
	python $scriptsDir"/histocatdirrenamer.py" --indir $histocatDir --pattern Tonsilcontrol --disease_sample TONSIL_SAMPLE_1
}

function MakeConfig(){
	# make the template that *MUST* be edited when the QC has been done and the markers that are 
	# successfull are known
	echo MakeConfig
	$test $scriptsDir/"writeConfig.py" --indir $histocatDir --outfile $markerFile
}

function MakeResultsTree() {
	# makes the output results tree and populates ROIs with histocat folder
	# which is a symlink
	# needs inDir and outDir arguments for 
	#inDir = "/t1-data/project/covidhyperion/staylor/covid2/histocat"
	#outDir = "/t1-data/project/covidhyperion/staylor/covid2/tree"
	echo MakeResultsTree
	$test $scriptsDir"/dirtree.py" --indir $histocatDir --outdir $treeDir
}

function DeepCell() {
	#deepcell
	echo DeepCell
	# need to go into the folder to get the basename
	#find $treeDir -name "histocat" -exec mkdir {}/../$deepCellDir \; 
	#find $treeDir -name "histocat" -exec $test $scriptsDir/"deepercell.py" --cyt $confDir/cyto.txt --nuc $confDir/nuc.txt --outdir {}/../deepcell/ --dir {} \; 
	find $treeDir -name "histocat" -exec $test $scriptsDir/"deepercell.py" --markerfile $markerFile --outdir {}/../deepcell/ --indir {} \; 
}


function ZegamiROI() {
	# generate all the images using the histocat directory so this must be up to date
	$test $scriptsDir/roi2zegami.py --indir $histocatDir --outdir $zegamiDir --yaml $zegamiDir/zegami.yaml --name $projectName
	# upload to cloud
	module load zegami-cli/
	zeg create collections --project XY6PsWre --config $zegamiDir/zegami.yaml
}


function ZegamiROIStacked() {
	# zegami stacks
	echo ZegamiROIStacked
	# copy any of the deepcell images to the zegami stack dir
	find $treeDir -name "deepcell" -exec $test $scriptsDir/"deepercellcopyimages.py" --zegamidir $zegamiStackDir --indir {} \; 
	$test $scriptsDir/flat2stacks.py -i $zegamiDir/img/ -o $zegamiStackDir -y $zegamiStackDir/"zegami.yaml" -n $projectName -d $projectName -m $markerFile -t $zegamiStackDir/"zegami.tsv"
	zeg create collections --project XY6PsWre --config $zegamiStackDir/zegami.yaml
}	



function SignalExtraction() {
	echo SignalExtraction
	#So if you want to use an existing one to make cell images (e.g., Ilastik probability mask, coregistered H&E slide, whatever else you might like) 
	#then you can just point the script at it. But if you call MAKE_NEW and pass a string of stain names then it'll generate a new RGB image 
	#(so you'd do "MAKE_NEW" and then set "DNA1_Irxxx,CD68_Whatever,PanK_Etc' and it'll split that up and use those histocat images to make an RGB image
	# not make sure add / on directory names in the cli 
	# todo use os.path.join rather than +
	# this runs Signal extraction initially at the ROI level
	# signalextraction_1
	#find $treeDir -name "deepcell" -type d -exec $test $scriptsDir/"signalextraction.py" {}/deepcell.tif {}/../histocat/ MAKE_NEW {}/../$signalExtractionDir/ --verbose --distribution_crop_percentile 1.0 \; 
	# signalextracttion
	find $treeDir -name "deepcell" -type d -exec $test $scriptsDir/"signalextraction.py" {}/deepcell.tif {}/../histocat/ MAKE_NEW {}/../$signalExtractionDir/ --verbose --distribution_crop_percentile 1.0 --doArcsinh \;
}

function ZegamiUploadCellImages() {
	echo ZegamiUploadCellImages
	find $treeDir -name "signalextraction" -exec $test $scriptsDir/"uploadcellstozegami.py" {}/ XY6PsWre --verbose \;

}

function MergeSignalExtractionAtDifferentScales() {
	echo MergeSignalExtractionAtDifferentScales
	# this merges the signalextraction directories at the SAMPLE level 
	find $treeDir -name "SAMPLE*" -exec $test $scriptsDir/"mergecelldata.py"  --indir {} \;
	# run at the condition level, there is no pattern to match so need to do it for each condition
	find $treeDir -name COVID -type d -exec $test $scriptsDir/"mergecelldata.py"  --indir {} \;
	find $treeDir -name HEALTHY -type d -exec $test $scriptsDir/"mergecelldata.py"  --indir {} \;
	find $treeDir -name TONSIL -type d -exec $test $scriptsDir/"mergecelldata.py"  --indir {} \;
}

function Phenograph() {
	echo Phenograph
	# this runs Phenograph at the SAMPLE level since uses maxdepth 3, this needs to be run after 
	#find $treeDir -name "signalextraction" -exec $test $scriptsDir/"phenographclustering.py" -i {}/cellData.tab -o clustering -m /t1-data/project/covidhyperion/staylor/covid2/config/markers_full.tsv clustering -a umap -k 30 -p 3d \;
	find $treeDir -name "signalextraction" -exec $test "$scriptsDir/"phenographclustering.py -i {}/cellData.tab -o clustering -m $markerFile clustering -a umap -k 30 -p 3d  --distribution_crop_percentile 0.99 --normalise_intensities \;

}

function SpatialStatistics() {
	echo SpatialStatistics 
	#$scriptsDir/"spatialstats.py" -i /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/clustering/cellData.tab -c /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/clustering/annotations.tab --output /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/spatialstats/ 
	$scriptsDir/"spatialstats.py" -i /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/clustering/cellData.tab -c /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/clustering/annotations.tab -f paircorrelationfunction --output /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/spatialstats/ 


	#find $treeDir -name "clustering" -exec $test $scriptsDir/"spatialstats.py" -i {}/cellData.tab -c /t1-data/project/covidhyperion/shared/data/panel2/tree/COVID/clustering/annotations.tab --output {}../$spatialStatsDir \;
}


#CreateDirs
#MCDToTiff
#Tiff2Histocat
#Renamer
#ZegamiROI
#MakeConfig
#MakeResultsTree
#DeepCell
ZegamiROIStacked
#SignalExtraction
#MergeSignalExtractionAtDifferentScales
#ZegamiUploadCellImages
#Phenograph
#SpatialStatistics
#TODO
#Harmony


