module load python-cbrg
module load imctools/2.1.0

test=$1
echo Pipeline start

# set up
incomingDir=/stopgap/hyperion/lho/incoming/
projectName=covid2
baseDir=/t1-data/project/covidhyperion/staylor/covid2/
treeDir=$baseDir"tree"
histocatDir=$baseDir"/histocat/"
omeDir=$baseDir"/ometiff/"
confDir=$baseDir"/config/"
zegamiDir=$baseDir"/zegami/"
zegamiStackDir=$baseDir"/zegamistacks/"

#names of dirs for subdirectries in the tree
deepCellDir="deepcell/"
signalExtractionDir="signalextraction/"
spatialStatsDir="spatialstats"

scriptsDir=/home/s/staylor/hyperion/scripts

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
	find $incomingDir -maxdepth 1 -type d -exec imctools mcdfolder-to-imcfolder {} $omeDir \;
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
	find $treeDir -name "histocat" -exec $test $scriptsDir/"deepercell.py" --cyt $confDir/cyto.txt --nuc $confDir/nuc.txt --outdir {}/../deepcell/ --dir {} \; 
}


function ZegamiROI() {
	# generate all the images
	python $scriptsDir/roi2zegami.py --indir $histocatDir --outdir $zegamiDir --yaml $zegamiDir/zegami.yaml --name $projectName
	# upload to cloud
	module load zegami-cli/
	zeg create collections --project XY6PsWre --config $zegamiDir/zegami.yaml
}


function ZegamiROIStacked() {
	# zegami stacks
	python $scriptsDir/flat2stacks.py -i $zegamiDir/img/ -o $zegamiStackDir -y $zegamiStackDir/"zegami.yaml" -n $projectName -d $projectName -m $markersToShow -t $zegamiStackDir/"zegami.tsv"
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
	find $treeDir -name "deepcell" -type d -exec $test $scriptsDir/"signalextraction.py" {}/deepcell.tif {}/../histocat/ MAKE_NEW {}/../$signalExtractionDir/ --verbose \;
}

function ZegamiUploadCellImages() {
	echo ZegamiUploadCellImages
	find $treeDir -name "signalextraction" -exec $test $scriptsDir/"uploadcellstozegami.py" {}/ XY6PsWre --verbose \;

}

function MergeSignalExtractionAtDifferentScales() {
	echo MergeSignalExtractionAtDifferentScales
	# uses signalextraction at the ROI level
	# this runs SignalExtraction at the SAMPLE level 
	#find $treeDir -name "SAMPLE*" -exec $test $scriptsDir/"mergecelldata2.py"  --indir {} \;
	# run at the condition level, there is no pattern to match so need to do it like this
	find $treeDir -maxdepth 1 -type d -exec $test $scriptsDir/"mergecelldata2.py"  --indir {} \;
	# on condition level
	#find $treeDir -maxdepth 2 -exec $test $scriptsDir/"mergecelldata2.py"  --indir {} \;
	# ./mergecelldata2.py  --indir /t1-data/project/covidhyperion/staylor/covid2/tree/COVID/SAMPLE_10
}

function Phenograph() {
	echo Phenograph
	# this runs Phenograph at the SAMPLE level since uses maxdepth 3, this needs to be run after 
	find $treeDir -name "signalextraction" -exec $test $scriptsDir/"phenographclustering.py" -i {}/cellData.tab -o clustering -m /t1-data/project/covidhyperion/staylor/covid2/config/markers_full.tsv clustering -a umap -k 30 -p 3d \;
}

function SpatialStatistics() {
	echo SpatialStatistics 
	find $treeDir -name "clustering" -exec $test $scriptsDir/"spatialstats.py" -i {}/cellData.tab -c $confDir/exampleclusterannotation.tab -f celllocationmap paircorrelationfunction --output {}../$spatialStatsDir
}


#CreateDirs
#MCDToTiff
#Tiff2Histocat
#Renamer
#ZegamiROI
#MakeResultsTree
#DeepCell
#ZegamiROIStacked
#SignalExtraction
#MergeSignalExtractionAtDifferentScales
#ZegamiUploadCellImages
#Phenograph
#TODO
#Harmony
#SpatialStatistics

