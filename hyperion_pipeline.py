'''
Spatial proteomics pipeline for Hyperion (imaging mass cytometry) from Fluidigm
usingup to 40 different marker antibodies with metal tags

input data format - mcd files (binary) one per folder

dependencies:
imctools (module)
zegami-cli (module)
'''
from ruffus import *
from cgatcore import pipeline as P
import sys
import os
from imctools.converters import mcdfolder2imcfolder
from imctools.converters import ome2histocat



print("*** Pipeline started ***")
# print current working directory
print("Current working directory: {}".format(os.getcwd()))

# import pipeline parameters from yaml format configuration file in working directory
PARAMS = P.get_parameters("pipeline.yml")

#make a log directory if it doesn't exist
logDir="log"
if not os.path.exists(logDir):
    os.makedirs(logDir)
    

# imctools works on folders not files but folders are poor for tracking as their modicfication time is updated on access
# To mitigate this we track a hidden dummy file (.ruffus) added to each directory of interest for tracking purposes
@originate(PARAMS["hyperion_dir"]+"/*/.ruffus")
def mark_input_folders(outfiles):
    statement = '''find %(hyperion_dir)s/* -maxdepth 1 -type d -exec touch {}/.ruffus \;;'''
    P.run(statement, without_cluster=True)

# mcd_to_tiff
#
@follows(mkdir ("ometiff"),mark_input_folders)
@transform (PARAMS["hyperion_dir"]+"/*/.ruffus", regex(PARAMS["hyperion_dir"]+r'/(.*)/.ruffus'), r'ometiff/\1/.ruffus')
def mcd_to_tiff (infile, outfile):
    indir = os.path.dirname(infile)
    statement = '''python %(scripts_dir)s/parse_mcd.py -i %(indir)s -o ometiff && touch %(outfile)s >> log/mcd_to_tiff.log 2>&1'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# tiff_to_histocat
@follows(mkdir ("histocat"))
@subdivide(mcd_to_tiff, regex(r'ometiff/(.*)/.ruffus'), r'histocat/\1_ROI_*/.ruffus')
def tiff_to_histocat (infile, outfiles):
    indir = os.path.dirname(infile)
    statement = '''python %(scripts_dir)s/parse_ome.py -i %(indir)s -o histocat >> log/tiff_to_histocat.log 2>&1'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# mark_histocat
# touching the file to get around ruffus not being able to handle directories
@follows(tiff_to_histocat)
@originate("histocat/.ruffus")
def mark_histocat(outfile):
    statement = '''touch %(outfile)s'''
    P.run(statement, without_cluster=True)

# removebadimages
# remove images that are too small (test ones) and ones that do not have biological content
@transform(mark_histocat, regex(r'histocat/*/.ruffus'), r'badimages/.ruffus')
def removebadimages(infile,outfile):
    indir = os.path.dirname(infile)
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/removebadimages.py -i %(indir)s -o %(outdir)s %(removebadimages_options)s  > log/removebadimages.log 2>&1;touch badimages/.ruffus'''
    P.run(statement, without_cluster=True)

# make_config file
@follows(tiff_to_histocat)
@originate(PARAMS["marker_file"])
def make_config(outfile):
    statement = '''python %(scripts_dir)s/writeconfig.py --indir histocat --outfile %(marker_file)s >> log/make_config.log  2>&1'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# deepcell
@transform(tiff_to_histocat, regex(r'histocat/(.*)/.ruffus'), r'deepcell/\1/deepcell.tif')
def deepcell (infile, outfile):
    indir = os.path.dirname(infile)
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/deepercell.py --markerfile %(marker_file)s 
                --indir %(indir)s --outdir %(outdir)s %(deepcell_options)s >> log/deepcell.log 2>&1'''
    P.run(statement, job_queue=PARAMS['batch_queue'])


# make pngs for visualisation
@follows(deepcell,mkdir("pngs"))
def roi2pngs():
    statement = '''python %(scripts_dir)s/roi2png.py 
                    --indir histocat --outdir pngs >> log/roi2pngs.log 2>&1'''
    P.run(statement, without_cluster=True)
    statement = '''python %(scripts_dir)s/deepcell2png.py 
                    --indir deepcell --outdir pngs/img >> log/deepcell2pngs.log 2>&1'''
    P.run(statement, without_cluster=True)

# signal_extraction
@transform(deepcell, regex(r'deepcell/(.*)/deepcell.tif'), r'signalextraction/\1/cellData.tab')
def signal_extraction (infile, outfile):
    indir = os.path.dirname(infile)
    analysis_name = indir.replace('deepcell/','')
    intensity_data_dir = indir.replace('deepcell','histocat') + '/'
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/signalextraction.py 
                %(infile)s
                %(intensity_data_dir)s
                MAKE_NEW
                %(outdir)s
                --analysisName %(analysis_name)s
                %(signal_extraction_options)s >> log/signalextraction.log 2>&1'''
    P.run(statement, job_queue=PARAMS['batch_queue'])


@follows(signal_extraction)
@originate("signalextraction/mergecellData.tab")
def mergecelldata(outfile):
    statement = '''python %(scripts_dir)s/make_metadata.py >> log/mergecelldata.log 2>&1'''
    P.run(statement, without_cluster=True)
    statement = '''python %(scripts_dir)s/mergecelldata.py >> log/mergecelldata.log 2>&1'''
    P.run(statement, without_cluster=True)
    

# phenoclustering
@follows(mergecelldata)
@originate("clustering/clusters.txt")
def phenoharmonycluster(outfile):
    statement = '''Rscript %(scripts_dir)s/R_scripts/RPhenoHarmonyCluster.R
                --panel_file %(marker_file)s
                --input_file signalextraction/mergecellData.tab
                --output_dir clustering
                %(phenograph_options)s >> log/clustering.log 2>&1'''
    P.run(statement,without_cluster=True)

# Generate a graph of the pipeline
#pipeline_printout_graph('flowchart.png', 'png', [phenoharmonycluster], no_key_legend = False)

### Zegami (optional) so commented out ######################################################################

# zegami_roi
# @follows(tiff_to_histocat, mkdir("zegami"))
# @transform (mark_histocat, regex(r'histocat/.ruffus'), r'zegami/zegami.yaml')
# def zegami_roi (infile, outfile):
#     statement = '''python %(scripts_dir)s/roi2zegami.py 
#                     --indir histocat --outdir zegami --yaml zegami/zegami.yaml --name %(zegami_project)s'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])

# # zegami_roi_stacked
# @follows(mkdir("zegamistack"))
# @transform (deepcell, regex(r'deepcell/(.*)/deepcell.tif'), r'zegamistack/zegami.yaml')
# def zegami_roi_stacked (infile, outfile):
#     indir = os.path.dirname(infile)
#     outdir = os.path.dirname(outfile)
#     statement = '''python %(scripts_dir)s/deepercellcopyimages.py 
#                     --zegamidir zegamistacks --indir %(indir)s'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])
    

# @transform (zegami_roi, regex(r'zegami/zegami.yaml'), r'zegami/.ruffus')
# def zegami_create (infile, outfile):
#     statement = '''zeg create collections --project %(zegami_collection)s --config zegami/zegami.yaml
#                     && touch zegami/.ruffus'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])


# @transform (signal_extraction, 
#             regex(r'signalextraction/(.*)/cellData.tab'), 
#             r'signalextraction/\1/.zegami_upload')
# def zegami_upload (infile, outfile):
#     indir = os.path.dirname(infile)
#     statement = '''python %(scripts_dir)s/uploadcellstozegami.py %(indir)s %(zegami_collection)s --verbose 
#                     && touch %(outfile)s'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])


# @follows(zegami_roi)
# @transform (zegami_roi_stacked, regex(r'deepcell/(.*)/.ruffus'), r'zegamistack/zegami.yaml')
# def zegami_roi_flat2stacked (infile, outfile):
#     indir = os.path.dirname(infile)
#     base_dir = indir.replace('deepcell','')
#     outdir = os.path.dirname(outfile)
#     statement = '''python %(scripts_dir)s/flat2stacks.py 
#                     -i zegami/img/ -o zegamistack -y zegamistack/zegami.yaml 
#                     -n %(zegami_project)s -d %(zegami_project)s 
#                     -m %(marker_file)s -t zegamistack/zegami.tsv'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])


# @transform (zegami_roi_flat2stacked, regex(r'zegamistack/zegami.yaml'), r'zegamistack/.ruffus')
# def zegami_roi_stacked_create (infile, outfile):
#     indir = os.path.dirname(infile)
#     base_dir = indir.replace('deepcell','')
#     outdir = os.path.dirname(outfile)
#     statement = '''zeg create collections --project %(zegami_collection)s 
#                     --config zegamistack/zegami.yaml  && touch %(outfile)s'''
#     P.run(statement, job_queue=PARAMS['batch_queue'])

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
