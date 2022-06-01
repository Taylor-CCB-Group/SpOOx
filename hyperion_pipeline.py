'''
Spatial transcriptomics pipeline for Hyperion data

imaging mass cytometry from fluidigm
up to 40 different marker antibodies with metal tags

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


# import pipeline parameters from yaml format configuration file in working directory
PARAMS = P.get_parameters("pipeline.yml")

# imctools works on folders not files but folders are poor for tracking as their modicfication time is updated on access
# To mitigate this we track a hidden dummy file (.ruffus) added to each directory of interest for tracking purposes
@originate(PARAMS["hyperion_dir"]+"/*/.ruffus")
def mark_input_folders(outfiles):
    statement = '''find %(hyperion_dir)s/* -maxdepth 1 -type d -exec touch {}/.ruffus \;;'''
    P.run(statement, without_cluster=True)

# mcd_to_tiff
@follows (mkdir ("ometiff"),mark_input_folders)
@transform (PARAMS["hyperion_dir"]+"/*/.ruffus", regex(PARAMS["hyperion_dir"]+r'/(.*)/.ruffus'), r'ometiff/\1/.ruffus')
def mcd_to_tiff (infile, outfile):
    indir = os.path.dirname(infile)
    statement = '''python %(scripts_dir)s/parse_mcd.py -i %(indir)s -o ometiff && touch %(outfile)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# tiff_to_histocat
@follows (mkdir ("histocat"))
@subdivide(mcd_to_tiff, regex(r'ometiff/(.*)/.ruffus'), r'histocat/\1_ROI_*/.ruffus')
def tiff_to_histocat (infile, outfiles):
    indir = os.path.dirname(infile)
    statement = '''python %(scripts_dir)s/parse_ome.py -i %(indir)s -o histocat'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# make_config
@follows(tiff_to_histocat)
@originate(PARAMS["marker_file"])
def make_config(outfile):
    statement = '''python %(scripts_dir)s/writeconfig.py --indir histocat --outfile %(marker_file)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# deepcell
@transform (tiff_to_histocat, regex(r'histocat/(.*)/.ruffus'), r'deepcell/\1/deepcell.tif')
def deepcell (infile, outfile):
    indir = os.path.dirname(infile)
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/deepercell.py --markerfile %(marker_file)s 
                --indir %(indir)s --outdir %(outdir)s --contrast 5'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# signal_extraction
@transform (deepcell, regex(r'deepcell/(.*)/deepcell.tif'), r'signalextraction/\1/cellData.tab')
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
                %(signal_extraction_options)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# make_metadata
@follows(signal_extraction)
@originate("metadata.tsv")
def make_metadata(outfile):
    statement = '''python %(scripts_dir)s/make_metadata.py'''
    P.run(statement, without_cluster=True)

# phenoclustering
@follows(make_metadata)
@originate("clustering/phenocluster_sceobj_k30.RData")
def phenocluster(outfile):
    statement = '''Rscript Rphenoclustering.R
                --panel_file %(marker_file)s
                --metadata_file metadata.tsv
                --analysisName phenocluster
                --out_dir clustering
                %(phenograph_options)s'''
    P.run(statement, without_cluster=True, job_condaenv="hyperion_R")



### Zegami ######################################################################


# mark_histocat
@follows(tiff_to_histocat)
@originate("histocat/.ruffus")
def mark_histocat(outfile):
    statement = '''touch %(outfile)s'''
    P.run(statement, without_cluster=True)

# zegami_roi
@follows(tiff_to_histocat, mkdir("zegami"))
@transform (mark_histocat, regex(r'histocat/.ruffus'), r'zegami/zegami.yaml')
def zegami_roi (infile, outfile):
    statement = '''python %(scripts_dir)s/roi2zegami.py 
                    --indir histocat --outdir zegami --yaml zegami/zegami.yaml --name %(zegami_project)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

# zegami_roi_stacked
@follows(mkdir("zegamistack"))
@transform (deepcell, regex(r'deepcell/(.*)/deepcell.tif'), r'zegamistack/zegami.yaml')
def zegami_roi_stacked (infile, outfile):
    indir = os.path.dirname(infile)
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/deepercellcopyimages.py 
                    --zegamidir zegamistacks --indir %(indir)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])
    


### other Zegami-related code ######################################################################

@transform (zegami_roi, regex(r'zegami/zegami.yaml'), r'zegami/.ruffus')
def zegami_create (infile, outfile):
    statement = '''zeg create collections --project %(zegami_collection)s --config zegami/zegami.yaml
                    && touch zegami/.ruffus'''
    P.run(statement, job_queue=PARAMS['batch_queue'])


@transform (signal_extraction, 
            regex(r'signalextraction/(.*)/cellData.tab'), 
            r'signalextraction/\1/.zegami_upload')
def zegami_upload (infile, outfile):
    indir = os.path.dirname(infile)
    statement = '''python %(scripts_dir)s/uploadcellstozegami.py %(indir)s %(zegami_collection)s --verbose 
                    && touch %(outfile)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])


@follows(zegami_roi)
@transform (zegami_roi_stacked, regex(r'deepcell/(.*)/.ruffus'), r'zegamistack/zegami.yaml')
def zegami_roi_flat2stacked (infile, outfile):
    indir = os.path.dirname(infile)
    base_dir = indir.replace('deepcell','')
    outdir = os.path.dirname(outfile)
    statement = '''python %(scripts_dir)s/flat2stacks.py 
                    -i zegami/img/ -o zegamistack -y zegamistack/zegami.yaml 
                    -n %(zegami_project)s -d %(zegami_project)s 
                    -m %(marker_file)s -t zegamistack/zegami.tsv'''
    P.run(statement, job_queue=PARAMS['batch_queue'])


@transform (zegami_roi_flat2stacked, regex(r'zegamistack/zegami.yaml'), r'zegamistack/.ruffus')
def zegami_roi_stacked_create (infile, outfile):
    indir = os.path.dirname(infile)
    base_dir = indir.replace('deepcell','')
    outdir = os.path.dirname(outfile)
    statement = '''zeg create collections --project %(zegami_collection)s 
                    --config zegamistack/zegami.yaml  && touch %(outfile)s'''
    P.run(statement, job_queue=PARAMS['batch_queue'])

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )