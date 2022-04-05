"""
Converts folder (or zipped folder) containing raw acquisition data (mcd and txt files) 
to IMC folder containing standardized files.

Parameters
----------
input
    Input folder / .zip file with  raw .mcd/.txt acquisition data files.
output_folder
    Path to the output folder.
create_zip
    Whether to create an output as .zip file.
parse_txt
    Always use TXT files if present to get acquisition image data.
"""

import glob
import re
import os
import logging
import zipfile
import timeit
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory 
from typing import Union
from imctools.io.imc.imcwriter import ImcWriter
from imctools.io.mcd.mcdparser import McdParser
from imctools.io.txt.txtparser import TXT_FILE_EXTENSION, TxtParser
from imctools.io.utils import MCD_FILENDING, SCHEMA_FILENDING, ZIP_FILENDING

logger = logging.getLogger(__name__)

def mcdfolder_to_imcfolder(
    input: Union[str, Path], 
    output_folder: Union[str, Path], 
    create_zip: bool = False,
    parse_txt: bool = False):

    # convert input string to path object 
    if isinstance(input, str):
        input = Path(input)
    # if input is zipped extract into temp directory
    tmpdir = None
    if input.is_file() and input.suffix == ZIP_FILENDING:
        tmpdir = TemporaryDirectory()
        with zipfile.ZipFile(input, allowZip64=True) as zip:
            zip.extractall(tmpdir.name)
        input_folder = Path(tmpdir.name)
    else:
        input_folder = input

    # Parse mcd file - expect only one per folder
    # Exclude filenames starting with .
    mcd_parser = None
    try:
        mcd_files = list(input_folder.rglob(f"*{MCD_FILENDING}"))
        mcd_files = [f for f in mcd_files if not f.name.startswith(".")]
        assert len(mcd_files) == 1
        input_folder = mcd_files[0].parent
        # Parse schema files in mcd directory, expect 0 or 1
        schema_files = glob.glob(str(input_folder / f"*{SCHEMA_FILENDING}"))
        schema_file = schema_files[0] if len(schema_files) > 0 else None
        try:
            mcd_parser = McdParser(mcd_files[0])
        except:
            if schema_file is not None:
                logging.error("MCD file is corrupted, trying to rescue with schema file")
                mcd_parser = McdParser(mcd_files[0], xml_metadata_filepath=schema_file)
            else:
                raise

        # Parse text files in mcd directory
        txt_files = glob.glob(str(input_folder / f"*[0-9]{TXT_FILE_EXTENSION}"))
        txt_acquisitions_map = {TxtParser.extract_acquisition_id(f): f for f in txt_files}

        # Change session name as this is used as basis of output file names
        #   note - use the name of the input mcd file to keep naming schemes consistent
        #          this step was used previously to remove spaces in session.metaname - now changing session.metaname entirely)
        pat = 'mcd\/.*\/(.*)\.mcd'
        match = re.search(pat, mcd_parser.mcd_filename)
        mcd_file_root = match.group(1)
        setattr(mcd_parser.session, 'name', mcd_file_root)
        print("NAME:"+str(mcd_parser.session.metaname)+"OUT:"+str(parse_txt))

        # Write output files
        imc_writer = ImcWriter(output_folder, mcd_parser, txt_acquisitions_map, parse_txt)
        imc_writer.write_imc_folder(create_zip=create_zip)
        
        # rename tiff files:
        out_files_dir = str(output_folder) + '/' + str(mcd_parser.session.metaname) + "/*.tiff"
        #out_files_dir = str(output_folder) + "/*.tiff"
        tiff_files = glob.glob(str(out_files_dir))
        for tf in tiff_files:
            pat = '_s0_a(\d+)_ac\.ome\.tiff'
            match = re.search(pat, tf)
            roi = int(match.group(1))
            new_tiff_filename = str(output_folder) + '/' + str(mcd_parser.session.metaname) + '/' + str(mcd_parser.session.metaname) + '_ROI_' + str(roi) + '.ome.tiff'
            #new_tiff_filename = str(output_folder) +  '/' + str(mcd_parser.session.metaname) + '_ROI_' + str(roi) + '.ome.tiff'
            logging.info("renaming tiff file, from: "+tf+" to: "+new_tiff_filename)
            os.rename(tf, new_tiff_filename)
        
        
    finally:
        if mcd_parser is not None:
            mcd_parser.close()
        if tmpdir is not None:
            tmpdir.cleanup()


if __name__ == "__main__":
    
    tic = timeit.default_timer()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser() 
    parser.add_argument('-o', '--output_folder', dest='output_folder', 
                        help='output folder name', default = 'ometiff')
    parser.add_argument('-i', '--input', dest='input_folder', 
                        help='Input folder name', default='.')
    parser.add_argument('-z', '--create_zip', dest='create_zip', 
                        help='Whether to create an output as .zip file', 
                        action='store_true', default=False)
    parser.add_argument('-t', '--parse_txt', dest='parse_txt', 
                        help='Always use TXT files if present to get acquisition image data', 
                        action='store_true', default=False)
    args = parser.parse_args()
    print(args)

    mcdfolder_to_imcfolder(
        Path(args.input_folder),
        Path(args.output_folder),
        create_zip=args.create_zip,
        parse_txt=args.parse_txt,
    )

    print(timeit.default_timer() - tic)