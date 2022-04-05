"""Converts all OME-TIFF files in input folder to a folder compatible with HistoCAT software.
Parameters
----------
input_folder
    Input folder.
output_folder
    Output folder.
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

from imctools.io.ometiff.ometiffparser import OmeTiffParser
from imctools.converters import ome2histocat


logger = logging.getLogger(__name__)

def omefolder_to_histocatfolder(
    input_folder: Union[str, Path],
    output_folder: Union[str, Path],
):
    if isinstance(input_folder, Path):
        input_folder = str(input_folder)

    ome_files = [os.path.basename(fn) for fn in glob.glob(os.path.join(input_folder, "*")) if fn.endswith(".ome.tiff")]

    for fn_ome in ome_files:
        len_suffix = len(".ome.tiff")
        basename_ome = fn_ome[:-len_suffix]
        mask_file = None
        dtype = None
        path_ome = os.path.join(input_folder, fn_ome)
        
        ome2histocat.omefile_to_histocatfolder(path_ome, output_folder, mask_file=mask_file, dtype=dtype)

        ruff = str(output_folder) + "/" + str(basename_ome) + "/.ruffus"
        Path(ruff).touch()


if __name__ == "__main__":
    import timeit

    tic = timeit.default_timer()
    
    # Parse command line arguments
    parser = argparse.ArgumentParser() 
    parser.add_argument('-i', '--input', dest='input_folder', 
                        help='Input folder name', default='ometiff')
    parser.add_argument('-o', '--output_folder', dest='output_folder', 
                        help='output folder name', default = 'histocat')
    args = parser.parse_args()
    print(args)

    omefolder_to_histocatfolder(
        Path(args.input_folder),
        Path(args.output_folder),
    )
    

    print(timeit.default_timer() - tic)
