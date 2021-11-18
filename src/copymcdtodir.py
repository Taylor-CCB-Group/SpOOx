#!/usr/bin/env python3
# -*- coding: utf-8 -*-
DESC = ("Makes a directory for the MCD files in the output (usually the 'mcd') directory because mcdtools seems to require it when running to create ometiff files.")

from shutil import copy
import os
import sys
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument(
     '-i', '--input', 
        help = 'Input mcd file', 
        required = True
)

parser.add_argument(
     '-d', '--dir', 
        help = 'Directory to write the dir/filename into', 
        required = True
)


args = parser.parse_args()
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)		   



pathToMCD = args.input
pathToFolder = args.dir
subDir = os.path.basename(os.path.splitext(pathToMCD)[0])
fullPathToDir = os.path.join(pathToFolder,subDir)
print("Making directory:",fullPathToDir)
print("Copying directory:",pathToMCD,"to",fullPathToDir)


path = Path(fullPathToDir)
path.mkdir(parents=True, exist_ok=True)
copy(pathToMCD, fullPathToDir)

