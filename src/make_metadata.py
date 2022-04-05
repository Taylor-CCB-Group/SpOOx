#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re

cwd = os.getcwd()

filedata = "sample_id	sample_name	condition	ROI	path\n"
for d in os.listdir('histocat'):
    pattern = '_ROI_'
    if (re.search(pattern, d)):
        (sample_name, ROI) = re.split(pattern, d)
        full_path = cwd + '/' + 'signalextraction' + '/' + d
        metadata = d + "\t" + sample_name + "\t-\tROI_" + ROI + "\t" + full_path + "\n"
        filedata += metadata
        
metadata_file = cwd + '/' + 'metadata.tsv'
f = open(metadata_file,'w')  
f.write(filedata)
f.close()