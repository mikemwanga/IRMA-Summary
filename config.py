#!/bin/python3

# ALL PACKAGES IMPORTATION HERE

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import argparse
import fnmatch
import glob
import math
import shutil
from pathlib import Path
from datetime import datetime, date

import kaleido
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

# MERGE CHANGE: prevent MathJax lookup issues in static image export.
pio.kaleido.scope.mathjax = None

today = datetime.today().strftime('%Y-%m-%d')
time = datetime.today().strftime('%H-%M-%S')
time = today + '_' + time

# MERGE CHANGE: shared plot style constants (from irma_dev, cleaned).
black_color = 'rgb(0,0,0)'
textcolor = '#212020'
gridcolor = '#ccd1d1'
tickfontcolor = '#ccd1d1'
fontfamily = 'arial'
barcolor = '#afb2b3'
bgcolor = 'rgba(0, 0, 0, 0)'
width = 3
grid_color = '#c5c3c3'
margin = dict(l=5, r=5, t=50, b=5)
font_color = '#191818'
textsize = 10

# DEFINING GENE IDS
geneids = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']

# DEFINE COLOR FOR THE GENE SEGMENTS MAPPING
color_map = {
    'HA': 'blue',
    'NA': '#075d1c',
    'MP': 'black',
    'PA': '#ff21b2',
    'PB2': '#5a189a',
    'PB1': '#9e5231',
    'NP': 'red',
    'NS': '#339dff',
}

# MERGE CHANGE: color map for assembled/not assembled read distribution plots.
gray = '#d1d1d1'
red = '#888686'
read_colormap = {
    'total_assembled': red,
    'notassembled': gray,
    'prop_assembled': red,
    'prop_notassembled': gray,
}

CURRDIR = sys.argv[1]

# CREATE OUTPUT DIRECTORY
dir_name = 'IRMA_Summary'
reference_folder = 'Read_Depth_Plot'
dirname = dir_name + '_' + time
if not os.path.exists(os.path.join(CURRDIR, dirname)):
    print('Directory creating')
    os.mkdir(os.path.join(CURRDIR, dirname))
    os.mkdir(os.path.join(CURRDIR, dirname, reference_folder))
else:
    print('Directory Exists')
    print('Generate summary plots')

ref_folder_path = os.path.join(CURRDIR, dirname, reference_folder)
irma_dir = os.path.join(CURRDIR, dirname)
