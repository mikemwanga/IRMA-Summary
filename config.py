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
black_color = "#060606"
textcolor = black_color
gridcolor = '#ccd1d1'
tickfontcolor = '#ccd1d1'
fontfamily = 'arial'
barcolor = '#afb2b3'
bgcolor = 'rgba(0, 0, 0, 0)'
width = 3
margin = dict(l=5, r=5, t=50, b=5)
textsize = 16

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
figures_folder = 'figures'
# CHANGE: dedicated folder for per-sample depth PDFs to reduce crowding in top-level figures folder.
individual_depth_folder = 'individual_sample_depth'
summary_folder = 'summary_data'
dirname = dir_name + '_' + time
if not os.path.exists(os.path.join(CURRDIR, dirname)):
    print('Directory creating')
    os.mkdir(os.path.join(CURRDIR, dirname))
    os.mkdir(os.path.join(CURRDIR, dirname, figures_folder))
    os.mkdir(os.path.join(CURRDIR, dirname, figures_folder, individual_depth_folder))
    os.mkdir(os.path.join(CURRDIR, dirname, summary_folder))
else:
    print('Directory Exists')
    print('Generate summary plots')

irma_dir = os.path.join(CURRDIR, dirname)
figures_dir = os.path.join(CURRDIR, dirname, figures_folder)
# CHANGE: central path used by plot_line_coverage_depth() for individual sample depth plot exports.
individual_depth_dir = os.path.join(CURRDIR, dirname, figures_folder, individual_depth_folder)
summary_dir = os.path.join(CURRDIR, dirname, summary_folder)

# CHANGE: guarantee subfolder exists when reusing an existing run directory.
os.makedirs(individual_depth_dir, exist_ok=True)
