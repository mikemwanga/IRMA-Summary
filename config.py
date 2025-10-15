#!/bin/python3

#ALL PACKAGES IMPORTATION HERE

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import os,sys, argparse
# import glob
import kaleido
import plotly.express as px
import plotly.graph_objects as go
# import plotly.io as pio

import numpy as np
import math
import fnmatch
from pathlib import Path

import shutil
#pio.kaleido.scope.mathjax = None

from datetime import datetime, date 
today = datetime.today().strftime('%Y-%m-%d')
time = datetime.today().strftime('%H-%M-%S')
time = today+'_' + time

#DEFINING GENE IDS
geneids = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
#DEFINE COLOR FOR THE GENE SEGMENTS MAPPING
color_map = {'HA': 'blue','NA': '#075d1c','MP': 'black','PA': '#ff21b2','PB2': '#5a189a',
            'PB1': '#9e5231','NP': 'red','NS': '#339dff'}

'GnBu','OrRd'

CURRDIR=sys.argv[1]

#CREATE OUTPUT DIRECTORY
dir_name = 'IRMA_Summary'
reference_folder='Read_Depth_Plot'
dirname = dir_name +'_'+ time
if not os.path.exists(os.path.join(CURRDIR,dirname)):
    print('Directory creating')
    os.mkdir(os.path.join(CURRDIR,dirname))
    os.mkdir(os.path.join(CURRDIR,dirname,reference_folder))
else:
    print('Directory Exists')
    print('Generate summary plots')

ref_folder_path = os.path.join(CURRDIR,dirname,reference_folder)
irma_dir = os.path.join(CURRDIR,dirname)

width = 3
grid_color = '#ddd9d9'
margin=dict(l=5, r=5, t=50, b=5) 
font_color='black'