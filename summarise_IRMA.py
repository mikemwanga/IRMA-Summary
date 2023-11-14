#!/bin/python3
import pandas as pd
import os
import plotly.express as px
import kaleido
import numpy as np

margin=dict(l=5, r=5, t=50, b=5) 


CURRDIR=os.getcwd()

#GET DIRECTORIES WITH NECESSARY DATA
selected_folder = []
for folder in os.listdir(CURRDIR):
    folder_path=os.path.join(CURRDIR, folder)
    if os.path.isdir(folder_path):
        table_folder = os.path.join(folder_path, 'tables')
        if os.path.exists(table_folder) and os.path.isdir(table_folder):
            table_files = [f for f in os.listdir(table_folder) if f.endswith('coverage.txt')]
            if table_files:
                selected_folder.append(folder)

#GET DATAFRAME FROM THE COVERAGE PLOT
datalist = []

for folder in selected_folder:
    table_path = os.path.join(CURRDIR, folder, 'tables')
    # print(table_path)
    for file in os.listdir(table_path):
        # print(file)
        if file.endswith('coverage.txt'):
            data = pd.read_table(os.path.join(table_path,file))
            genelength = len(data)
            name = file.split('-')[0]
            genename = name.split('_')[1]
            maxv = data['Coverage Depth'].max()
            minv = data['Coverage Depth'].min()
            meanv= round(data['Coverage Depth'].mean(),1)
            df = {'sample':folder, 'Gene':genename,'seq_len':genelength,'coverage':meanv}
            datalist.append(df)

dataframe = pd.DataFrame(datalist)

#PROCESS DATAFRAME EXTRACT SEQUENCE LENGTH AND COVERAGE
df_lenth = dataframe[dataframe.columns.difference(['coverage'])]
df_coverage = dataframe[dataframe.columns.difference(['seq_len'])]
df_lenth= df_lenth.pivot(index='sample',columns='Gene',values='seq_len').reset_index()
columns = ['HA','NA','MP','PB1','PB2','NP','NS','PA']
df_lenth = df_lenth[['sample'] + columns]
length_counts = {'HA':1701,'NA':1410,'MP':982,'PB1':2274,'PB2':2280,'NP':1497,'NS':838,'PA':2151}
for col in columns:
    df_lenth[col] = round(100*df_lenth[col]/length_counts[col],1)
df_lenth.set_index('sample', inplace=True)


df_coverage = df_coverage.pivot(index='sample',columns='Gene',values='coverage').reset_index()
df_coverage = df_coverage.set_index('sample')

df_coverage_log = round(np.log(df_coverage[df_coverage.columns]+1), 1)
# dataframe.set_index('sample', inplace=True)

# print(df_coverage_log)
#FUNCTION TO PLOT HEATMAP TO DEPTH
dirname = 'IRMA_Summary'
if not os.path.exists(os.path.join(CURRDIR,dirname)):
    print('Directory creating')
    os.mkdir(os.path.join(CURRDIR,dirname))
else:
    print('Directory Exists')
    print('Generate summary plots')
    ...

def plot_coverage(data, header):
    '''Function to plot coverage figure'''
    fig = px.imshow(data,text_auto=True,aspect='auto',color_continuous_scale='OrRd')
    fig.update_layout(margin = margin,title = header)
    fig.update_yaxes(tickfont=dict(size=10), title = None, linecolor='gray')
    fig.update_xaxes(title = None, linecolor='gray')
    return fig


depth_coverage_plot = plot_coverage(df_coverage, 'Read Coverage - Absolute Counts' )
depth_coverage_plot.write_image('IRMA_Summary/coverage_depth_absolute_counts.pdf')

depth_coverage_plot = plot_coverage(df_coverage_log, 'Read Coverage - Log Counts' )
depth_coverage_plot.write_image('IRMA_Summary/coverage_depth_log_counts.pdf')


depth_length_plot = plot_coverage(df_lenth, 'Coverage - Length')
depth_length_plot.write_image('IRMA_Summary/coverage_length_plot.pdf')

