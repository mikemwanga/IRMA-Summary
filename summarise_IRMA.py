#!/bin/python3

'''
SCRIPT TO SUMMARISE IRMA ASSEMBLER OUTPUT.
REQUIRES INSTALLATION OF THE BELOW LIBRARIES IN CONDA ENVIROMENT
RUN THE SCRIPT IN DIRECTORY WHERE IRMA OUTPUTS ARE LOCATED

USAGE python3 summarise_IRMA.py

'''
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import os,sys, argparse
# import glob
import plotly.express as px
import plotly.io as pio
import kaleido
import numpy as np
import math
import fnmatch
from pathlib import Path
# import re
import shutil
pio.kaleido.scope.mathjax = None
margin=dict(l=5, r=5, t=50, b=5) 
from datetime import datetime, date 
today = datetime.today().strftime('%Y-%m-%d')
time = datetime.today().strftime('%H-%M-%S')
time = today+'_' + time


CURRDIR=sys.argv[1]
# CURRDIR=os.getcwd()


#DECLARING SOME FUNCTIONS
# def plot_sample_coverage(data):
#     maxvalue = data['Coverage Depth'].max(numeric_only=True).max()
#     '''Plots gene segment coverage based on primary reference used'''
#     fig = px.line(data,x='Position',y='Coverage Depth',facet_col="Reference_Name",facet_col_wrap=4,
#                 facet_row_spacing = 0.2, facet_col_spacing = 0.05,
#                 range_y=(0,maxvalue ))
#     fig.update_yaxes(matches=None, showticklabels=True,title=None)
#     fig.update_xaxes(matches=None, showticklabels=True)
#     fig.update_layout(margin=margin)
#     return fig

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

#GET DIRECTORIES WITH COVERAGE DATA
selected_folder = []
for folder in os.listdir(CURRDIR):
    folder_path=os.path.join(CURRDIR, folder)
    if os.path.isdir(folder_path): #check whether folder exists
        table_folder = os.path.join(folder_path, 'tables')
        if os.path.exists(table_folder) and os.path.isdir(table_folder):
            table_files = [f for f in os.listdir(table_folder) if f.endswith('coverage.txt')]
            if table_files:
                selected_folder.append(folder)


#------------------------------------------------------------------------------------------
#GET AVERAGE READ DEPTH FOR EACH GENE SEGMENT BY SUMMING THE COVERAGE AND DEVIDE BY GENE LENGTH.
#CREATE A DATA FRAME WITH SAMPLE NAME, GENE NAME AND AVERAGE COVERAGE
datalist = []
readsdf = []


file_data = {} #COLLECT COVERAGE DATA, GENENAME AND COLOR FOR COVERAGE PLOT
color_map = {
                'HA': 'blue','NA': '#075d1c','MP': 'black','PA': '#ff21b2','PB2': '#5a189a','PB1': '#9e5231',
                'NP': 'red','NS': '#339dff',
            }

print('Get Subtype gene segment lengthe and average coverage')

for folder in selected_folder:
    #GET GENOTYPE INFORMATION
    H_name = 'Hx'
    N_name='Nx'
    path = os.path.join(CURRDIR,folder)
    table_path = os.path.join(CURRDIR, folder, 'tables') #GET PATH FOR THE COVERAGE DATAFRAME THE FOLDERS

    #EXTRACT LENGTH OF REFERENCE SEQUNCE USED IN FIRST ITERATION***************************
    reference_path = Path(path)/'intermediate/0-ITERATIVE-REFERENCES'
    reference_files = reference_path.glob('R0-*.ref')
    coverage_path = Path(table_path)
    coverage_files = coverage_path.glob('*coverage.txt')
    #create temp folder

    # print(file for file in coverage_files)

    if os.path.exists(os.path.join(CURRDIR,folder,'tmp')):
        shutil.rmtree(os.path.join(CURRDIR,folder,'tmp'))
        os.mkdir(os.path.join(CURRDIR,folder,'tmp'))
    else:
        os.mkdir(os.path.join(CURRDIR,folder,'tmp'))

    #COPY FILES TO tmp FOLDER
    tmp_folder = os.path.join(CURRDIR,folder,'tmp')
    genenames = []
    for filename in coverage_files:
        shutil.copy(filename, tmp_folder)
        # print(filename)
    for file in reference_files:
        shutil.copy(file, tmp_folder)

    for coverage_file in os.listdir(tmp_folder):
        if coverage_file.endswith('coverage.txt'):
            name = coverage_file.split('-')[0]
            genename = name.split('_')[1]
            genenames.append(genename)
            cov_data = pd.read_table(os.path.join(tmp_folder,coverage_file))[['Reference_Name','Position','Coverage Depth']]

            # print(cov_data)
            cov_data['Reference_Name'] = cov_data['Reference_Name'].apply(lambda name:name.split('_')[1])
            cov_data.to_csv(f'{tmp_folder}/{genename}.cov.txt',sep='\t',index=False)

        if coverage_file.endswith('.ref'):
            reference_file = Path(os.path.join(tmp_folder,coverage_file))
            gene_name = reference_file.name.split('_')[1].split('.')[0]

            seq_data_frame = []
            with open(reference_file,'r') as file:
                sequence = file.readlines()[1].strip('\n')
                sequence_length=len(sequence)
                for pos, nt in enumerate(sequence, start=1):
                    seq_data_frame.append({'Reference_Name':gene_name, 'Position':pos,'Nucleotide':nt})
            ref_data = pd.DataFrame(seq_data_frame)
            ref_data.to_csv(f'{tmp_folder}/{gene_name}.ref.txt',sep='\t',index=False)
    seqfiles = [seq_file for seq_file in os.listdir(path) if seq_file.endswith('.fasta')]

    merged_coverage_data = []
    for gene_segment in genenames:
        genefile = os.path.join(tmp_folder,f'{gene_segment}.ref.txt')
        reference_data = pd.read_table(genefile, keep_default_na=False)
        genefile = os.path.join(tmp_folder,f'{gene_segment}.cov.txt')
        # print(pd.read_table(genefile))
        coverage_data = pd.read_table(genefile, usecols=['Position','Coverage Depth'],keep_default_na=False)
        coverage_data['sample_name'] = folder
        # print(coverage_data)
        read_cov = pd.merge(reference_data, coverage_data, on='Position', how='left')
        read_cov['Coverage Depth'].fillna(0,inplace=True)
        merged_coverage_data.append(read_cov)

    merged_coverage_data = pd.concat(merged_coverage_data)
    # print(merged_coverage_data)
    # coverage_plot = plot_sample_coverage(merged_coverage_data)
    # coverage_plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    # coverage_plot.write_image(f'{ref_folder_path}/{folder}-read-coverage-ref.pdf')
    # merged_coverage_data.to_csv(f'{ref_folder_path}/{folder}-read-coverage-table.tsv',index=False,sep='\t')

    # reference_lengths = {} #create dictionary to HOLD GENE NAME AND LENGTH
    seq_data_frame = []
    for gene_reference_file in reference_files:
        gene_name = gene_reference_file.name.split('_')[1].split('.')[0]
        with open(gene_reference_file,'r') as file:
            sequence = file.readlines()[1].strip('\n')
            sequence_length=len(sequence)

    #EXTRACT SUBTYPE INFORMATION*************************************************************
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*_HA_*.bam') and os.path.isfile(os.path.join(path, file)):
            H_name = file.split('_')[-1].split('.')[0]
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*_NA_*.bam') and os.path.isfile(os.path.join(path, file)):
            N_name = file.split('_')[-1].split('.')[0]
    subtype = H_name+N_name #merge the subtypes into a single name

    #EXTRACT COVERAGE INFORMATION FOR EACH GENE SEGMENT********************************************
    # table_path = os.path.join(CURRDIR, folder, 'tables') #GET PATH FOR THE COVERAGE DATAFRAME THE FOLDERS
    # print(folder)
    gene_cov_data_file=[]
    for file in os.listdir(table_path):
        # print(file)
        # 
        if file.endswith('coverage.txt'):
            data = pd.read_table(os.path.join(table_path,file))
            genelength = len(data) #SEQUENCE LENGTH IS RELATIVE TO DATA LENGTH
            name = file.split('-')[0]
            genename = name.split('_')[1]
            meanv= round(data['Coverage Depth'].mean(),1)
            df = {'sample':folder+'_' + subtype, 'Gene':genename,'seq_len':genelength,'coverage':meanv}
            datalist.append(df)
            #pass this to a file data and call a separate function to do the new plot type
            gene_name = file.split('_')[1].split('-')[0]
            color = color_map.get(gene_name, 'gray') 
            gene_cov_data_file.append((data[['Reference_Name','Position','Coverage Depth']], gene_name,color)) #>>>>>> UPDATED ON 1-09-2024
    # print(gene_cov_data_file)
    file_data[folder] = gene_cov_data_file






    #**********************************************************************************************
    #GET SUMMARY OF READS FROM READS_COUNT.xt
    reads_data = pd.read_table(table_path + '/READ_COUNTS.txt')
    reads_data = reads_data[['Record','Reads']]
    reads_data['Record'] = reads_data['Record'].str.rsplit('_').str[-1].str.split('-').str[-1]
    reads_data.rename(columns={'Reads':folder+'_' + subtype}, inplace=True)
    reads_data = reads_data.set_index('Record')
    reads_data = reads_data.loc[['initial','match']] #select relevant rows
    readsdf.append(reads_data.T)
    
dataframe = pd.DataFrame(datalist)
#create a dataframe of total reads with assembled reads and pencentage assembled
readsdata = pd.concat(readsdf) #returns a data frame with total reads and assembled reads
readsdata['assembled(%)'] = round(readsdata['match']/readsdata['initial']*100,1)
readsdata.rename(columns = {'initial':'total_reads','match':'total_assembled'},inplace=True)
readsdata.reset_index(inplace=True)


#-----------------------------------------------------------------------------------------------------------------
#FUNCTION TO PLOT HEATMAP TO DEPTH

#PROCESS DATAFRAME EXTRACT SEQUENCE LENGTH AND COVERAGE
gene_lenth = dataframe[dataframe.columns.difference(['coverage'])]
gene_coverage = dataframe[dataframe.columns.difference(['seq_len'])]
gene_lenth= gene_lenth.pivot(index='sample',columns='Gene',values='seq_len').reset_index()
# print(gene_lenth)
geneids = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
columns = [col for col in geneids if col in gene_lenth.columns]
gene_lenth = gene_lenth[['sample'] + columns]


for name in gene_lenth['sample']:
    gene_lenth['sample_id'] = [name.split('_H')[0] for name in gene_lenth['sample']]
#Export raw length data
# gene_lenth.to_csv(f'{irma_dir}/gene_absolute_length_coverage_data.tsv', sep="\t", index=True)

##**************************************************************************************************
# length_counts = {'HA':1701,'NA':1416,'MP':982,'PB1':2276,'PB2':2280,'NP':1497,'NS':838,'PA':2151} #PB1 initially et at 2274 but should be 2276
gene_lenth = gene_lenth.set_index('sample_id')

# print(gene_lenth)
reference_lengths_gene = dict()
for sample_name in gene_lenth['sample']:
    foldername = sample_name.split('_H')[0]
    reference_path = Path(os.path.join(CURRDIR,foldername))/'intermediate/0-ITERATIVE-REFERENCES'
    reference_files = reference_path.glob('R0-*.ref')
    reference_lengths = {} #create dictionary to HOLD GENE NAME AND LENGTH
    for gene_reference_file in reference_files:
        gene_name = gene_reference_file.name.split('_')[1].split('.')[0]
        with open(gene_reference_file,'r') as file:
            sequence = file.readlines()[1].strip('\n')
            sequence_length=len(sequence)
            reference_lengths[gene_name] = sequence_length
    reference_lengths_gene[foldername]=reference_lengths
reference_lengths = pd.DataFrame(reference_lengths_gene).T
# print(reference_lengths)

for col in columns:
    # if col in gene_lenth.columns:
    gene_lenth[col] = round(gene_lenth[col]/reference_lengths[col] *100,1)
gene_lenth.reset_index(inplace=True)

#CREATE A SUBTYPING FILE
type_dict = list()
for sample in gene_lenth['sample']:
    sampleid = sample.split('_')[0]
    stype = sample.split('_')[1]
    type_dict.append({'SampleID':sampleid,'Subtype':stype})

pd.DataFrame(type_dict).to_csv(f'{irma_dir}/influenza_subtypes.tsv', sep="\t", index=True)

gene_lenth=gene_lenth[['sample']+ columns].set_index('sample') #WHAT IF SOME OF THE SEGMENTS ARE MISSING. NEED TO CONDITION THIS

##**************************************************************************************************
#NEW PRINT FOR FILE DATA. CONSIDER RE-EDITING THE CODES

from line_plots import plot_line_coverage_depth
# print(file_data)
plot_line_coverage_depth(file_data,irma_dir, reference_folder)
# figg.show()


##**************************************************************************************************
#GENERATE DATA FOR COVERAGE BY DEPTH
df_coverage = gene_coverage.pivot(index='sample',columns='Gene',values='coverage').reset_index()
df_coverage = df_coverage.set_index('sample')
df_coverage = df_coverage[columns] #['PB1','PB2','PA','HA','NP','NA','MP','NS']]
df_coverage_log = round(np.log10(df_coverage[df_coverage.columns]+1), 1)
#EXPORT THE DATAFRAME
print('Output data to files')
# df_coverage.to_csv(f'{irma_dir}/depth_coverage_data.tsv', sep="\t", index=True)
# df_coverage_log.to_csv(f'{irma_dir}/depth_coverage_data_log_10.tsv', sep="\t", index=True)
# gene_lenth.to_csv(f'{irma_dir}/gene-percentage-length-coverage-data.tsv', sep="\t", index=True)
# readsdata.to_csv(f'{irma_dir}/assembled-reads-data.tsv', sep="\t", index=True)

#PLOT FUNCTIONS
'GnBu','OrRd'
font_color='black'
def plot_coverage_depth(data, header,bartitle,max_value,color_scale):
    '''Function to plot coverage by depth'''
    fig = px.imshow(data,text_auto=True,aspect='auto',color_continuous_scale=color_scale,zmin=0, zmax=max_value)
    fig.update_layout(margin = margin,title = header,font_color=font_color,
                      plot_bgcolor = "white",#height=300,
                      coloraxis_colorbar=dict(title=bartitle,orientation='v',#y=-0.3,
                                              titlefont=dict(size=12),
                                              titleside='top',thickness=13))
    fig.update_yaxes(tickfont=dict(size=14), title = None, linecolor='black',
                     showticklabels=True,ticks='outside',tickcolor='black')
    fig.update_xaxes(tickfont=dict(size=14),title = 'Gene_Name (Ordered from longest to shortest)', linecolor='black',
                     showticklabels=True,ticks='outside',tickcolor='black')
    fig.update_traces(ygap=0.5,xgap=0.5)
    return fig


#*************************************
print('Plotting figures for coverage by depth')

def plot_figures_log(data):
    max_value = data.max(numeric_only=True).max()
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset_{i+1}'
        fig = plot_coverage_depth(subset_df, f'Read Coverage - Log Counts_{name}','Log10', max_value,'OrRd')
        fig.write_image(f'{irma_dir}/coverage-depth-log-counts-{name}.pdf')

def plot_figures_length(data):
    max_value = data.max(numeric_only=True).max()
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset-{i+1}'
        fig = plot_coverage_depth(subset_df, f'Gene Coverage - Length_{name}','Coverage (%)',max_value,'GnBu' )
        fig.write_image(f'{irma_dir}/coverage_length_{name}.pdf')


#**************************************************************************************************************
#Plot total reads and assembled 
def plot_bar(data,header):
    '''Function to plot reads assembled relative to total number of reads'''
    fig = px.bar(data, y='index',x = 'assembled(%)', text="total_reads", range_x=(0,100))
    fig.update_layout(margin = margin,title = header,plot_bgcolor = "white", font_color=font_color)
    fig.update_yaxes(titlefont=dict(size=14),tickfont=dict(size=14),tickcolor='black',title=None,
                     linecolor='black',showticklabels=True,ticks='outside')
    fig.update_xaxes(linecolor='black',titlefont=dict(size=14),tickfont=dict(size=14),tickcolor='black',
                     title = '% Reads Assembled',gridcolor='#ddd9d9',
                     showticklabels=True,ticks='outside')
    fig.update_traces(marker_color = '#367588',width=0.8,
                      textfont_size=12, textangle=0, textposition="inside", cliponaxis=False)
    return fig

def plot_assembled_reads(data):
    data.sort_values('index', inplace=True)
    binsize=20
    num_bins = (len(data) + binsize -1) // binsize #// RETURNS A ROUNDED OFF VALUE AND NOT FLOAT
    for i in range(num_bins):
        start = i * binsize
        stop = min((i + 1) * binsize,len(data))
        subset_df = data.iloc[start:stop, :]
        name =f'Subset-{i+1}'
        fig = plot_bar(subset_df, f'{name}_total-reads_vs_%-assembled' )
        fig.write_image(f'{irma_dir}/assembled_reads_{name}.pdf')

print('Plotting figures for coverage by depth logwise')
plot_figures_log(df_coverage_log)
print('Plotting figures for coverage by gene length')
plot_figures_length(gene_lenth)
print('Generate Read Data')
plot_assembled_reads(readsdata)

print('Voila!')