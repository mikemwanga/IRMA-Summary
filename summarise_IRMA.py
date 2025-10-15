#!/bin/python3

'''
SCRIPT TO SUMMARISE IRMA ASSEMBLER OUTPUT.
REQUIRES INSTALLATION OF THE BELOW LIBRARIES IN CONDA ENVIROMENT
RUN THE SCRIPT IN DIRECTORY WHERE IRMA OUTPUTS ARE LOCATED

USAGE python3 summarise_IRMA.py

'''
import config
from config import *
from figures import plot_figures_length, plot_figures_log, plot_assembled_reads, plot_line_coverage_depth

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
            cov_data = pd.read_table(os.path.join(tmp_folder,coverage_file))[['Reference_Name','Position','Consensus_Count']]
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
        coverage_data = pd.read_table(genefile, usecols=['Position','Consensus_Count'],keep_default_na=False)
        coverage_data['sample_name'] = folder
        # print(coverage_data)
        read_cov = pd.merge(reference_data, coverage_data, on='Position', how='left')
        read_cov['Consensus_Count'].fillna(0,inplace=True)
        merged_coverage_data.append(read_cov)

    merged_coverage_data = pd.concat(merged_coverage_data)
 
    merged_coverage_data.to_csv(f'{ref_folder_path}/{folder}-read-coverage-table.tsv',index=False,sep='\t')

    seq_data_frame = []
    for gene_reference_file in reference_files:
        gene_name = gene_reference_file.name.split('_')[1].split('.')[0]
        with open(gene_reference_file,'r') as file:
            sequence = file.readlines()[1].strip('\n')
            sequence_length=len(sequence)

    #EXTRACT SUBTYPE INFORMATION*************************************************************
    HSubtypes=[]
    NSubtypes=[]
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*_HA_*.fasta') and os.path.isfile(os.path.join(path, file)):
            H_name = file.split('_')[-1].split('.')[0]
            HSubtypes.append(H_name)
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*_NA_*.fasta') and os.path.isfile(os.path.join(path, file)):
            N_name = file.split('_')[-1].split('.')[0]
            NSubtypes.append(N_name)
    
    HSubtypes_join = ''.join(map(str, HSubtypes))
    NSubtypes_join = ''.join(map(str, NSubtypes))
    subtype = HSubtypes_join+NSubtypes_join #merge the subtypes into a single name

    #EXTRACT COVERAGE INFORMATION FOR EACH GENE SEGMENT********************************************

    gene_cov_data_file=[]
    for file in os.listdir(table_path):
        if file.endswith('coverage.txt'):
            data = pd.read_table(os.path.join(table_path,file))
            genelength = len(data) #SEQUENCE LENGTH IS RELATIVE TO DATA LENGTH
            name = file.split('-')[0]
            genename = name.split('_')[1]
            meanv= round(data['Consensus_Count'].mean(),1)
            df = {'sample':folder+'_' + subtype, 'Gene':genename,'seq_len':genelength,'coverage':meanv}
            datalist.append(df)
            #pass this to a file data and call a separate function to do the new plot type
            gene_name = file.split('_')[1].split('-')[0]
            # gene_name = file.split('-')[0]#.split('-')[0]
            color = color_map.get(gene_name, 'gray') 
            gene_cov_data_file.append((data[['Reference_Name','Position','Consensus_Count']], gene_name,color)) #>>>>>> UPDATED ON 1-09-2024

    file_data[folder] = gene_cov_data_file

    #*******************************************************************************************************************************
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
readsdata['assembled'] = round(readsdata['match']/readsdata['initial']*100,1)
readsdata.rename(columns = {'index':'sample', 'initial':'total_reads','match':'total_assembled'},inplace=True)
readsdata.reset_index(inplace=True)
readsdata.rename(columns = {'index':'sample'}, inplace=True)

#*******************************************************************************************************************************

#FUNCTION TO PLOT HEATMAP TO DEPTH

#PROCESS DATAFRAME EXTRACT SEQUENCE LENGTH AND COVERAGE
def arrange_duplicates(data):
    "append suffix in samplesIds in coinfections. adding suffixes. Ensure data has the included columns here."
    data['suffix'] = data.groupby(['sample', 'Gene']).cumcount().map(lambda x: '' if x == 0 else x+1)
    data['suffix'] = data['suffix'].apply(lambda x:'-'+ str(x) if pd.notna(x) and x != '' else x)
    data['sample'] = data['sample'] + data['suffix']
    data.drop(columns=['suffix'], inplace=True)
    return data

gene_lenth = dataframe[dataframe.columns.difference(['coverage'])]
gene_lenth = arrange_duplicates(gene_lenth)


gene_lenth= gene_lenth.pivot(index='sample',columns='Gene',values='seq_len').reset_index()
# print(gene_lenth)
columns = [col for col in geneids if col in gene_lenth.columns]
gene_lenth = gene_lenth[['sample'] + columns]

for name in gene_lenth['sample']:
    gene_lenth['sample_id'] = [name.split('_H')[0] for name in gene_lenth['sample']]
#Export raw length data
gene_lenth.to_csv(f'{irma_dir}/gene_absolute_length_coverage_data.tsv', sep="\t", index=True)

##**************************************************************************************************
# length_counts = {'HA':1701,'NA':1416,'MP':982,'PB1':2276,'PB2':2280,'NP':1497,'NS':838,'PA':2151} #PB1 initially et at 2274 but should be 2276
gene_lenth = gene_lenth.set_index('sample_id')

# print(gene_lenth)
reference_lengths_gene = dict()
for sample_name in gene_lenth['sample']:
    foldername = sample_name.split('_H')[0]
    # foldername=sample_name
    reference_path = Path(os.path.join(CURRDIR,foldername))/'intermediate/0-ITERATIVE-REFERENCES'
    # print(reference_path)
    reference_files = reference_path.glob('R0-*.ref')
    reference_lengths = {} #create dictionary to HOLD GENE NAME AND LENGTH
    for gene_reference_file in reference_files:
        # print(gene_reference_file)
        gene_name = gene_reference_file.name.split('_')[1].split('.')[0]
        # gene_name = gene_reference_file.name.split('-')[1]#.split('.')[0]
        with open(gene_reference_file,'r') as file:
            sequence = file.readlines()[1].strip('\n')
            sequence_length=len(sequence)
            # print(sequence_length)
            reference_lengths[gene_name] = sequence_length
    reference_lengths_gene[foldername]=reference_lengths
reference_lengths = pd.DataFrame(reference_lengths_gene).T
# print(reference_lengths)

for col in columns:
#     if col in gene_lenth.columns:
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


##**************************************************************************************************
#GENERATE DATA FOR COVERAGE BY DEPTH
gene_coverage = dataframe[dataframe.columns.difference(['seq_len'])]
gene_coverage = arrange_duplicates(gene_coverage)
df_coverage = gene_coverage.pivot(index='sample',columns='Gene',values='coverage').reset_index()

df_coverage = df_coverage.set_index('sample')
df_coverage = df_coverage[columns] #['PB1','PB2','PA','HA','NP','NA','NS','MP']]
df_coverage_log = round(np.log10(df_coverage[df_coverage.columns]+1), 1)
#EXPORT THE DATAFRAME
print('Output data to files')
df_coverage.to_csv(f'{irma_dir}/depth_coverage_data.tsv', sep="\t", index=True)
df_coverage_log.to_csv(f'{irma_dir}/depth_coverage_data_log_10.tsv', sep="\t", index=True)
gene_lenth.to_csv(f'{irma_dir}/gene-percentage-length-coverage-data.tsv', sep="\t", index=True)
readsdata.to_csv(f'{irma_dir}/assembled-reads-data.tsv', sep="\t", index=False)


#for every segment get counts of gene numbers
sample_list = []
for folder in os.listdir(CURRDIR):
    filepath=os.path.join(CURRDIR,folder,'tables','READ_COUNTS.txt')
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0: #if not folder.startswith('IRMA'):
        try:
            df = pd.read_table(filepath,usecols=['Record','Reads'])
            df['Record'] = df['Record'].apply(lambda value: value.split('-')[1])
            df = df[~df['Record'].isin(['initial','failQC','passQC','chimeric','nomatch','match','altmatch'])]
            df['Record'] = df['Record'].apply(lambda value: value.split('_')[1])#.dropna(subset='Reads',inplace=True)
            df.dropna(subset='Reads',inplace=True)
            df = df.sort_values(by='Reads').drop_duplicates(subset='Record',keep='last')
            df['sample'] = folder
            df = df.pivot(index='sample',values='Reads',columns='Record').reset_index()
            sample_list.append(df)
            #print(df)
        except Exception as e:
            print(f"failed to preocess {filepath}: {e}")
sample_data = pd.concat(sample_list).fillna(0)#.set_index('Sample')
# # print(sample_data
sample_data.to_csv(os.path.join(irma_dir,'gene_count.tsv'),sep='\t',index=False)


plot_line_coverage_depth(file_data,irma_dir, reference_folder)
print('Plotting figures for coverage by depth logwise')
plot_figures_log(df_coverage_log)
# print('Plotting figures for coverage by gene length')
plot_figures_length(gene_lenth)
# print('Generate Read Data')
plot_assembled_reads(readsdata)
print('Voila!')