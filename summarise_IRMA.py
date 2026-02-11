#!/bin/python3

'''
SCRIPT TO SUMMARISE IRMA ASSEMBLER OUTPUT.
REQUIRES INSTALLATION OF THE BELOW LIBRARIES IN CONDA ENVIROMENT
RUN THE SCRIPT IN DIRECTORY WHERE IRMA OUTPUTS ARE LOCATED

USAGE python3 summarise_IRMA.py
'''

import config
from config import *
from figures import (
    export_data,
    plot_all_line_plots,
    plot_assembled_reads,
    plot_figures_length,
    plot_figures_log,
    plot_line_coverage_depth,
    plot_read_distribution,
)


# GET DIRECTORIES WITH COVERAGE DATA
selected_folder = []
for folder in os.listdir(CURRDIR):
    folder_path = os.path.join(CURRDIR, folder)
    if os.path.isdir(folder_path):
        table_folder = os.path.join(folder_path, 'tables')
        if os.path.exists(table_folder) and os.path.isdir(table_folder):
            table_files = [f for f in os.listdir(table_folder) if f.endswith('coverage.txt')]
            if table_files:
                selected_folder.append(folder)


# GET AVERAGE READ DEPTH FOR EACH GENE SEGMENT BY SUMMING THE COVERAGE AND DIVIDE BY GENE LENGTH.
datalist = []
readsdf = []
file_data = {}

print('Get subtype, gene segment length and average coverage')

for folder in selected_folder:
    path = os.path.join(CURRDIR, folder)
    table_path = os.path.join(CURRDIR, folder, 'tables')

    # EXTRACT LENGTH OF REFERENCE SEQUENCE USED IN FIRST ITERATION
    reference_path = Path(path) / 'intermediate/0-ITERATIVE-REFERENCES'
    reference_files = list(reference_path.glob('R0-*.ref'))
    coverage_files = list(Path(table_path).glob('*coverage.txt'))

    # Build temp folder
    tmp_folder = os.path.join(CURRDIR, folder, 'tmp')
    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)
    os.mkdir(tmp_folder)

    # COPY FILES TO tmp FOLDER
    genenames = []
    for filename in coverage_files:
        shutil.copy(filename, tmp_folder)
    for filename in reference_files:
        shutil.copy(filename, tmp_folder)

    for filename in os.listdir(tmp_folder):
        fullpath = os.path.join(tmp_folder, filename)
        if filename.endswith('coverage.txt'):
            name = filename.split('-')[0]
            genename = name.split('_')[1]
            genenames.append(genename)
            cov_data = pd.read_table(fullpath)[['Reference_Name', 'Position', 'Consensus_Count']]
            cov_data['Reference_Name'] = cov_data['Reference_Name'].apply(lambda value: value.split('_')[1])
            cov_data.to_csv(f'{tmp_folder}/{genename}.cov.txt', sep='\t', index=False)

        if filename.endswith('.ref'):
            reference_file = Path(fullpath)
            gene_name = reference_file.name.split('_')[1].split('.')[0]
            seq_data_frame = []
            with open(reference_file, 'r') as handle:
                sequence = handle.readlines()[1].strip('\n')
                for pos, nt in enumerate(sequence, start=1):
                    seq_data_frame.append({'Reference_Name': gene_name, 'Position': pos, 'Nucleotide': nt})
            ref_data = pd.DataFrame(seq_data_frame)
            ref_data.to_csv(f'{tmp_folder}/{gene_name}.ref.txt', sep='\t', index=False)

    merged_coverage_data = []
    for gene_segment in genenames:
        ref_file = os.path.join(tmp_folder, f'{gene_segment}.ref.txt')
        cov_file = os.path.join(tmp_folder, f'{gene_segment}.cov.txt')
        if not os.path.exists(ref_file) or not os.path.exists(cov_file):
            print(f'Coverage/reference file missing for {folder}:{gene_segment}. Skipping.')
            continue

        reference_data = pd.read_table(ref_file, keep_default_na=False)
        coverage_data = pd.read_table(cov_file, usecols=['Position', 'Consensus_Count'], keep_default_na=False)
        coverage_data['sample_name'] = folder
        read_cov = pd.merge(reference_data, coverage_data, on='Position', how='left')
        read_cov['Consensus_Count'].fillna(0, inplace=True)
        merged_coverage_data.append(read_cov)

    if merged_coverage_data:
        pd.concat(merged_coverage_data).to_csv(
            f'{ref_folder_path}/{folder}-read-coverage-table.tsv', index=False, sep='\t'
        )

    # EXTRACT SUBTYPE INFORMATION
    ha_subtypes = []
    na_subtypes = []
    for filename in os.listdir(path):
        if fnmatch.fnmatch(filename, '*_HA_*.fasta') and os.path.isfile(os.path.join(path, filename)):
            h_name = filename.split('_')[-1].split('.')[0].strip()
            ha_subtypes.append(h_name if h_name else 'Hx')
        if fnmatch.fnmatch(filename, '*_NA_*.fasta') and os.path.isfile(os.path.join(path, filename)):
            n_name = filename.split('_')[-1].split('.')[0].strip()
            na_subtypes.append(n_name if n_name else 'Nx')

    subtype = ''.join(map(str, ha_subtypes)) + ''.join(map(str, na_subtypes))
    if not subtype:
        subtype = 'HxNx'

    # EXTRACT COVERAGE INFORMATION FOR EACH GENE SEGMENT
    gene_cov_data_file = []
    for filename in os.listdir(table_path):
        if filename.endswith('coverage.txt'):
            data = pd.read_table(os.path.join(table_path, filename))
            genelength = len(data)
            name = filename.split('-')[0]
            genename = name.split('_')[1]
            meanv = round(data['Consensus_Count'].mean(), 1)
            datalist.append({'sample': folder + '_' + subtype, 'Gene': genename, 'seq_len': genelength, 'coverage': meanv})

            gene_name = filename.split('_')[1].split('-')[0]
            color = color_map.get(gene_name, 'gray')
            gene_cov_data_file.append((data[['Reference_Name', 'Position', 'Consensus_Count']], gene_name, color))

    file_data[folder] = gene_cov_data_file

    # GET SUMMARY OF READS FROM READ_COUNTS.txt
    reads_count_path = os.path.join(table_path, 'READ_COUNTS.txt')
    reads_data = pd.read_table(reads_count_path)
    reads_data = reads_data[['Record', 'Reads']]
    reads_data['Record'] = reads_data['Record'].str.rsplit('_').str[-1].str.split('-').str[-1]
    reads_data.rename(columns={'Reads': folder + '_' + subtype}, inplace=True)
    reads_data = reads_data.set_index('Record')
    reads_data = reads_data.loc[['initial', 'match']]
    readsdf.append(reads_data.T)


dataframe = pd.DataFrame(datalist)

# Create a dataframe of total reads with assembled reads and percentage assembled
readsdata = pd.concat(readsdf)
readsdata['prop_assembled'] = round(readsdata['match'] / readsdata['initial'] * 100, 1)
readsdata['assembled'] = readsdata['prop_assembled']  # compatibility with legacy plotting path
readsdata.rename(columns={'initial': 'total_reads', 'match': 'total_assembled'}, inplace=True)
readsdata['notassembled'] = readsdata['total_reads'] - readsdata['total_assembled']
readsdata['prop_notassembled'] = round(100 - readsdata['prop_assembled'], 1)
readsdata.reset_index(inplace=True)
readsdata.rename(columns={'index': 'sample'}, inplace=True)


# FUNCTION TO PLOT HEATMAP TO DEPTH
# PROCESS DATAFRAME EXTRACT SEQUENCE LENGTH AND COVERAGE
def arrange_duplicates(data):
    """Append suffix in sample IDs for coinfections per gene."""
    data = data.copy()
    data['suffix'] = data.groupby(['sample', 'Gene']).cumcount().map(lambda x: '' if x == 0 else x + 1)
    data['suffix'] = data['suffix'].apply(lambda x: '-' + str(x) if pd.notna(x) and x != '' else x)
    data['sample'] = data['sample'] + data['suffix']
    data.drop(columns=['suffix'], inplace=True)
    return data


# MERGE CHANGE: split absolute length and proportional length cleanly.
gene_absolute_length = dataframe[dataframe.columns.difference(['coverage'])]
gene_absolute_length = arrange_duplicates(gene_absolute_length)
gene_absolute_length = gene_absolute_length.pivot(index='sample', columns='Gene', values='seq_len').reset_index()
columns = [col for col in geneids if col in gene_absolute_length.columns]
gene_absolute_length = gene_absolute_length[['sample'] + columns]
gene_absolute_length['sample_id'] = gene_absolute_length['sample'].apply(lambda name: name.split('_H')[0])

# Export raw length data
export_data(gene_absolute_length, 'gene_absolute_length_coverage_data.tsv', index=True)

gene_proportion_length = gene_absolute_length.copy()
gene_proportion_length['sample_id'] = gene_proportion_length['sample_id'].apply(lambda value: value.split('_')[0])
gene_proportion_length = gene_proportion_length.set_index('sample_id').fillna(0)

reference_lengths_gene = {}
for sample_name in gene_proportion_length['sample']:
    foldername = sample_name.split('_')[0]
    reference_path = Path(os.path.join(CURRDIR, foldername)) / 'intermediate/0-ITERATIVE-REFERENCES'
    reference_files = reference_path.glob('R0-*.ref')
    reference_lengths = {}
    for gene_reference_file in reference_files:
        gene_name = gene_reference_file.name.split('_')[1].split('.')[0]
        with open(gene_reference_file, 'r') as handle:
            sequence_length = len(handle.readlines()[1].strip('\n'))
        reference_lengths[gene_name] = sequence_length
    reference_lengths_gene[foldername] = reference_lengths

reference_lengths = pd.DataFrame(reference_lengths_gene).T.fillna(0)
reference_lengths = reference_lengths.replace(0, np.nan)
for col in columns:
    if col in reference_lengths.columns:
        gene_proportion_length[col] = round(gene_proportion_length[col] / reference_lengths[col] * 100, 1)
gene_proportion_length = gene_proportion_length.fillna(0)
gene_proportion_length.reset_index(inplace=True)

# CREATE SUBTYPING FILE
subtype_rows = []
for sample in gene_proportion_length['sample']:
    sample_parts = sample.split('_')
    sampleid = sample_parts[0]
    stype = sample_parts[1] if len(sample_parts) > 1 else 'HxNx'
    subtype_rows.append({'SampleID': sampleid, 'Subtype': stype})
subtypes_data = pd.DataFrame(subtype_rows)

# Plot-ready proportional length matrix
gene_proportion_length = gene_proportion_length[['sample'] + columns].set_index('sample')

# GENERATE DATA FOR COVERAGE BY DEPTH
gene_coverage = dataframe[dataframe.columns.difference(['seq_len'])]
gene_coverage = arrange_duplicates(gene_coverage)
df_coverage = gene_coverage.pivot(index='sample', columns='Gene', values='coverage').reset_index()
df_coverage = df_coverage.set_index('sample')
df_coverage = df_coverage[columns]
df_coverage_log = round(np.log10(df_coverage[df_coverage.columns] + 1), 1)


# For every segment get counts of gene numbers
sample_list = []
for folder in os.listdir(CURRDIR):
    filepath = os.path.join(CURRDIR, folder, 'tables', 'READ_COUNTS.txt')
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        try:
            df = pd.read_table(filepath, usecols=['Record', 'Reads'])
            df['Record'] = df['Record'].apply(lambda value: value.split('-')[1] if '-' in value else value)
            df = df[~df['Record'].isin(['initial', 'failQC', 'passQC', 'chimeric', 'nomatch', 'match', 'altmatch'])]
            df['Record'] = df['Record'].apply(lambda value: value.split('_')[1] if '_' in value else value)
            df.dropna(subset='Reads', inplace=True)
            df = df.sort_values(by='Reads').drop_duplicates(subset='Record', keep='last')
            df['sample'] = folder
            df = df.pivot(index='sample', values='Reads', columns='Record').reset_index()
            sample_list.append(df)
        except Exception as err:
            print(f'failed to process {filepath}: {err}')
sample_data = pd.concat(sample_list).fillna(0) if sample_list else pd.DataFrame()


# MERGE CHANGE: function to process generated coverage tables for all-sample faceted plot.
def process_data(fullpath):
    """Process coverage data from *-read-coverage-table.tsv files."""
    files = glob.glob(os.path.join(fullpath, '*-read-coverage-table.tsv'))
    depth_data = []

    for filepath in files:
        depth = pd.read_table(filepath, usecols=['Reference_Name', 'Position', 'Consensus_Count', 'sample_name'])
        depth['Reference_Name'].fillna('-NA', inplace=True)
        depth.dropna(subset='sample_name', inplace=True)
        depth['depth_log'] = round(np.log10(depth['Consensus_Count'] + 1), 2)
        depth_data.append(depth)

    if not depth_data:
        return pd.DataFrame(columns=['Reference_Name', 'Position', 'Consensus_Count', 'sample_name', 'depth_log'])

    return pd.concat(depth_data)


gene_coverage_data_all = process_data(ref_folder_path)

# EXPORT DATAFRAMES
print('Output data to files')
export_data(df_coverage, 'depth_coverage_data.tsv', index=True)
export_data(df_coverage_log, 'depth_coverage_data_log_10.tsv', index=True)
# MERGE CHANGE: keep legacy filename and correctly export proportional (%), not absolute lengths.
export_data(gene_proportion_length, 'gene-percentage-length-coverage-data.tsv', index=True)
export_data(gene_proportion_length, 'gene-proportion-length-coverage-data.tsv', index=True)
export_data(readsdata, 'assembled-reads-data.tsv', index=False)
export_data(sample_data, 'gene_count.tsv', index=False)
export_data(subtypes_data, 'subtypes.tsv', index=False)
# Compatibility file name from legacy pipeline
export_data(subtypes_data, 'influenza_subtypes.tsv', index=True)


print('Plotting assembled vs non-assembled reads')
plot_read_distribution(readsdata, 'prop_assembled', 'prop_notassembled', 'Read Proportion(%)')
plot_read_distribution(readsdata, 'total_assembled', 'notassembled', 'Read Frequency')

print('Plotting per-sample depth curves')
plot_line_coverage_depth(file_data, irma_dir, reference_folder)

print('Plotting figures for coverage by depth logwise')
plot_figures_log(df_coverage_log)
print('Plotting figures for coverage by gene length')
plot_figures_length(gene_proportion_length)
print('Generate read data')
plot_assembled_reads(readsdata)
print('Plotting all coverage lines')
plot_all_line_plots(gene_coverage_data_all)

print('Voila!')
