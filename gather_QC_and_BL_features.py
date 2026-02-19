"""
This script reads in tables and statistics generated for the features
for the S/M/L-QSD datasets. 
It uses the same utils as used for the feature generation and expects to find 
the output files from the different tools in the "features" following a certain 
scheme as defined in the utils funtion "get_feature_file_paths()".
The script uses the meta data for the QSD datasets to iterate over all
samples to gather all feature values.
"""

from genericpath import exists
import pandas as pd
import numpy as np

# utils implemented for the data generation of the QSD datasets
import data.utils.utils as utils

data_dir = './data/'
samples_meta = pd.read_csv('%smeta/fastq_samples_meta.csv'%(data_dir))
samples_meta['Organism'] = [donor.split('/')[1].split('-')[0] for donor in samples_meta['Donor']]
data_dir = './data/'

# collect the feature names for the different sets
feat_names_RAW = utils.get_FastQC_feature_names()
feat_names_RAW.remove('Basic_Statistics')
feat_names_MAP = ['0times', '1time' , 'multi', 'overall']
feat_names_LOC = utils.get_LOC_feature_names()
feat_names_TSS = utils.get_TSS_feature_names()

# prepare lists that are considered as columns for the final datasets
qc_cols  = [ 'RAW_'+fn for fn in feat_names_RAW ]
qc_cols += [ 'MAP_'+fn for fn in feat_names_MAP ]
qc_cols += feat_names_LOC
qc_cols += feat_names_TSS
qc_values = []

# initially, this script was supposed to run for a fixed list of different ratios 
# for the liftOver blocklists. However, we now derive the counts for all blocklists 
# we receive with a minMatch ratio of 0.10.
# From these counts, datasets for different ratios can be created (see create_BLn_dataset.py)
blocklist_fp = './data/utils/blocklists/liftover/min_ratio_0_10/GRCh38.bed'
blocklist = pd.read_csv(blocklist_fp, sep='\t', names=['chr', 'start', 'end', 'ID'])
bl_cols = list(blocklist['ID'])
bl_values = []

# iterate over all accessions that could potentially be in the dataset
possible_accessions = samples_meta['Accession']

# this is just for testing to show how the script has been used to gather the features
# and assemble the different datasets
possible_accessions = ['ENCFF001NAO', 'ENCFF001NFW']

accessions = []
for index, accession in enumerate(possible_accessions):
    file_paths = utils.get_feature_file_paths(accession, data_dir)

    # get the values for the QC features
    feature_vals = {}
    raw_vals, map_vals, loc_vals, tss_vals = None, None, None, None

    # read out the FastQC report to receive the RAW features
    try:
        value_map = utils.get_FastQC_features(file_paths['RAW']) 
        raw_vals = [ value_map.get(fn, np.nan) for fn in feat_names_RAW ]
    except:
        print("Couldn't get the RAW features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue

    # read out the Bowtie2 statsitics to to receive the MAP features
    try:
        map_stats, file_content = utils.read_Bowtie_stats(file_paths['MAP_stats'])
        map_vals = [ map_stats['perc_'+name] for name in feat_names_MAP ]
    except:
        print("Couldn't get the MAP features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue
    
    # read out the ChIPseeker table to receive the LOC features
    try:
        value_map = utils.get_LOC_features(file_paths['LOC'])
        loc_vals = [ value_map.get(fn, 0.0) for fn in feat_names_LOC ]
    except:
        print("Couldn't get the LOC features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue

    # read out the ChIPpeakAnno table to receive the TSS features
    try:
        value_map = utils.get_TSS_features(file_paths['TSS']) 
        tss_vals = [ value_map.get(fn, 0.0) for fn in feat_names_TSS ]
    except:
        print("Couldn't get the TSS features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue

    # get the values for the two different blocklists (M/L-QSD)
    bl_row = None
    try:
        counts_fp = '%sfeatures/05_BLF/ratio_0_10/%s.tsv'%(data_dir, accession)
        counts = pd.read_csv(counts_fp)
        count_map = dict(zip(counts['blID'], counts['count']))
        bl_row = [ count_map.get(col, 0.0) for col in bl_cols ]
    except:
        print("Couldn't get the blocklist features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue
    
    qc_values.append( raw_vals + map_vals + loc_vals + tss_vals )
    bl_values.append(bl_row)
    accessions.append(accession)

# create the QC-QSD dataframe and save to CSV file
data = pd.DataFrame( data=qc_values, columns=qc_cols )
data.insert(0, 'accession', accessions)
data.to_csv('./QSD_datasets/QSD-QC.csv', index=False)
print('QSD-QC:')
print(data)
print()

# create the M/L-QSD dataframes and save to CSV file
data = pd.DataFrame( data=bl_values, columns=bl_cols )
data.insert(0, 'accession', accessions)
data.to_csv('./QSD_datasets/QSD-BL-all.csv', index=False)
print('QSD-BL-all:')
print(data)
print()
