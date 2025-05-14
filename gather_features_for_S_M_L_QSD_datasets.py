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

feat_names_RAW = utils.get_FastQC_feature_names()
feat_names_RAW.remove('Basic_Statistics')
feat_names_MAP = ['0times', '1time' , 'multi', 'overall']
feat_names_LOC = utils.get_LOC_feature_names()
feat_names_TSS = utils.get_TSS_feature_names()

s_cols  = [ 'RAW_'+fn for fn in feat_names_RAW ]
s_cols += [ 'MAP_'+fn for fn in feat_names_MAP ]
s_cols += feat_names_LOC
s_cols += feat_names_TSS
s_values = []

ml_cols, ml_values = {}, {}
liftOver_ratios = ['0_50', '0_25']
for ratio in liftOver_ratios:
    blocklist_fp = './data/utils/blocklists/liftover/min_ratio_%s/GRCh38.bed'%(ratio)
    blocklist = pd.read_csv(blocklist_fp, sep='\t', names=['chr', 'start', 'end', 'ID'])
    ml_cols[ratio] = list(blocklist['ID'])
    ml_values[ratio] = []

accessions = []
for index, accession in enumerate(samples_meta['Accession']):
    file_paths = utils.get_feature_file_paths(accession, data_dir)

    # get the values for the S-QSD features
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
        table = pd.read_csv(file_paths['LOC'], sep='\t')
        loc_vals = list(table['Frequency'])
    except:
        print("Couldn't get the LOC features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue

    # read out the ChIPpeakAnno table to receive the TSS features
    try:
        value_map = utils.get_TSS_features(file_paths['TSS']) 
        tss_vals = [ value_map.get(fn, np.nan) for fn in feat_names_TSS ]
    except:
        print("Couldn't get the TSS features for %s"%(accession))
        print('Therefore, sample %s is skipped entirely.'%(accession))
        continue

    s_values.append( raw_vals + map_vals + loc_vals + tss_vals )

    # get the values for the two different blocklists (M/L-QSD)
    successful = True
    for ratio in liftOver_ratios:
        try:
            counts_fp = '%sfeatures/05_BLF/ratio_%s/%s.tsv'%(data_dir, ratio, accession)
            counts = pd.read_csv(counts_fp)
            count_map = dict(zip(counts['blID'], counts['count']))
            row = [ count_map.get(col,0.0) for col in ml_cols[ratio] ]
            ml_values[ratio].append(row)
        except:
            print("Couldn't get the blocklist (%s) features for %s"%(ratio, accession))
            print('Therefore, sample %s is skipped entirely.'%(accession))
            successful = False
    if not successful:
        continue
        
    accessions.append(accession)

# create the S-QSD dataframe and save to CSV file
data = pd.DataFrame( data=s_values, columns=s_cols )
data.insert(0, 'accession', accessions)
data.to_csv('./QSD_datasets/S-QSD.csv', index=False)
print('S-QSD:')
print(data)
print()

# create the M/L-QSD dataframes and save to CSV file
for qsd_size, ratio in zip(['M','L'], liftOver_ratios):
    data = pd.DataFrame( data=ml_values[ratio], columns=ml_cols[ratio] )
    data.insert(0, 'accession', accessions)
    data.to_csv('./QSD_datasets/%s-QSD.csv'%(qsd_size), index=False)
    print('%s-QSD:'%(qsd_size))
    print(data)
    print()
