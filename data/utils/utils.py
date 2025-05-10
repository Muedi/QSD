import pandas as pd
import numpy as np
import os
import subprocess

# define the prioritization of different assays and organisms
assay_priority = ['ChIP-seq','RNA-seq','polyA plus RNA-seq','DNase-seq','eCLIP'] #,'Flow-FISH CRISPR screen','single-cell RNA sequencing assay']
orga_priority = ['human','mouse'] #,'fly','worm','manatee','NaN']

# create maps for ranking assay, organism and status
status_pmap = {'revoked':1,'released':2} #,'archived':3}
orga_pmap = dict( (orga,orga_priority.index(orga)+1) for orga in orga_priority )
assay_pmap = dict( (assay,assay_priority.index(assay)+1) for assay in assay_priority )

nhr_data_dir = '/lustre/project/nhr-qsd/data/'

bin_sizes = ['500bp', '1kb', '2kb', '5kb']
fastqc_values = {'FAIL':0, 'WARN':1, 'PASS':2}
bl_type_map = {'LM':1, 'HSR':2}

def priority_sort(files):
    files['prio1'] = [ status_pmap.get(status,99) for status in files['Status'] ]
    files['prio2'] = [ orga_pmap.get(orga,99) for orga in files['organism'] ]
    files['prio3'] = [ assay_pmap.get(assay,99) for assay in files['Assay term name'] ]

    files = files.loc[ [assay in assay_pmap for assay in files['Assay term name']] ]
    files = files.loc[ [orga in orga_pmap for orga in files['organism']] ]
    files = files.loc[ [status in status_pmap for status in files['Status']] ]

    sorted_files = files.copy()
    sorted_files.sort_values(by=['prio1','prio2','prio3'],ignore_index=True,inplace=True)

    return sorted_files

def strict_selection(files):
    files = files.loc[files['Paired end identifier'] != 2]
    files = files.dropna(subset=['Accession'])
    files = files.loc[files['No file available'] == False]
    return files

def get_files_meta_and_maps():
    files = pd.read_csv('./meta/files.csv', low_memory=False)
    
    # get the organism from the Donor ID
    files['organism'] = ['NaN' if type(donor) != str else donor.split('/')[1].split('-')[0] for donor in files['Donor']]
    
    # preprocess the paired-end identifier
    new_peIDs = []
    for peID in files['Paired end identifier']:
        if type(peID) == str:
            if ',' in peID:
                new_peIDs.append( 0 )
            else:
                new_peIDs.append( int(float(peID)) )
        elif type(peID) == float:
            if np.isnan(peID):
                new_peIDs.append( -1 )
    files['Paired end identifier'] = new_peIDs
    
    maps = {}
    maps['assay'] = dict(zip(files['Accession'], files['Assay term name']))
    maps['url'] = dict(zip(files['Accession'], files['Download URL']))
    maps['peID'] = dict(zip(files['Accession'], files['Paired end identifier']))
    maps['organism'] = dict(zip(files['Accession'], files['organism']))
    maps['status'] = dict(zip(files['Accession'], files['Status']))
    maps['size'] = dict(zip(files['Accession'], files['File size']))
    print('Read in the files meta with dimensions:', files.shape)
    print('\tProviding maps for the fields:', ', '.join(map(str,maps.keys())), '\n')
    return files, maps

def get_file_length(fp):
    wc = None
    try:
        call = 'wc -l %s'%(fp)
        res = subprocess.check_output(call, shell=True, text=True)
        wc = int(res.split()[0])
    except:
        print('Could not get the wc -l here:',fp)
    return wc

def get_feature_file_paths(acc, given_dir=None):
    data_dir = nhr_data_dir
    if given_dir != None:
        data_dir = given_dir

    paths = {}
    paths['fastq'] = data_dir + 'fastq/' + acc + '.fastq.gz'

    # folder for the FastQC output
    paths['RAW_dir'] = data_dir + 'features/01_RAW/' + acc + '/'
    
    # files for the mapping statistics, the BAM file and BED files (full and 1M reads randomly sampled)
    paths['MAP_stats'] = data_dir + 'features/02_MAP/stats/' + acc + '.txt'

    paths['MAP_bam'] = data_dir + 'features/02_MAP/bams/' + acc + '.bam'

    paths['MAP_bed_unsorted'] = data_dir + 'features/02_MAP/beds_full/' + acc + '_unsorted.bed'
    paths['MAP_bed'] = data_dir + 'features/02_MAP/beds_full/' + acc + '.bed'

    paths['MAP_bed1M_unsorted'] = data_dir + 'features/02_MAP/beds_oneM/' + acc + '_unsorted.bed'
    paths['MAP_bed1M'] = data_dir + 'features/02_MAP/beds_oneM/' + acc + '.bed'

    # file paths for the LOC and TSS features
    paths['LOC'] = data_dir + 'features/03_LOC/' + acc + '.tsv'
    paths['TSS'] = data_dir + 'features/04_TSS/' + acc + '.tsv'

    paths['BIN_BL'] = data_dir + 'features/05_BIN/v2_BL/' + acc + '.tsv'
    for bs in bin_sizes:
        paths['BIN_%s'%(bs)] = data_dir + 'features/05_BIN/%s/%s.tsv'%(bs,acc)

    return paths

def get_FastQC_feature_names():
    names = ['Basic_Statistics', 'Per_base_sequence_quality', 'Per_tile_sequence_quality',
        'Per_sequence_quality_scores', 'Per_base_sequence_content', 'Per_sequence_GC_content',
        'Per_base_N_content', 'Sequence_Length_Distribution', 'Sequence_Duplication_Levels', 
        'Overrepresented_sequences', 'Adapter_Content', 'Kmer_Content']
    return names

def get_LOC_feature_names():
    names = ['Promoter', '5_UTR', '3_UTR', '1st_Exon', 'Other_Exon', 
             '1st_Intron', 'Other_Intron', 'Downstream', 'Distal_Intergenic']
    names = ['LOC_' + n for n in names]
    return names

def get_TSS_feature_names():
    names  = [ 'TSS_m%d'%d for d in [4500, 3500, 2500, 1500, 500] ]
    names += [ 'TSS_p%d'%d for d in [500, 1500, 2500, 3500, 4500] ]
    return names

def get_FastQC_features(fp):
    try:
        value_map = {}
        with open(fp,'r') as f:
            for line in f:
                line = line.strip().split('\t')
                value_map[line[1].replace(' ','_')] = fastqc_values[line[0]]
        return value_map
    except:
        return None

def read_Bowtie_stats(fp):
    if not os.path.exists(fp):
        return None
    lines = open(fp,'r').read().split('\n')
    if len(lines) != 7:
        return None
    stats = {'total': int(lines[0].split()[0]) }
    stats['unpaired'] = int(lines[1].split()[0])
    stats['0times'] = int(lines[2].split()[0])
    stats['1time'] = int(lines[3].split()[0])
    stats['multi'] = int(lines[4].split()[0])
    stats['overall'] = stats['1time'] + stats['multi']
    percentages = {}
    for key in stats:
        if key != 'total':
            percentages['perc_'+key] = stats[key] / stats['total'] * 100.0
    stats.update(percentages)
    return stats

def read_blacklist(bl_file):
    df_blacklist = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','ID'])
    blacklist = {}
    for index, row in df_blacklist.iterrows():
        if not row['chr'] in blacklist:
            blacklist[row['chr']] = []
        bl_type = row['ID'].split('_')[0][2:]
        blacklist[row['chr']].append( (index+1, row['chr'], row['start'], row['end'],
            bl_type_map[bl_type], row['ID']) )
    return blacklist



