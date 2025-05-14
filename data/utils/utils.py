import pandas as pd
import os
import subprocess

fastqc_values = {'FAIL':0, 'WARN':1, 'PASS':2}
bl_type_map = {'LM':1, 'HSR':2}


def get_file_length(fp):
    """
    Use the linux function to get the number of lines in a file.

    Args:
        fp (str): File path

    Returns:
        int: number of lines
    """
    wc = None
    try:
        call = 'wc -l %s'%(fp)
        res = subprocess.check_output(call, shell=True, text=True)
        wc = int(res.split()[0])
    except:
        print('Could not get the wc -l here:',fp)
    return wc

def get_feature_file_paths(acc, data_dir):
    """
    Creates file paths for all the output files relevant to the feature generation.

    Args:
        acc (str): ENCODE accession for fast sample.
        data_dir (str): path to the directory containing the data.
                        must contain a "feature" folder.

    Returns:
        dict: dictionary with all file paths.
    """
    paths = {}
    paths['fastq'] = data_dir + 'fastq/' + acc + '.fastq.gz'

    # folder for the FastQC output
    paths['RAW_dir'] = data_dir + 'features/01_RAW/' + acc + '/'
    paths['RAW'] = data_dir + 'features/01_RAW/%s/%s_fastqc/summary.txt'%(acc, acc)
    
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

    # file paths for the blocklist features 
    for ratio in ['0_25','0_50']:
        path_name = '05_BLF_%s'%(ratio)
        paths[path_name] = '%sfeatures/05_BLF/ratio_%s/%s.tsv'%(data_dir, ratio, acc)
    
    return paths

# returns the RAW feature names
def get_FastQC_feature_names():
    names = ['Basic_Statistics', 'Per_base_sequence_quality', 'Per_tile_sequence_quality',
        'Per_sequence_quality_scores', 'Per_base_sequence_content', 'Per_sequence_GC_content',
        'Per_base_N_content', 'Sequence_Length_Distribution', 'Sequence_Duplication_Levels', 
        'Overrepresented_sequences', 'Adapter_Content', 'Kmer_Content']
    return names

# returns the LOC feature names
def get_LOC_feature_names():
    names = ['Promoter', '5_UTR', '3_UTR', '1st_Exon', 'Other_Exon', 
             '1st_Intron', 'Other_Intron', 'Downstream', 'Distal_Intergenic']
    names = ['LOC_' + n for n in names]
    return names

# returns the TSS feature names
def get_TSS_feature_names():
    names  = [ 'TSS_m%d'%d for d in [4500, 3500, 2500, 1500, 500] ]
    names += [ 'TSS_p%d'%d for d in [500, 1500, 2500, 3500, 4500] ]
    return names

# format the TSS feature names correctly
def get_TSS_features(fp):
    table = pd.read_csv(fp, sep='\t')
    table['name'] = ['TSS_'+str(tsd) if '-' in str(tsd) else 'TSS_p'+str(tsd) for tsd in table['tss_dist']]
    table['name'] = [name.replace('-','m') for name in table['name']]
    return dict(zip(table['name'], table['perc']))

def get_FastQC_features(fp):
    """
    Reads FastQC report summary and returns the feature values in a dict.

    Args:
        fp (str): File path to FastQC summary.txt.

    Returns:
        dict: a map of RAW feature names and feature values.
    """
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
    """
    Reads a Bowtie2 report and returns the feature values in a dict.

    Args:
        fp (str): File path to FastQC summary.txt.

    Returns:
        dict: a map of MAP feature names and feature values.
        str: unparsed file content or "not exist" message.
    """
    if not os.path.exists(fp):
        return None, 'File does not exist.'
    lines = open(fp,'r').read().split('\n')
    if len(lines) != 7:
        return None, lines
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
    return stats, lines

def read_blocklist(bl_file):
    """
    Reads in a blocklist file and uses a certain format that allows for 
    using the regions efficiently in read counting procedure

    Args:
        bl_file (str): File path to the blocklist regions file.

    Returns:
        dict: dictionary containing lists of regions for each chromosome used as key
    """
    df_blocklist = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','ID'])
    blocklist = {}
    for index, row in df_blocklist.iterrows():
        if not row['chr'] in blocklist:
            blocklist[row['chr']] = []
        bl_type = row['ID'].split('_')[0][2:]
        blocklist[row['chr']].append( (index+1, row['chr'], row['start'], row['end'],
            bl_type_map[bl_type], row['ID']) )
    return blocklist

def count_reads_in_regions(summits, regions, chrom_size_map):
    """
    For each region, counting the summits within the region. 
    An overlap is only given if the summit is in the region, hence, 
    if more than half of the read overlaps with the blocklist region
    Args:
        summits (dict): dictionary with chromosome as key. The values
                        are lists of summits for the chromosome. The summits 
                        describe the center of the reads in this application.
        regions (dict): dictionary with chromosome as key. The values 
                        are lists of blocklist regions. 
        chrom_size_map (dict): dictionary with chromosome as key. The values
                                describe the chromosome length.
    Returns:
        pd.DataFrame: 
    """
    bincov = {'binID':[], 'chr':[], 'start':[], 'end':[], 'count':[],
        'blID':[], 'blType':[]}

    for chrom in chrom_size_map:
        if not chrom in summits or not chrom in regions:
            continue

        # sort summits to allow for an early exit when counting the summits in the regions
        chr_summits, rp = sorted(summits[chrom]), 0

        for binID, reg_chrom, reg_start, reg_end, reg_type, reg_ID in regions[chrom]:
            if chrom != reg_chrom:
                print('!!! Something is WRONG here: ', reg_chrom, chrom)
                return None
            count = 0
            if rp < len(chr_summits):
                while chr_summits[rp] < reg_end:
                    if chr_summits[rp] > reg_start:
                        count += 1
                    rp += 1
                    if rp >= len(chr_summits):
                        break

            if count != 0:
                bincov['binID'].append( binID )
                bincov['chr'].append( reg_chrom )
                bincov['start'].append( reg_start )
                bincov['end'].append( reg_end )
                bincov['count'].append( count )
                bincov['blType'].append( reg_type )
                bincov['blID'].append( reg_ID )
                    
    return pd.DataFrame(bincov)

