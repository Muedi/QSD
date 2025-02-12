from sys import *
import pandas as pd
import time

def find_overlap(start, end, regions):
    overlap = None
    for reg in regions:
        reg_start, reg_end = reg[2], reg[3]
        if reg_end > start:
            if reg_start < end:
                overlap = reg
                break
        if reg_start > end:
            break
    return overlap

def bincov_for_bins(summits, chrom_size_map, bin_size, reads_total,
    reads_mapped, blacklist=None):

    total = reads_total / 1e6
    mapped = reads_mapped / 1e6

    # start the binning and count summits within the bin ranges
    binID = 1
    bincov = {'binID':[], 'chr':[], 'start':[], 'end':[], 'count':[],
        'cRelTotal':[], 'cRelMapped':[]}
    if blacklist != None:
        bincov['blacklist'] = []

    for chrom, chr_size in chrom_size_map.items():
        chr_summits, rp = sorted(summits[chrom]), 0
        bin_start = 0
        while bin_start < chr_size:
            bin_end = bin_start + bin_size
            bin_end = bin_end if bin_end < chr_size else chr_size

            count = 0
            if rp < len(chr_summits):
                while chr_summits[rp] < bin_end:
                    count += 1
                    rp += 1
                    if rp >= len(chr_summits):
                        break

            if count != 0:
                bincov['binID'].append( binID )
                bincov['chr'].append( chrom )
                bincov['start'].append( bin_start )
                bincov['end'].append( bin_end )
                bincov['count'].append( count )
                bincov['cRelTotal'].append( count/total )
                bincov['cRelMapped'].append( count/mapped )
                if blacklist != None:
                    annotation = 0
                    if chrom in blacklist:
                        overlap = find_overlap(bin_start, bin_end, blacklist[chrom])
                        if overlap != None:
                            annotation = overlap[-1]
                    bincov['blacklist'].append( annotation )

            bin_start += bin_size
            binID += 1

    return pd.DataFrame(bincov)

def bincov_for_BL(summits, regions, reads_total, reads_mapped):

    total = reads_total / 1e6
    mapped = reads_mapped / 1e6

    # start the binning and count summits within the bin ranges
    bincov = {'binID':[], 'chr':[], 'start':[], 'end':[], 'count':[],
        'cRelTotal':[], 'cRelMapped':[], 'blacklist':[]}

    for chrom in chrom_size_map:
        chr_summits, rp = sorted(summits[chrom]), 0

        for binID, reg_chrom, reg_start, reg_end, reg_type in regions[chrom]:
            if chrom != reg_chrom:
                print('!!! Something is WRONG here: ', reg_chrom, chrom)
                exit(-1)
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
                bincov['cRelTotal'].append( count/total )
                bincov['cRelMapped'].append( count/mapped )
                bincov['blacklist'].append( reg_type )
                    
    return pd.DataFrame(bincov)

# example:
# python 02_binit.py ENCFF000VQC GRCh38 3310365 1480616

accession = argv[1]
genome = argv[2]
reads_total = int(argv[3])
reads_mapped = int(argv[4])

# read in the chromosome sizes
fp_sizes = './chromosome_sizes/%s.tsv'%(genome)
chrom_sizes = pd.read_csv(fp_sizes, sep='\t', names=['chr', 'size'])
chrom_sizes = chrom_sizes.loc[ [not c in ['chrX', 'chrY'] for c in chrom_sizes['chr']] ]
chrom_size_map = dict(zip(chrom_sizes['chr'], chrom_sizes['size']))

# first get the summits of the reads
bed_file = './02_MAP/%s.bed'%(accession)
summits = dict((chrom, []) for chrom in chrom_size_map.keys())
with open(bed_file, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        chrom = line[0]
        if chrom in summits:
            summits[chrom].append( int((int(line[1])+int(line[2]))/2.0) )

# load and preprocess blacklists
bl_type_map = {'Low Mappability':1, 'High Signal Region':2}
bl_file = './blacklists/v2_%s.bed'%(genome)
df_blacklist = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','type'])
blacklist = {}
for index, row in df_blacklist.iterrows():
    if not row['chr'] in blacklist:
        blacklist[row['chr']] = []
    blacklist[row['chr']].append( (index+1, row['chr'], row['start'], row['end'],
        bl_type_map[row['type']]) )

for bin_size in ['500bp', '1kb', '2kb', '5kb']:
    # the bin size can be wither in bp or kb
    start_time = time.time()
    bin_size_int = int(bin_size.replace('bp','').replace('kb','000'))
    bincov = bincov_for_bins(summits, chrom_size_map, bin_size_int,
                reads_total, reads_mapped, blacklist=blacklist)
    
    out_file = './05_BIN/%s/%s.csv'%(bin_size, accession)
    bincov.to_csv(out_file, index=False)
    time_needed = time.time() - start_time 
    print('Done with %s (%.2f seconds)'%(bin_size, time_needed))

BL_bincov = bincov_for_BL(summits, blacklist, reads_total, reads_mapped)
    
out_file = './05_BIN/BL/%s.csv'%(accession)
BL_bincov.to_csv(out_file, index=False)

