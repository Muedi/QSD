"""
Script for generating the features for the S/M/L-QSD datasets.
It takes one command line argument that is the ENCODE accession for a 
NGS samples in fastq format. This script has access to the full meta 
information required for the QSD data and automatically derives all features.
Use "-sim" to simulate the feature generation of a sample. If the simulation 
is switched on, the linux calls are only printed, not executed.

After most processing steps, there are additional steps to check if
output files exist, the content is expected and feature values can be read out.
"""

import pandas as pd
from sys import argv
from os.path import exists
from os import system, mkdir, stat

# utils implemented for the data generation of the QSD datasets
import data.utils.utils as utils

data_dir = './data/'
samples_meta = pd.read_csv('%smeta/fastq_samples_meta.csv'%(data_dir))
samples_meta['Organism'] = [donor.split('/')[1].split('-')[0] for donor in samples_meta['Donor']]
samples_meta = samples_meta.set_index('Accession')

# first receive the fastq sample accession from the command line
# and get the organism (human or mouse) from the meta data
accession = argv[1]
simulate = '-sim' in argv
organism  = samples_meta.at[accession,'Organism']
n_cores = 1

print('Starting to generate features for sample %s - the organism is %s'%(accession, organism))

file_paths = utils.get_feature_file_paths(accession, data_dir)

# these are the genome builts for the Bowtie2 mapper
# such genomes are large and are expected to be built by the user
# according to their requirements
bowtie_genomes = {}
bowtie_genomes['mouse'] = './data/utils/genomes/mouse_mm10/bowtie2_index/mm10'
bowtie_genomes['human'] = './data/utils/genomes/human_GRCh38/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
# these are the assembly names required by the Bioconductor scripts 
assembly = {'human':'GRCh38', 'mouse':'GRCm38'}

# ==== Download fastq sample ====
if not exists(file_paths['fastq']):
	print('Downloading the file into:', file_paths['fastq'])
	download_url = samples_meta.at[accession,'Download URL']
	wget =  'wget https://www.encodeproject.org%s -P %sfastq/'%(download_url, data_dir)
	print('wget:', wget)
	if not simulate:
		system(wget)

if exists(file_paths['fastq']):
	fastq_size = stat(file_paths['fastq']).st_size
	if fastq_size != samples_meta.at[accession,'File size']:
		print('Download failed, file size does not match.')
		exit(1)
else:
	print('Download failed, file does not exist.')
	exit(1)

# ==== Derive RAW features using FastQC ====
feats_RAW = utils.get_FastQC_features(file_paths['RAW'])
if feats_RAW == None:
	print('FastQC report could not be read in.')
	mkdir(file_paths['RAW_dir']) if not exists(file_paths['RAW_dir']) else None
	fastqc = 'fastqc %s --extract -d %s -o %s'%(file_paths['fastq'], file_paths['RAW_dir'], file_paths['RAW_dir'])
	print('Running FastQC to create the RAW features.')
	print('FastQC:', fastqc)
	if not simulate:
		system(fastqc)

feats_RAW = utils.get_FastQC_features(file_paths['RAW'])
if feats_RAW == None:
	print('Something seems to be wrong with the FastQC report')
	exit(1)

# ==== Reads Mapping with Bowtie2 ====
map_stats, file_content = utils.read_Bowtie_stats(file_paths['MAP_stats'])
if map_stats == None or not exists(file_paths['MAP_bam']):
	print('Mapping does not exists.')
	bowtie2  = 'bowtie2 -p %d -x %s '%(n_cores, bowtie_genomes[organism])
	bowtie2 += '-U %s 2> %s | '%(file_paths['fastq'], file_paths['MAP_stats'])
	bowtie2 += 'samtools view -@ %d -Sb -o %s'%(n_cores, file_paths['MAP_bam'])
	print('Running Bowtie2 to create the MAP features and receive the Mapping in BAM file format.')
	print('Bowtie2:', bowtie2)
	if not simulate:
		system(bowtie2)
map_stats, file_content = utils.read_Bowtie_stats(file_paths['MAP_stats'])
if map_stats == None or not exists(file_paths['MAP_bam']):
	print('Something seems to be wrong with the MAP features or the mapping itself.')
	print('Details from the mapping stats file:')
	print(file_content)
	exit(1)

# ==== Convert BAM file and sample 1M reads ====
if not exists(file_paths['MAP_bed_unsorted']):
	print('Converting BAM file to BED.')
	bamtobed = 'bedtools bamtobed -i %s > %s'%(file_paths['MAP_bam'], file_paths['MAP_bed_unsorted'])
	print('BEDtools:', bamtobed)
	if not simulate:
		system(bamtobed)

# check if file exists, then check if it has the right content 
# by comparing the number of reads in the BED file (number of lines)
# with the expected number of mapped reads from the mapping statistics 
n_mapped_reads = map_stats['1time'] + map_stats['multi']
if exists(file_paths['MAP_bed_unsorted']):
	bed_lines = utils.get_file_length(file_paths['MAP_bed_unsorted'])
	if bed_lines != n_mapped_reads:
		system('rm ' + file_paths['MAP_bed_unsorted'])
		print('BED file is not complete. Expected are %d reads, has %d.'%(n_mapped_reads, bed_lines))
		exit(1)
else:
	print('Failed to convert BAM to BED')
	exit(1)

# create and sort the 
if not exists(file_paths['MAP_bed1M']):
	shuf = 'shuf -n 1000000 %s > %s'%(file_paths['MAP_bed_unsorted'], file_paths['MAP_bed1M_unsorted'])
	print('Sample 1M reads from full BED:', shuf)
	if not simulate:
		system(shuf)
	bedsort = 'bedtools sort -i %s > %s'%(file_paths['MAP_bed1M_unsorted'], file_paths['MAP_bed1M'])
	print('Sort downsampled BED:', bedsort)
	if not simulate:
		system(bedsort)
		# after the subsampled BEd has been sorted, the unsorted BED is not kept
		system('rm ' + file_paths['MAP_bed1M_unsorted'])

if not exists(file_paths['MAP_bed1M']):
	print('Failed to create the subsampled BED file.')
	exit(1)

# check if the subsamples BED file has been created successfully
# this is done by comparing the number of reads (lines in the BED)
# with the expected number of either 1M reads or the number of mapped reads
bed_lines = utils.get_file_length(file_paths['MAP_bed1M'])
if bed_lines != n_mapped_reads and bed_lines != 1000000:
	system('rm ' + file_paths['MAP_bed1M'])
	print('oneM BED has the unexpected content with %d lines'%(bed_lines))
	exit(1)

# ==== Run ChIPseeker to derive the LOC features ====
feats_LOC, feats_TSS = None, None 
try:
	feats_LOC = pd.read_csv(file_paths['LOC'], sep='\t')
except:
	print('Running ChIPseeker to create the LOC features.')
	chipseeker = 'Rscript ./data/utils/get_LOC_features.R %s %s %s'%(file_paths['MAP_bed1M'],
						assembly[organism], file_paths['LOC'])
	print('ChIPseeker:', chipseeker)
	if not simulate: 
		system(chipseeker)

# double check if the LOC features were created successfully
try:
	feats_LOC = pd.read_csv(file_paths['LOC'], sep='\t')
except:
	print('Something went wrong while creating the LOC features.')
	system('rm ' + file_paths['LOC'])
	exit(1)

# ==== Run ChIPpeakAnno to derive the TSS features ====
try:
	feats_TSS = pd.read_csv(file_paths['TSS'], sep='\t')
except:
	print('Running ChIPpeakAnno to create the TSS features.')
	chippeakanno = 'Rscript ./data/utils/get_TSS_features.R %s %s %s'%(file_paths['MAP_bed1M'],
				assembly[organism], file_paths['TSS'])
	print('ChIPpeakAnno:', chippeakanno )
	if not simulate:
		system(chippeakanno)

# double check if the TSS features were created successfully
try:
	feats_TSS = pd.read_csv(file_paths['TSS'], sep='\t')
except:
	print('Something went wrong while creating the TSS features.')
	system('rm ' + file_paths['TSS'])
	exit(1)

# ==== Derive the blocklist features ====

# first read in the chomosome sizes. the sizes are not needed, 
# however, it provides a list of relevant chromosomes.
fp_sizes = '%sutils/chromosome_sizes/%s.tsv'%(data_dir, assembly[organism])
chrom_sizes = pd.read_csv(fp_sizes, sep='\t', names=['chr', 'size'])
chrom_sizes = chrom_sizes.loc[ [not c in ['chrX', 'chrY'] for c in chrom_sizes['chr']] ]
chrom_size_map = dict(zip(chrom_sizes['chr'], chrom_sizes['size']))

# the summits describe the center of the mapped read locations in the full BED
summits = dict((chrom, []) for chrom in chrom_size_map.keys())
with open(file_paths['MAP_bed_unsorted'], 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        chrom = line[0]
        if chrom in summits:
            summits[chrom].append( int((int(line[1])+int(line[2]))/2.0) )

# count the reads overlapping with the blocklisted regions
# using the summits, a overlap of a read is only given if the summit is within the blocklist
# region. hence, if more than half of the read overlaps with the blocklist region
for ratio in ['0_25','0_50']:
	try:
		bl_file = '%sutils/blocklists/liftover/min_ratio_%s/%s.bed'%(data_dir, ratio, assembly[organism])
		blocklist = utils.read_blocklist(bl_file)
		count_BL_reg = utils.count_reads_in_regions(summits, blocklist, chrom_size_map)
		count_BL_reg.to_csv(file_paths['05_BLF_%s'%(ratio)], index=False)
	except:
		print('Failed to create the blocklist features for ratio', ratio)
		exit(1)

print('Successfully processed the sample:', accession)
exit(0)



