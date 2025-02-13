import os
from sys import argv

def read_stats(fp):
	lines = open(fp,'r').read().split('\n')
	stats = {'total': int(lines[0].split()[0]) }
	stats['unpaired'] = int(lines[1].split()[0])
	stats['0times'] = int(lines[2].split()[0])
	stats['1time'] = int(lines[3].split()[0])
	stats['multi'] = int(lines[4].split()[0])
	stats['overall'] = stats['1time'] + stats['multi']
	percentages = {}
	for key in stats:
		if key == 'total':
			continue
		percentages['perc_'+key] = stats[key] / stats['total'] * 100.0
	stats.update(percentages)
	return stats

genome = 'GRCh38' # could also be mm10 
# for seqQscorer there was this incosisetnecy using GRCh38 but then mm10 and not GRCm??
# I think we had some problems with this ...

bowtie2_idx_HS = '../../seqQscorer_proj/get_data/Bowtie/genomes/Homo_sapiens/GRCh38/index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
bowtie2_idx = '../../../secondPD/Alignment_stuff/bowtie2_indices/Mouse/Mus_musculus.GRCm38'
# on the local machine, I just used the mm10 as test
if os.path.exists('/lustre/project/m2_jgu-cbdm/salbrec/'):
	bowtie2_idx = bowtie2_idx_HS

accession = argv[1]

# for this test, the fastq files are stored here:
# /lustre/project/m2_jgu-cbdm/salbrec/anno_detect/data/raw_FASTQ/
fastq_file = '../data/raw_FASTQ/%s.fastq.gz'%(accession)
map_uns_bam_file = './02_MAP/%s.unsorted.bam'%(accession)

# unsorted bed files
map_uns_bed_file = './02_MAP/%s.unsorted.bed'%(accession)
map_uns_bed1M_file = './02_MAP/%s_1M.unsorted.bed'%(accession)

# this will be the input for the binning
map_bed_file = './02_MAP/%s.bed'%(accession)

# this has always been the input for the R scripts deriving the LOC and TSS features
map_bed1M_file = './02_MAP/%s_1M.bed'%(accession)

map_stats_file = './02_MAP/%s_stats.txt'%(accession)

bowtie2  = 'time bowtie2 -p 1 -x %s '%(bowtie2_idx)
bowtie2 += '-U %s 2> %s | '%(fastq_file, map_stats_file)
bowtie2 += 'samtools view -@ 1 -Sb -o %s'%(map_uns_bam_file)
print('\n' + bowtie2 + '\n')
if not os.path.exists(map_uns_bam_file):
	os.system(bowtie2)

print('Mapping is done, these are the stats:')
map_stats = read_stats(map_stats_file)
for k, v in map_stats.items():
	print('\t', k,'\t', v)
print()

bamtobed = 'time bedtools bamtobed -i %s > %s'%(map_uns_bam_file, map_uns_bed_file)
print(bamtobed, '\n')
os.system(bamtobed)

shuf = 'time shuf -n 1000000 %s > %s'%(map_uns_bed_file, map_uns_bed1M_file)
print(shuf, '\n')
os.system(shuf)

# sort the 2 bed files:
bedsort = 'time bedtools sort -i %s > %s'%(map_uns_bed_file, map_bed_file)
print(bedsort, '\n')
os.system(bedsort)

bedsort = 'time bedtools sort -i %s > %s'%(map_uns_bed1M_file, map_bed1M_file)
print(bedsort, '\n')
os.system(bedsort)

binning  = 'time python 02_binit.py %s %s %d %d'%(accession, genome, map_stats['total'], map_stats['overall'])
print(binning, '\n')
os.system(binning)






