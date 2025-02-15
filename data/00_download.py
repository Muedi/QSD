import pandas as pd
from sys import *
import os

seqQscorer_fastq = '/lustre/project/m2_jgu-cbdm/salbrec/seqQscorer_proj/data/raw_fastq/zipped/'

files = pd.read_csv('./meta/files_essential.tsv', sep='\t')
print(files)

not_skip = int(argv[1])

files['organism'] = ['NaN' if type(donor) != str else donor.split('/')[1].split('-')[0] for donor in files['Donor']]

assay_priority = ['ChIP-seq','polyA plus RNA-seq','RNA-seq','DNase-seq','eCLIP','Flow-FISH CRISPR screen','single-cell RNA sequencing assay']
assay_pmap = dict( (assay,assay_priority.index(assay)+1) for assay in assay_priority )
orga_priority = ['human','mouse','fly','worm','manatee','NaN']
orga_pmap = dict( (orga,orga_priority.index(orga)+1) for orga in orga_priority )

assay_map = dict(zip(files['Accession'], files['Assay term name']))
url_map = dict(zip(files['Accession'], files['Download URL']))
peID_map = dict(zip(files['Accession'], files['Paired end identifier']))
orga_map = dict(zip(files['Accession'], files['organism']))
status_map = dict(zip(files['Accession'], files['Status']))

files['prio1'] = [ orga_pmap.get(orga,99) for orga in files['organism'] ]
files['prio2'] = [ assay_pmap.get(assay,99) for assay in files['Assay term name'] ]

files_sorted = files.copy()
files_sorted.sort_values(by=['prio1','prio2'],ignore_index=True,inplace=True)

for index, acc in enumerate(files_sorted['Accession']):
    peID = peID_map[acc]
    print(index, acc, assay_map[acc], orga_map[acc], peID)
    
    if index % 4 != not_skip:
        print('\tSkip because of "modulo(index) != %d"'%(not_skip))
        continue
    
    if peID == 2.0:
        print('\tSkip because the paired-end identifier is', peID)
        continue

    if status_map[acc] == 'archived':
        print('\tSkip because it is archived', peID)
        continue
    
    fp = './fastq/%s.fastq.gz'%(acc)
    if os.path.exists(fp):
        print('\tSkip because it does exist already:', fp)
        continue

    seqQscorer_fp = '%s%s.fastq.gz'%(seqQscorer_fastq, acc)
    if os.path.exists(seqQscorer_fp):
        print('\tCopy from seqQscorer!')
        print(seqQscorer_fp)
        cp = 'cp %s ./fastq/'%(seqQscorer_fp)
        os.system(cp)
        continue

    print('Downloading...')
    wget =  'wget https://www.encodeproject.org' + url_map[acc]
    wget += ' -P ./fastq/'
    print(wget)
    os.system(wget)

    print()




