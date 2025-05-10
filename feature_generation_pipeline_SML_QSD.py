import pandas as pd
from sys import *
import data.utils.utils as utils
import os

samples_meta = pd.read_csv('./data/meta/fastq_samples_meta.csv')
samples_meta['Organism'] = [donor.split('/')[1].split('-')[0] for donor in samples_meta['Donor']]
print(set(samples_meta['Organism']))

accession = argv[1]
organsim  = dict(zip(samples_meta['Accession'], samples_meta['Organism']))[accession]

print('Starting to generate features for sample %s - the organism is %s'%(accession, organsim))




