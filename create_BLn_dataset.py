import pandas as pd
from sys import *

def get_assay(assay_title):
    if 'ChIP-seq' in assay_title:
        return 'ChIP-seq'
    elif 'RNA-seq' in assay_title:
        return 'RNA-seq'
    elif 'eCLIP' in assay_title:
        return 'eCLIP'
    else:
        return assay_title


# add assay and organism
meta = pd.read_csv('./data/meta/fastq_samples_meta.csv')
meta['assay'] = [ get_assay(at) for at in meta['Assay title'] ]
meta['organism'] = [ donor.split('/')[1].split('-')[0] for donor in meta['Donor'] ]
meta['status'] = list(meta['Status'])

ratio = float(argv[1])

blocklists_meta = pd.read_csv('./data/meta/blocklist_minMatch_ratio.csv')
blocklists_meta.dropna(inplace=True)

selection = blocklists_meta[blocklists_meta['ratio'] >= ratio]
IDs = list(selection['ID'])

print('\nStarting with %d possible features for minMatch ratio = %.2f\n'%(len(IDs), ratio))

counts = None
try:
    parts = [ pd.read_csv('./data/blocklist_read_counts/data_part%d.csv'%(d)) for d in range(10) ]
    counts = pd.concat(parts)
except:
    counts = pd.read_csv('./data/blocklist_read_counts/blocklist_read_counts_TEST.csv')
counts = counts[['accession'] + IDs]

# quickly check for features that are all zeros
# the minimum number of non-zero values for a feature to be included in the dataset
keep_cols, kick_cols = [], []
for col in counts.columns:
    if col == 'accession':
        continue
    if counts[col].sum() == 0:
        kick_cols.append(col)
    else:
        keep_cols.append(col)
message = '\n1 feature has been exlcuded because it is all zeros:'
if len(kick_cols) > 0:
    if len(kick_cols) > 1:
        message = '%d features have been exlcuded because they are all zeros:'%(len(kick_cols))
    print(message)
    print(', '.join(kick_cols), '\n')

# add assay and organism information to the dataset
counts = counts[['accession'] + keep_cols]
for col in ['status','organism','assay']:
    temp_map = dict(zip(meta['Accession'], meta[col]))
    counts.insert(1, col, counts['accession'].map(temp_map))

print('\nThe final dataset has %d features.'%(len(keep_cols)))
print('(and %d samples)\n'%(counts.shape[0]))

# write the blocklist counts to a csv file
file_name = './QSD_datasets/BL-%s'%(len(keep_cols))
if counts.shape[0] < 100:
    file_name = 'TEST_' + file_name
counts.to_csv(file_name + '.csv', index=False)

print('Double check for missing values:')
print(counts.shape, '\t', counts.dropna().shape, '\n')






