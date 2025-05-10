import pandas as pd
from sys import *

vol = argv[1]

columns = ['chr','start','end','ID']

origs = ['GRCh38','GRCm38']
los = ['mm_to_hs','hs_to_mm']
same_IDs = {'GRCh38':'hs_to_mm', 'GRCm38':'mm_to_hs'}

lo_data = {}
for lo in los:
    lo_data[lo] = pd.read_csv( './%s/lo_output/%s.bed'%(vol, lo), sep='\t', names=columns )
    print(lo)
    print(lo_data[lo])

orig_data = {}
for orig in origs:
    temp = pd.read_csv( '../v2_original/new_IDs/%s.bed'%(orig), sep='\t', names=columns )

    counter_IDs = set(lo_data[same_IDs[orig]]['ID'])
    temp = temp.loc[ [ ID in counter_IDs for ID in temp['ID'] ] ]

    orig_data[orig] = temp
    print(orig)
    print(orig_data[orig])

for orig, lo in zip(origs, los):
    print(orig, lo)
    new = pd.concat( [orig_data[orig], lo_data[lo]] )
    new.to_csv('./%s/%s.bed'%(vol, orig), sep='\t', index=False, header=False)