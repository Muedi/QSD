import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

cis = pd.read_csv('CistromeDB_ChIPseq_ENCODE.csv')
cistrome_quality_flags = ["PeaksFoldChangeAbove10","FRiP", "FastQC","PeaksUnionDHSRatio", "PBC"]

rel = cis.loc[cis['ENCODE_status'] == 'released'].reset_index(drop = True)
rev = cis.loc[cis['ENCODE_status'] == 'revoked'].reset_index(drop = True)

p_values = []
for quality_flag in cistrome_quality_flags:

    rel_sub = rel[["ENCODE_status", quality_flag]].dropna() 
    rev_sub = rev[["ENCODE_status", quality_flag]].dropna()

    rel_sub = rel_sub.reset_index(drop = True) 
    rev_sub = rev_sub.reset_index(drop = True) 

    stat, p = mannwhitneyu(rel_sub[quality_flag], rev_sub[quality_flag], alternative='two-sided')

    p_values.append(p)

cistrome_flags_p_values = dict(zip(cistrome_quality_flags, p_values))

p_corrected = multipletests(p_values, alpha = 0.05, method = "holm")
print(cistrome_flags_p_values, p_corrected)

# print('\n\nNumber of ChIP-seq samples (human and mouse): %d'%(cis.shape[0]))
# print('Percentage of revoked: %.2f\n'%(rev.shape[0] / cis.shape[0] * 100.0))

print('Median PeaksFoldChangeAbove10 for released:', rel['PeaksFoldChangeAbove10'].median())
print('Median PeaksFoldChangeAbove10 for revoked: ', rev['PeaksFoldChangeAbove10'].median())
print()

print('Median FRiP for released:', rel['FRiP'].median())
print('Median FRiP for revoked: ', rev['FRiP'].median())
print()

print('Median PBC for released:', rel['PBC'].median())
print('Median PBC for revoked: ', rev['PBC'].median())
print()

print('Median FastQC for released:', rel['FastQC'].median())
print('Median FastQC for revoked: ', rev['FastQC'].median())
print()

print('Median PeaksUnionDHSRatio for released:', rel['PeaksUnionDHSRatio'].median())
print('Median PeaksUnionDHSRatio for revoked: ', rev['PeaksUnionDHSRatio'].median())
print()