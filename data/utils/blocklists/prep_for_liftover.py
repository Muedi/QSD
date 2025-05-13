from sys import *

bed_file = argv[1]
out_file = argv[2]

prefix = None
if 'GRCh38' in bed_file:
    prefix = 'hs'
if 'GRCm38' in bed_file:
    prefix = 'mm'

new_bed = ''
with open(bed_file, 'r') as f:
    hsr_ID = 1; lm_ID = 1;
    for line in f:
        chrom, start, end, name = tuple(line.strip().split('\t'))
        start = int(start); end = int(end);

        if name.startswith('High'):
            name = prefix + 'HSR_' + str(hsr_ID)
            hsr_ID += 1
        if name.startswith('Low'):
            name = prefix + 'LM_' + str(lm_ID)
            lm_ID += 1

        new_bed += '%s\n'%('\t'.join(map(str, [chrom, start, end, name])))

open(out_file,'w').write(new_bed)

