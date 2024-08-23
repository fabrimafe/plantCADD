import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument('gff_file', help='path to GFF file')
parser.add_argument('out_file', type=str,help='path to output file')
parser.add_argument('--pos', choices=['1', '2', '3'], default='2',
                    help='which position in the codon to extract (default: 2)')
args = parser.parse_args()

codon_pos = int(args.pos) - 1  # convert from 1-based to 0-based index
features = []
out_file=str(args.out_file)
print(out_file)
with open(args.gff_file, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if fields[2] == 'CDS':
                attributes = dict(re.findall('(\S+)\s*=\s*("[^"]*"|\S+)', fields[8]))
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                phase = int(attributes.get('phase', 0))
                length = end - start + 1
                inCDS=True
#                print(fields)
                newstart=start
                newend=end
                while inCDS:
                    if strand == '+':
                         codon_base = newstart + codon_pos - phase #) % 3
                         newstart=newstart+3
                         if (newstart+3)>end:
                            inCDS=False
                    else:
                         codon_base = newend - codon_pos + phase #) % 3
                         newend=newend-3
                         if (newend-3)<start:
                              inCDS=False
                    features.append((fields[0], codon_base-1, codon_base))

df = pd.DataFrame(features, columns=['chromosome', 'codon_start', 'codon_end'])
df.to_csv(args.out_file, sep="\t", header=False,index=False)



