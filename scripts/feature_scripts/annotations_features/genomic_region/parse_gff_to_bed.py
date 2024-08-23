from collections import defaultdict
import sys

def expand_gff_to_feature_bed(input_file, output_file):
    features = set()
    records = defaultdict(lambda: defaultdict(int))
    gene_count = defaultdict(int)

    # Read the unique features and store the ranges
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue # Skip comments
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue # Skip invalid rows
            chr, start, end, region = int(cols[0]), int(cols[3]) - 1, int(cols[4]), cols[2]
            features.add(region)
            for i in range(start, end):
                records[i][region] = 1
                if region == 'gene':
                    gene_count[i] += 1

    # Write the BED file with binary columns for each feature and gene count
    sorted_features = sorted(list(features))
    with open(output_file, 'w') as outfile:
        #header = ['chr','start', 'end'] + sorted_features + ['gene_count']
        #outfile.write('\t'.join(header) + '\n')
        for position, feature_dict in sorted(records.items()):
            row = [chr,position, position + 1]
            row += [feature_dict[feature] for feature in sorted_features]
            row.append(gene_count[position])
            outfile.write('\t'.join(map(str, row)) + '\n')

input_gff=sys.argv[1]
output_bed=sys.argv[2]

expand_gff_to_feature_bed(input_gff, output_bed)
