#!/bin/bash

# Input file containing the Newick tree
input_file=$1

# Output file where the modified Newick tree will be stored

# Regular expression pattern to match numbers, colons, commas, and semicolons
pattern="[0-9:,.;]"

# Remove numbers, colons, commas, and semicolons from the input file and save the result to the output file
sed "s/$pattern//g" "$input_file"
