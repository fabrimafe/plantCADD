#!/bin/bash

xls_file=$1
zcat ${xls_file} | \
# This fucking program gives windows style cartridge returns (\r) that needs to be fucking removed
sed 's/\r$//' | \
# Get relevant columns
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$11"\t"$9"\t"$12"\t"$13"\t"$14"\t"$15"\t"$17}' | \

bgzip > SIFT_annotations.xls 