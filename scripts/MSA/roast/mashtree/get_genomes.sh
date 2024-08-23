#!/bin/bash

name_file=$1

while IFS= read -r line; do
    cp /home/labs/alevy/omerbar/expanded_alignment_sl5/genome_fastas/${line}.fa . &

done < ${name_file}