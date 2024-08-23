#!/bin/bash
# First I append values to stop loss and gain to the sift score
module load tabix

input=$1
output=$2

cat ${input} | \
awk '{
    if ($8=="NONCODING") {next}
    if ($8=="STOP-LOSS" || $8=="STOP-GAIN") {
        $10=0
        $11=last_median
        $12=last_numseqs
        $13="DELETERIOUS"
    } 
    last_median=$11
    last_numseqs=$12
    if (last_median == "NA") {last_median=-1}  
    if (last_numseqs == "NA") {last_numseqs=1} 

    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' | \
bgzip >  ${output}.gz
