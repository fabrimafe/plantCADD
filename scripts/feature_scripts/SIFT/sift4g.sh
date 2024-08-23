#!/bin/bash
module load jdk/11.0.6

# Args
sift4g_jar_path="/home/labs/alevy/omerbar/bin/sift4gannotator/SIFT4G_Annotator.jar"
db_folder=$1
vcf_file=$2

# Run sift4g
mkdir -p annotations
java -jar ${sift4g_jar_path} -c -i ${vcf_file} -d ${db_folder} -r annotations -t 