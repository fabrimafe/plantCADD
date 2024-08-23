#!/bin/bash

hal_file=$1 #/home/labs/alevy/omerbar/bin/cactus/cactus-bin-v2.7.0/bin/maf2hal
bed_to_convert=$2
output=$3


/home/labs/alevy/omerbar/bin/hal/bin/halLiftover ${hal_file} SL31 ${bed_to_convert} SL5 ${output}