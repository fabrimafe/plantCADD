#!/bin/bash

input=$1
output=$2

sed "s/:[0-9\.-]*//g" ${input} | sed "s/,/ /g"> ${output}