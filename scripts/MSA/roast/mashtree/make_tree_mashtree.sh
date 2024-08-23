#!/bin/bash

fasta_loc=$1

mashtree --mindepth 0 --numcpus 12 ${fasta_loc}/*.fa > mashtree.dnd
