#!/bin/bash
# make a list of excluded sample, in this case non-European samples

# command args:
panel=$1 # 1000G panel file
output=$2

# non-European samples from 1000 Genomes:
grep -v EUR $panel | cut -f1 > $output


