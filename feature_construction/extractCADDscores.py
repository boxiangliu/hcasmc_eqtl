#!/usr/bin/env python

# author: Emily Tsang
# modified by Boxiang Liu to run for HCASMC eQTL

"""
Modeled after script provided on CADD website written by Martin Kircher

Reads a bed file from std in with the following columns:
chr pos0 pos1 MAF allele1 allele2

Gets the CADD score for each allele (only exists the non reference ones) and records the max of the two

Prints:
chr pos0 pos1 MAF CADD

Run as follows:
cat input.bed | python extractCADDscores.py 
"""

filename = "/mnt/lab_data/montgomery/shared/CADD/whole_genome_SNVs.tsv.gz"

import pysam
import sys

# BED FIELDS
fchr = 0
fpos0 = 1
fpos1 = 2
fid = 3
fallele1 = 4
fallele2 = 5
fmaf = 6

# CADD FIELDS
alt = 3
rawscore = 4
phredscore = 5

sys.stderr.write("Opening %s...\n"%(filename))
regionTabix = pysam.Tabixfile(filename,'r')

for line in sys.stdin:
  fields = line.rstrip().split('\t')
  chrom = fields[fchr].replace('chr','')
  # strip "chr"
  #chrom = chrom[3:]
  pos = int(fields[fpos1])
  alleles = set([fields[fallele1],fields[fallele2]])

  raw = []
  phred = []
  for allele in alleles:
    for CADDline in regionTabix.fetch(chrom,pos-1,pos):
      caddFields = CADDline.rstrip().split('\t')
      if (caddFields[alt] == allele):
        raw.append(float(caddFields[rawscore]))
        phred.append(float(caddFields[phredscore]))
        break
  
  # print max of scores or NA if scores is empty
  if not raw:
    sys.stdout.write('chr' + "\t".join(fields + ['NA','NA']) + '\n')
  else:
    maxRaw = max(raw) # need max in case of multiple alternate alleles
    maxPhred = max(phred)
    sys.stdout.write("\t".join(fields + [str(maxRaw),str(maxPhred)]) + '\n')


