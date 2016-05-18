#!/bin/bash 
awk 'BEGIN {FS="\t";OFS="\t"} {if ($7!="") print $1,$7}' ../processed_data/rna_wgs_match.reduced_050616.txt
awk 'BEGIN {FS="\t";OFS=" "} {if ($7!="") print $7,"\\"}' ../processed_data/rna_wgs_match.reduced_050616.txt


# create directories:
src=../data/rnaseq/alignments/
dst=../data/rnaseq2/alignments/
mkdir $dst
cd $dst
mkdir 1020301 \
102901 \
1042702 \
1051601 \
1060602 \
10705 \
112201 \
1278 \
1346 \
1347 \
1369 \
1386 \
1401 \
1448 \
1483 \
1522 \
1559 \
1576 \
1587 \
1596 \
177089 \
1795 \
200212 \
2030801 \
2040401 \
20805 \
2105 \
1508 \
2228 \
2282 \
2135 \
2356 \
2435 \
2463 \
2510 \
2139 \
2989 \
2999 \
3003 \
3100203 \
3101801 \
317155 \
59386145 \
59885590 \
7103002 \
8072501 \
8100901 \
9052004 \
# 9052004 \
9070202 \
9071501 \
9090701 \
2305


# create hard link from old to new files:
cd ../../../scripts
ln $src/1020301_26417_4_5_6_GTGAAA_4_5_GTGAAA/* $dst/1020301
ln $src/102901.8_26439_GATCAG/* $dst/102901
ln $src/1042702_26421_4_5_6_AGTTCC_4_5_AGTTCC/* $dst/1042702
ln $src/105106101_25499_1_2_3_ACTGAT_2_3_ACTGAT/* $dst/1051601
ln $src/1060602_26428_CAGATC/* $dst/1060602
ln $src/10705_25500_1_2_3_ATTCCT_2_3_ATTCCT/* $dst/10705
ln $src/122018_25497_1_2_3_CGTACG_2_3_CGTACG/* $dst/112201
ln $src/1278_25493_1_2_3_TAGCTT_2_3_TAGCTT/* $dst/1278
ln $src/1346_26418_4_5_6_GTCCGC_4_5_GTCCGC/* $dst/1346
ln $src/1347_25501_1_2_3_CGATGT_2_3_CGATGT/* $dst/1347
ln $src/1369_26422_4_5_6_AGTCAA_4_5_AGTCAA/* $dst/1369
ln $src/1386_25492_1_2_3_GATCAG_2_3_GATCAG/* $dst/1386
ln $src/1401.1_26437_TTAGGC/* $dst/1401
ln $src/1448_25507_1_2_3_AGTTCC_2_3_AGTTCC/* $dst/1448
ln $src/1483_26431_AGTTCC/* $dst/1483
ln $src/1522_26435_GTGAAA/* $dst/1522
ln $src/1559_25489_1_2_3_ATCACG_2_3_ATCACG/* $dst/1559
ln $src/1576_25504_1_2_3_CAGATC_2_3_CAGATC/* $dst/1576
ln $src/1587_26440_TAGCTT/* $dst/1587
ln $src/1596_25508_1_2_3_ATGTCA_2_3_ATGTCA/* $dst/1596
ln $src/177089_25506_1_2_3_AGTCAA_2_3_AGTCAA/* $dst/177089
ln $src/1795_26432_ATGTCA/* $dst/1795
ln $src/200212_26429_CTTGTA/* $dst/200212
ln $src/2030801_26423_4_5_6_CTTGTA_4_5_CTTGTA/* $dst/2030801
ln $src/2040401_26414_4_5_6_GAGTGG_4_5_GAGTGG/* $dst/2040401
# ln $src/20805.4_dase 20805
ln $src/2105_25491_1_2_3_ACTTGA_2_3_ACTTGA/* $dst/2105
ln $src/2108_26433_CCGTCC/* $dst/1508
ln $src/2228_26412_4_5_6_GTTTCG_4_5_GTTTCG/* $dst/2228
ln $src/2282_26411_4_5_6_GTGGCC_4_5_GTGGCC/* $dst/2282
ln $src/2315_26438_ACTTGA/* $dst/2135
ln $src/2356_26424_4_5_6_CAGATC_4_5_CAGATC/* $dst/2356
ln $src/2435_26408_4_5_6_GATCAG_4_5_GATCAG/* $dst/2435
ln $src/2463_25502_1_2_3_TGACCA_2_3_TGACCA/* $dst/2463
ln $src/2510_26419_4_5_6_CCGTCC_4_5_CCGTCC/* $dst/2510
ln $src/2913_26442_GCCAAT/* $dst/2139
ln $src/2989_26405_4_5_6_ATCACG_4_5_ATCACG/* $dst/2989
ln $src/2999_25495_1_2_3_GTGGCC_2_3_GTGGCC/* $dst/2999
ln $src/3003_25490_1_2_3_TTAGGC_2_3_TTAGGC/* $dst/3003
ln $src/3100203_26415_4_5_6_ACTGAT_4_5_ACTGAT/* $dst/3100203
ln $src/3101801.2_26434_GTCCGC/* $dst/3101801
ln $src/317155_26407_4_5_6_ACTTGA_4_5_ACTTGA/* $dst/317155
ln $src/59386143_26427_ACAGTG/* $dst/59386145
ln $src/59885590_26425_CGATGT_2/* $dst/59885590
ln $src/7103002_26416_4_5_6_ATTCCT_4_5_ATTCCT/* $dst/7103002
ln $src/8072501_26426_TGACCA/* $dst/8072501
ln $src/8100901_25503_1_2_3_ACAGTG_2_3_ACAGTG/* $dst/8100901
# ln $src/9052004_26406_4_5_6_TTAGGC_4_5_TTAGGC/* $dst/9052004
# ln $src/9052004.4_dase 9052004
# ln $src/90702_26413_4_5_6_CGTACG_4_5_CGTACG/* $dst/9070202
ln $src/9071501.8_26436_ATCACG/* $dst/9071501
ln $src/9090701_25498_1_2_3_GAGTGG_2_3_GAGTGG/* $dst/9090701
# ln $src/CA2305_FBS2_dase $dst/2305


# transfer 9052004, 9070202, 2305 and 20802 from valk: 
rsync -azvh bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/HCASMC/Mapped_Files/HCASMC_RNASEQ_60lines_mapping_part1/90702_26413_4_5_6_CGTACG_4_5_CGTACG/NEWPARAM/90702_26413_4_5_6_CGTACG_4_5_CGTACG/Pass2/{Aligned.out.bam,Log.final.out,Log.out,Log.progress.out} $dst/9070202 &
rsync -azvh bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/HCASMC/Mapped_Files/HCASMC_RNASEQ_60lines_mapping_part1/REMAPPING_3_SAMPLES_WITH_LOW_MAPPED_READS/NEW_PARAM/9052004_26406_4_5_6_TTAGGC_4_5_TTAGGC/Pass2/{Aligned.out.bam,Log.final.out,Log.out,Log.progress.out} $dst/9052004 &
rsync -azvh bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase/FBS2_S4_merged_R1_001.fastq.gz_FBS2_S4_merged_R2_001.fastq.gz/Pass2/{Aligned.out.sam,Log.final.out,Log.out,Log.progress.out} $dst/2305 &
rsync -azvh bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase/S7_run0002_lane5_index7_1.fastq.gz_S7_run0002_lane5_index7_2.fastq.gz/Pass2/{Aligned.out.sam,Log.final.out,Log.out,Log.progress.out} $dst/20805 &
