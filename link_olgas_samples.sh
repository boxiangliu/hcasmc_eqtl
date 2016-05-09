#!/bin/bash 
# bosh liu
# 2016/04/14
# durga
# make symbolic links of RNAseq bam files from dase/ to HCASMC_eQTL/

# Samples in dase/data/RNAseq_HCASMC/: 
src=/srv/persistent/bliu2/dase/data/RNAseq_HCASMC/alignment/orignal/
dst=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/

# declare -A sample_hash=(["1051601.3"]=mS12.Aligned.out.sorted.bam \
# 	["8072501.3"]=mS6.Aligned.out.sorted.bam \
# 	["9071501.8"]=mS14.Aligned.out.sorted.bam \
# 	["200212"]=S2_run0002_lane4_index4.Aligned.out.sorted.bam \
# 	["177089"]=S16_run0002_lane6_index7.Aligned.out.sorted.bam \
# 	["20805.4"]=S8_run0002_lane5_index12.Aligned.out.sorted.bam \
# 	["9090701.3"]=mS4.Aligned.out.sorted.bam \
# 	["8100901.2"]=mS10.Aligned.out.sorted.bam \
# 	["9052004.4"]=mS18.Aligned.out.sorted.bam)

declare -A sample_hash=(["1051601.3"]=S11_run0002_lane5_index4.Aligned.out.sorted.bam \
	["8072501.3"]=S5_run0002_lane4_index6.Aligned.out.sorted.bam \
	["9071501.8"]=S13_run0002_lane6_index5.Aligned.out.sorted.bam \
	["200212"]=S1_run0002_lane4_index2.Aligned.out.sorted.bam \
	["177089"]=S15_run0002_lane6_index6.Aligned.out.sorted.bam \
	["20805.4"]=S7_run0002_lane5_index7.Aligned.out.sorted.bam \
	["8100901.2"]=S9_run0002_lane5_index2.Aligned.out.sorted.bam \
	["9052004.4"]=pS17.Aligned.out.sorted.bam)

for key in ${!sample_hash[@]}; do
	value=${sample_hash[$key]}
	echo "sample: $key"
	echo "bam file: $value"
	[[ ! -d $dst/${key}_dase/ ]] && mkdir $dst/${key}_dase/
	ln $src/${value} $dst/${key}_dase/Aligned.out.sorted.bam
	ln $src/${value/Aligned.out.sorted.bam/Log.out} $dst/${key}_dase/Log.out
	ln $src/${value/Aligned.out.sorted.bam/Log.final.out} $dst/${key}_dase/Log.final.out	
	ln $src/${value/Aligned.out.sorted.bam/Log.progress.out} $dst/${key}_dase/Log.progress.out
done 


# samples in dase/data/RNAseq_CA2305_NextSeq_20150512/:
src=/srv/persistent/bliu2/dase/data/RNAseq_CA2305_NextSeq_20150512/alignment/merged/original
dst=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/
# declare -A sample_hash=(["CA2305_SF4"]=SF4_S1_merged_R1_001.sorted.bam \
# ["CA2305_SF5"]=SF5_S2_merged_R1_001.sorted.bam \
# ["CA2305_SF6"]=SF6_S3_merged_R1_001.sorted.bam)


declare -A sample_hash=(["CA2305_FBS2"]=FBS2_S4_merged_R1_001.sorted.bam \
["CA2305_FBS4"]=FBS4_S5_merged_R1_001.sorted.bam \
["CA2305_FBS6"]=FBS6_S6_merged_R1_001.sorted.bam)

for key in ${!sample_hash[@]}; do
	value=${sample_hash[$key]}
	echo "sample: $key"
	echo "bam file: $value"
	[[ ! -d $dst/${key}_dase/ ]] && mkdir $dst/${key}_dase/
	ln $src/${value} $dst/${key}_dase/Aligned.out.sorted.bam
done 

touch .link_olgas_samples.sh.done