#!/bin/bash

# paths:
wd=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
cd $wd 

# make array of sample names: 
samples=($(ls -d */))
samples=(${samples[@]///}) # remove the "/"

# add read group:
AddOrReplaceReadGroups=/software/picard-tools/1.92/AddOrReplaceReadGroups.jar
n=0
for sample in ${samples[@]}; do 
echo $sample
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
if [[ ! -f $sample/Aligned.out.sorted.rg.bam ]]; then 
	java -Xmx2g -jar $AddOrReplaceReadGroups I=$sample/Aligned.out.sorted.bam O=$sample/Aligned.out.sorted.rg.bam RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> AddOrReplaceReadGroups.$sample.log &
fi 
done
wait


# index:
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.bam.bai ]]; then 
	samtools index $sample/Aligned.out.sorted.rg.bam &
fi
done
wait

# filter for uniquely mapped reads:
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 15 ]];then
	wait
	n=0
fi

echo $sample

if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.bam ]]; then 
	samtools view -h -q 255 -b -o $sample/Aligned.out.sorted.rg.uniq.bam $sample/Aligned.out.sorted.rg.bam &
fi 
done
wait

# index:
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi 
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.bam.bai ]]; then 
	samtools index $sample/Aligned.out.sorted.rg.uniq.bam &
fi
done

# mark duplicates: 
MarkDuplicates=/software/picard-tools/1.92/MarkDuplicates.jar 
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 15 ]];then
	wait
	n=0
fi
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.dup.bam ]]; then 
java -Xmx2g -jar $MarkDuplicates I=$sample/Aligned.out.sorted.rg.uniq.bam O=$sample/Aligned.out.sorted.rg.uniq.dup.bam  M=$sample/marked_dup_metrics.txt 2> MarkDuplicates.${sample///}.log &
fi 
done
wait 

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.dup.bam.bai ]]; then 
samtools index $sample/Aligned.out.sorted.rg.uniq.dup.bam &
fi
done
wait 

# make sample file for RNA-seQC:
ls -d */ > sample_file.tmp
cat sample_file.tmp | sed "s:/::" | awk 'BEGIN {OFS="\t"; print "Sample ID","Bam File","Notes"} {print $1,$1"/Aligned.out.sorted.rg.uniq.dup.bam",$1}' > sample_file.txt
rm sample_file.tmp


# get gtex v6p gencode gene models: 
# this annotation collapses exons into genes:
cd /srv/persistent/bliu2/shared/annotation/
mkdir gtex
ln -s /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz . 
gunzip -c gencode.v19.genes.v6p.patched_contigs.gtf.gz > gencode.v19.genes.v6p.patched_contigs.gtf
cat gencode.v19.genes.v6p.patched_contigs.gtf | sed -e "/^[^#]/ s/^/chr/" -e "s/MT/M/" > gencode.v19.genes.v6p.hg19.gtf


# run RNA-seQC: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
rnaseqc=/srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar
gencode19=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf
hg19=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa
rRNA=/srv/persistent/bliu2/shared/genomes/rRNA/human_all_rRNA.fasta

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi
sample=${sample///} # remove the backslash

if [[ $(grep "Finished Successfully" rnaseqc.$sample.log) == "" ]]; then 
echo $sample
java -Xmx6g -jar $rnaseqc -n 1000 -s "$sample|$sample/Aligned.out.sorted.rg.uniq.dup.bam|$sample" -t $gencode19 -r $hg19 -o $sample/report -noDoC -strictMode > rnaseqc.$sample.log &
fi
done
