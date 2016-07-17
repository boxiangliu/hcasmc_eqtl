#!/bin/bash

# paths:
wd=$1
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
if [[ n -gt 12 ]];then
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
if [[ n -gt 12 ]];then
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
if [[ n -gt 12 ]];then
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
if [[ n -gt 12 ]];then
	wait
	n=0
fi 
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.bam.bai ]]; then 
	samtools index $sample/Aligned.out.sorted.rg.uniq.bam &
fi
done
wait 

# mark duplicates: 
MarkDuplicates=/software/picard-tools/1.92/MarkDuplicates.jar 
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 12 ]];then
	wait
	n=0
fi
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.dup.bam ]]; then 
java -Xmx2g -jar $MarkDuplicates I=$sample/Aligned.out.sorted.rg.uniq.bam O=$sample/Aligned.out.sorted.rg.uniq.dup.bam  M=$sample/marked_dup_metrics.txt 2> MarkDuplicates.${sample///}.log &
fi 
done
wait 

# index:
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 12 ]];then
	wait
	n=0
fi 
echo $sample
if [[ ! -f $sample/Aligned.out.sorted.rg.uniq.dup.bam.bai ]]; then 
samtools index $sample/Aligned.out.sorted.rg.uniq.dup.bam &
fi
done
wait 


# run RNA-seQC: 
rnaseqc=/srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar
gencode19=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf
hg19=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa
rRNA=/srv/persistent/bliu2/shared/genomes/rRNA/human_all_rRNA.fasta

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 12 ]];then
	wait
	n=0
fi

if [[ $(grep "Finished Successfully" rnaseqc.$sample.log) == "" ]]; then 
echo $sample
java -Xmx6g -jar $rnaseqc -n 1000 -s "$sample|$sample/Aligned.out.sorted.rg.uniq.dup.bam|$sample" -t $gencode19 -r $hg19 -o $sample/report -noDoC -strictMode > rnaseqc.$sample.log &
fi
done
