#!/bin/bash
dst=../data/rnaseq2/alignments/
cd $dst
sample=9052004

# sort 
samtools sort -o $sample/Aligned.out.sorted.bam -O bam -@8 $sample/Aligned.out.sam 
samtools index $sample/Aligned.out.sorted.bam 
# rm $sample/Aligned.out.sam

# check bam file integrity:
samtools quickcheck $sample/Aligned.out.sorted.bam

# add read group: 
AddOrReplaceReadGroups=/software/picard-tools/1.92/AddOrReplaceReadGroups.jar
java -Xmx2g -jar $AddOrReplaceReadGroups I=$sample/Aligned.out.sorted.bam O=$sample/Aligned.out.sorted.rg.bam RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> AddOrReplaceReadGroups.$sample.log 
samtools index $sample/Aligned.out.sorted.rg.bam 

# select uniquely mapped reads:
samtools view -h -q 255 -b -o $sample/Aligned.out.sorted.rg.uniq.bam $sample/Aligned.out.sorted.rg.bam 
samtools index $sample/Aligned.out.sorted.rg.uniq.bam


# mark duplicates:
MarkDuplicates=/software/picard-tools/1.92/MarkDuplicates.jar 
java -Xmx2g -jar $MarkDuplicates I=$sample/Aligned.out.sorted.rg.uniq.bam O=$sample/Aligned.out.sorted.rg.uniq.dup.bam  M=$sample/marked_dup_metrics.txt 2> MarkDuplicates.$sample.log 
samtools index $sample/Aligned.out.sorted.rg.uniq.dup.bam 


# calculate rnaseqc:
rnaseqc=/srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar
gencode19=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf
hg19=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa
rRNA=/srv/persistent/bliu2/shared/genomes/rRNA/human_all_rRNA.fasta
java -jar $rnaseqc -n 1000 -s "$sample|$sample/Aligned.out.sorted.rg.uniq.dup.bam|$sample" -t $gencode19 -r $hg19 -o $sample/report -noDoC -strictMode > rnaseqc.$sample.log 