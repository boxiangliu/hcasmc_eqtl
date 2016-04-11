################# SCG settings ################### 
# Job Name 
#$ -N Xplaceholder1
# 
# Request Large Memory Machine  
# -P large_mem
# 
# memory usage 
#$ -l h_vmem=4G 
#
# maximum run time 
# -l h_rt=168:00:00  
# 
# set queue:
#$ -q extended
#
# check for errors in the job submission options
#$ -w e
# 
# run on multiple threads                     
#$ -pe shm 8                                
#
# run job in current working directory      
#$ -cwd                                    
#
# set the output file
# -o  template.log
#
# merge stdout and stderr                                         
#$ -j y
# 
# email: 
#$ -m ea 
#$ -M jollier.liu@gmail.com
####################################################


picard=/srv/gsfs0/software/picard-tools/1.92
gatk=/srv/gsfs0/software/gatk/gatk-3.4.0/GenomeAnalysisTK.jar
hg19=/srv/gsfs0/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk_bundle=/srv/gsfs0/projects/montgomery/bliu2/shared/gatk_bundle_2.8_hg19
BWAIndex=/srv/gsfs0/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
cutadapt=/srv/gsfs0/software/cutadapt/1.9.dev1/bin/cutadapt
python=/srv/gsfs0/software/python/2.7/bin/python
bwa=/srv/gsfs0/software/bwa/bwa-0.7.12/bin/bwa
bgzip=/srv/gsfs0/software/samtools/samtools-1.2/bin/bgzip
tabix=/srv/gsfs0/software/samtools/samtools-1.2/bin/tabix
module load java

sample=placeholder1
reads1=placeholder2
reads2=placeholder3


# change wd: 
cd /srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/data/$sample

# trim adapters:
$python $cutadapt -m 30 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${reads1/gz/trimmed.gz} -p ${reads2/gz/trimmed.gz} $reads1 $reads2 
touch ${reads1/gz/trimmed.gz.done}
touch ${reads2/gz/trimmed.gz.done}

# align, sort and mark duplicates: 
$bwa mem -M -t 8 -R '@RG\tID:placeholder1\tSM:placeholder1\tPL:illumina\tLB:lib1\tPU:unit1' $BWAIndex ${reads1/gz/trimmed.gz} ${reads2/gz/trimmed.gz} > aln.sam 
touch aln.sam.done
rm ${reads1/gz/trimmed.gz}
rm ${reads2/gz/trimmed.gz}

java -Xmx2g -jar $picard/SortSam.jar INPUT=aln.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate 
touch sorted_reads.bam.done 
rm aln.sam

java -Xmx2g -jar $picard/MarkDuplicates.jar INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt
touch dedup_reads.bam.done 
rm sorted_reads.bam

java -Xmx2g -jar $picard/BuildBamIndex.jar INPUT=dedup_reads.bam
touch dedup_reads.bai.done 

# indel realignment: 
java -Xmx16g -jar $gatk -T RealignerTargetCreator -nt 8 -R $hg19 -I dedup_reads.bam -known $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -known $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o realignment_targets.list
touch realignment_targets.list.done

java -Xmx4g -jar $gatk -T IndelRealigner -R $hg19 -I dedup_reads.bam -targetIntervals realignment_targets.list -known $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -known $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o realigned_reads.bam 
touch realigned_reads.bam.done
rm dedup_reads.bam

# Base recalibration: 
java -Xmx4g -jar $gatk -T BaseRecalibrator -nct 8 -R $hg19 -I realigned_reads.bam -knownSites $gatk_bundle/dbsnp_138.hg19.vcf -knownSites $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -knownSites $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o recal_data.table
touch recal_data.table.done
# java -Xmx2g -jar $gatk -T BaseRecalibrator -R $hg19 -I realigned_reads.bam -knownSites $gatk_bundle/dbsnp_138.hg19.vcf -knownSites $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -knownSites $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -BQSR recal_data.table -o post_recal_data.table 
# java -Xmx2g -jar $gatk -T AnalyzeCovariates  -R $hg19 -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf
java -Xmx4g -jar $gatk -T PrintReads -nct 8 -R $hg19 -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam 
touch recal_reads.bam.done 
rm realigned_reads.bam

# Variant discovery (single sample):
## java -Xmx2g -jar $gatk -T HaplotypeCaller -R $hg19 -I recal_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o raw_variants.vcf 

# Joint variant discovery: 
java -Xmx16g -jar $gatk -T HaplotypeCaller -nct 16 -R $hg19 -I recal_reads.bam --genotyping_mode DISCOVERY --emitRefConfidence GVCF -o raw_variants.g.vcf 
touch raw_variants.g.vcf.done

$bgzip raw_variants.g.vcf
$tabix -p vcf raw_variants.g.vcf.gz 
