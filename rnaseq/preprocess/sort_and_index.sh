cd ../data/rnaseq2/alignments/
# sort and index some bam and sam files: 
samtools sort -o 2305/Aligned.out.sorted.bam -O bam -@8 2305/Aligned.out.sam &
samtools sort -o 9070202/Aligned.out.sorted.bam -O bam -@8 9070202/Aligned.out.bam &
samtools sort -o 9052004/Aligned.out.sorted.bam -O bam -@8 9052004/Aligned.out.bam &
samtools sort -o 20805/Aligned.out.sorted.bam -O bam -@8 20805/Aligned.out.sam &
samtools index 2305/Aligned.out.sorted.bam & 
samtools index 9070202/Aligned.out.sorted.bam & 
samtools index 9052004/Aligned.out.sorted.bam & 
samtools index 20805/Aligned.out.sorted.bam & 
rm 2305/Aligned.out.sam 9070202/Aligned.out.bam 9052004/Aligned.out.bam 20805/Aligned.out.sam
