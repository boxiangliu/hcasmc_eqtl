# command args:
IN_FILE=$1
MAPQ=$2
DEPTH=$3
OUT_PREFIX=$4
# IN_FILE='../data/rnaseq2/alignments/1020301/Aligned.out.sorted.rg.uniq.bam'
# MAPQ=255
# DEPTH=1
# OUT_PREFIX='../data/rnaseq2/bigwig/1020301/1020301'


# constants: 
CHROMS='/srv/persistent/bliu2/shared/genomes/hg19.chrom.sizes'


# variables:
multiplier=`echo "(1 / $DEPTH) * 1000000" | bc -l`
bedgraph=$OUT_PREFIX.begGraph
bigwig=$OUT_PREFIX.bw


# generate bedGraph
samtools view -bh -q $MAPQ $IN_FILE | \
genomeCoverageBed -ibam stdin -g $CHROMS -bga -split -scale $multiplier \
> $bedgraph


#BedGraph -> BigWig
/software/ucsc_tools/3.0.9/bedGraphToBigWig $bedgraph $CHROMS $bigwig
