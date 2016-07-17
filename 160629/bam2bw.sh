# command args:
bam_dir=$1


# constants:
wd=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/


# global:
bam_dir=$wd/$bam_dir
bam=$bam_dir/Aligned.out.sorted.rg.uniq.bam


# Create output directory:
out_prefix=$(echo $bam | sed -e "s/alignments/bigwig/" -e "s/\.bam//")
echo $out_prefix
out_dir=$(dirname $out_prefix)
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi


# calculate read depth:
depth=$(samtools flagstat $bam | head -n1 | cut -d" " -f1)
echo depth: $depth


# call bam2bw: 
bash $scripts/160629/bam2bw.core.sh $bam 255 $depth $out_prefix
