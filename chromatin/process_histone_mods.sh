# Process histone modification data
# Need to run on Nandi!
set -e 
## Variables:
pipeline_dir=/users/bliu2/TF_chipseq_pipeline/
in_dir=/srv/scratch/bliu2/HCASMC_eQTL/chromatin/process_histone_mods/fastq
out_dir=/srv/scratch/bliu2/HCASMC_eQTL/chromatin/process_histone_mods/out
if [[ ! -d $in_dir ]]; then mkdir -p $in_dir; fi
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

## Main: 
# Move data from Durga to Nandi: 
if [[ ! -e $in_dir/h3k4me1.fastq.gz ]]; then 
	scp bliu2@durga:/srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/* $in_dir
fi 


# Process H3K4me1: 
[[ ! -d $out_dir/h3k4me1 ]] && mkdir $out_dir/h3k4me1/
bds $pipeline_dir/chipseq.bds -histone -species hg19 -nth 4 -fastq1 $in_dir/h3k4me1.fastq.gz -out_dir $out_dir/h3k4me1/


# Process H3K4me3: 
[[ ! -d $out_dir/h3k4me3 ]] && mkdir $out_dir/h3k4me3/
bds $pipeline_dir/chipseq.bds -histone -species hg19 -nth 4 -fastq1 $in_dir/h3k4me3.fastq.gz -out_dir $out_dir/h3k4me3/


# Process H3K27me3: 
[[ ! -d $out_dir/h3k27me3 ]] && mkdir $out_dir/h3k27me3/
bds $pipeline_dir/chipseq.bds -histone -species hg19 -nth 4 -fastq1 $in_dir/h3k27me3.fastq.gz -out_dir $out_dir/h3k27me3/


# Process H3K27ac: 
[[ ! -d $out_dir/h3k27ac ]] && mkdir $out_dir/h3k27ac/
bds $pipeline_dir/chipseq.bds -species hg19 -nth 4 -fastq1_1 $in_dir/h3k27ac_1.fastq.gz -fastq1_2 $in_dir/h3k27ac_2.fastq.gz -ctl_fastq1_1 $in_dir/h3k27ac_IgG_1.fastq.gz -ctl_fastq1_2 $in_dir/h3k27ac_IgG_2.fastq.gz -out_dir $out_dir/h3k27ac/
