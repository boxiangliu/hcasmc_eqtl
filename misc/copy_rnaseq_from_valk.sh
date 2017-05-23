dst_dir=$1
sample_list=$2
while read dir; do
	echo $dir
	[[ ! -d $dst_dir/$dir ]] && mkdir -p $dst_dir/$dir
	rsync -avzh --progress bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/HCASMC/Mapped_Files/$dir/Pass2/{Log.out,Log.final.out,Log.progress.out,Aligned.out.bam,Aligned.out.bam.bai} $dst_dir/$dir
done < $sample_list

touch .copy_rnaseq_from_valk.done