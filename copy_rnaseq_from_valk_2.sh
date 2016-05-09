dst_dir=$1
sample_list=$2
while read dir; do
	echo $dir
	[[ ! -d $dst_dir/$dir ]] && mkdir -p $dst_dir/$dir
	rsync -avzh --progress bosh@valkyr.stanford.edu:/home/mpjanic/tmp_rnaseq/$dir/Pass2/{Log.out,Log.final.out,Log.progress.out,Aligned.out.bam.sort.bam,Aligned.out.bam.sort.bam.bai} $dst_dir/$dir
done < $sample_list

touch .copy_rnaseq_from_valk_2.done