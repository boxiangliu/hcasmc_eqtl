dst_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments
sample_list=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments/sample_list.txt

while read dir; do
	echo $dir
	[[ ! -d $dst_dir/$dir ]] && mkdir -p $dst_dir/$dir
	rsync -avzh --progress bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase/$dir/Pass2/{Log.out,Log.final.out,Log.progress.out,SJ.out.tab,Aligned.out.sam.bam.sorted.bam} $dst_dir/$dir
done < $sample_list

touch .copy_rnaseq_from_valk.done