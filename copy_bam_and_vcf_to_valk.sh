sample_list=$1
while read dir; do
	echo $dir
	# ssh bosh@valkyr.stanford.edu "[[ ! -d /home/diskstation/wgs/WGS_HCASMC2/$dir/ ]] && mkdir -p /home/diskstation/wgs/WGS_HCASMC2/$dir/"
	# # rsync -avzh --progress bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/HCASMC/Mapped_Files/$dir/Pass2/{Log.out,Log.final.out,Log.progress.out,Aligned.out.bam,Aligned.out.bam.bai} $dst_dir/$dir
	if [[ $dir != "joint/" ]]; then
		rsync -avzh --progress /mnt/data/WGS_HCASMC/$dir/{raw_variants.g.vcf.gz,raw_variants.g.vcf.gz.tbi,realignment_targets.list,recal_data.table,recal_reads.bai,recal_reads.bam} bosh@valkyr.stanford.edu:/home/diskstation/wgs/WGS_HCASMC2/$dir/
	fi
	if [[ $dir == "joint/" ]]; then
		rsync -avzh --progress /mnt/data/WGS_HCASMC/$dir/{raw_variants.vcf,raw_variants.vcf.idx,recalibrated_variants.vcf.gz,recalibrated_variants.vcf.gz.tbi} bosh@valkyr.stanford.edu:/home/diskstation/wgs/WGS_HCASMC2/$dir/
	fi
done < $sample_list

touch .copy_bam_and_vcf_to_valk.sh.done