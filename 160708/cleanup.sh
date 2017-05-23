mkdir /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/wasp
sample_dirs=($(ls -d /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/*/))
echo ${sample_dirs[*]}

for src_dir in ${sample_dirs[@]}; do
	dst_dir=${src_dir/alignments/wasp}
	echo creating $dst_dir
	mkdir -p $dst_dir

	echo "moving $src_dir/wasp.* to $dst_dir"
	ln $src_dir/wasp.* $dst_dir

	echo "removing $dst_dir/{wasp.bam,wasp.sort.bam,wasp.remapped.*,wasp.remap.fq*.gz,wasp.to.remap.*}"
	rm $dst_dir/{wasp.bam,wasp.sort.bam,wasp.remapped.*,wasp.remap.fq*.gz,wasp.to.remap.*}

done 