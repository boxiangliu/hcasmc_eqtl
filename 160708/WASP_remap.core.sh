sample_list=$1
while read line; do
	# echo $line
	if [[ $line == "Path_to_alignment_on_valk"* ]]; then continue; fi
	line=($line)
	valk_dir=${line[0]}
	durga_dir=${line[1]}


	cd $durga_dir

	echo "rsync -vzh bosh@valkyr.stanford.edu:$valk_dir/GenomeForPass2/ $durga_dir"
	rsync -vzh bosh@valkyr.stanford.edu:$valk_dir/GenomeForPass2/* $durga_dir/GenomeForPass2/
	rsync -vzh bosh@valkyr.stanford.edu:$valk_dir/Pass2/Log.progress.out $durga_dir/GenomeForPass2/
	if ! cmp $durga_dir/GenomeForPass2/Log.progress.out $durga_dir/Log.progress.out >/dev/null 2>&1; then 
		echo "$durga_dir and $valk_dir differ" >> /srv/persistent/bliu2/HCASMC_eQTL/scripts/160708/WASP_remap.log
		continue 
	fi 
	CommonPars="--runThreadN 20 --outSAMattributes All --genomeLoad NoSharedMemory"
	Reads="wasp.remap.fq1.gz wasp.remap.fq2.gz --readFilesCommand zcat"
	echo "STAR $CommonPars --genomeDir GenomeForPass2 --readFilesIn $Reads --outFileNamePrefix wasp.remapped."
	STAR $CommonPars --genomeDir GenomeForPass2 --readFilesIn $Reads --outFileNamePrefix wasp.remapped. --outSAMtype BAM Unsorted
	rm -r GenomeForPass2

done < $sample_list
