# convert bam to junction:
bam_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/wasp/
if [[ -f $bam_dir/juncfiles.txt ]]; then rm $bam_dir/juncfiles.txt; fi
n=0
for bamfile in `ls $bam_dir/*/wasp.keep.merged.sorted.bam`; do
	n=$((n+1))
	if [[ n -ge 10 ]]; then 
		wait
		n=0
	fi
	juncfile=$(echo $bamfile | sed -e "s/wasp/leafcutter_wasp/" -e "s/\/wasp\.keep\.merged\.sorted\.bam//").junc
	echo Converting $bamfile to $juncfile
	sh /srv/persistent/bliu2/tools/leafcutter/scripts/bam2junc.sh $bamfile $juncfile &
	echo $juncfile >> $bam_dir/juncfiles.txt
done
wait