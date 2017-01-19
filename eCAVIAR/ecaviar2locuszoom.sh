col_file=$1
gwas_file=$2
eqtl_file=$3
tmp_dir=$4
out_dir=$5


# col_file=/srv/persistent/bliu2/HCASMC_eQTL/processed_data//eCAVIAR/eCAVIAR_output_rasqual/ENSG00000165895.13_ARHGAP42.ecaviar_col
# gwas_file=/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/cad.add.160614.website.txt
# eqtl_file=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output/ENSG00000165895.13_ARHGAP42.pval.txt 
# tmp_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data//eCAVIAR/tmp
# out_dir=../figures/locuszoom


fid=$(basename $col_file)
fid=${fid/.ecaviar_col/}
cat $col_file | awk '{if ($2>0.01) print $0}' > $tmp_dir/$fid.tmp

while read line; do
	split_line=($line)
	chr=$(echo ${split_line[0]} | cut -f1 -d_)
	pos=$(echo ${split_line[0]} | cut -f2 -d_)
	rsid=$(grep -P "$chr\t$pos" $gwas_file | cut -f1)
	echo ${split_line[0]} $chr $pos $rsid
 
	# make locuszoom for eqtl: 
	locuszoom --metal $eqtl_file --pvalcol pval --markercol rsid --refsnp $rsid --flank 1MB --source 1000G_March2012 --build hg19 ymax=5 --no-date --pop EUR title="eQTL ($fid)" --prefix $out_dir/${fid}_eqtl

	# make locuszoom for gwas: 
	locuszoom --metal $gwas_file --pvalcol p_dgc --markercol markername --refsnp $rsid --flank 1MB --source 1000G_March2012 --build hg19 ymax=5 --no-date --pop EUR title="GWAS ($fid)" --prefix $out_dir/${fid}_gwas
done < $tmp_dir/$fid.tmp