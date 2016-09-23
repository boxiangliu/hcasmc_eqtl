# count the number of egenes that passes FDR 0.05 for all GTEx tissue:
echo -e "tissue\tnum_egenes" > ../processed_data/egenes_vs_sample_size/num_egenes_subsampled_52.txt
for tissue in `cat /srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gtex.v6p.eqtl.tissues.txt`;do
	n_egenes=$(zcat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/$tissue/${tissue}_52.egenes.txt.gz | awk 'BEGIN{n=0} {if ($(NF-1)<0.05) n=n+1} END{print n}')
	echo -e "$tissue\t$n_egenes" >> ../processed_data/egenes_vs_sample_size/num_egenes_subsampled_52.txt
done

# count the number of egenes for HCASMC: 
n_egenes=$(cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data//160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt | awk 'BEGIN{n=0} {if ($(NF)<0.05) n=n+1} END{print n}')
echo -e "HCASMC\t$n_egenes" >> ../processed_data/egenes_vs_sample_size/num_egenes_subsampled_52.txt