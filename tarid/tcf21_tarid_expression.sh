# TARID
grep "ENSG00000227954" ../data/gtex/v6p/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct > ../processed_data/tarid/tarid_rpkm.gtex.txt
grep "ENSG00000227954" ../processed_data/160527/combined.filter.rpkm > ../processed_data/tarid/tarid_rpkm.hcasmc.txt


# TCF21:
grep "ENSG00000118526" ../data/gtex/v6p/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct > ../processed_data/tarid/tcf21_rpkm.gtex.txt
grep "ENSG00000118526" ../processed_data/160527/combined.filter.rpkm > ../processed_data/tarid/tcf21_rpkm.hcasmc.txt


# headers: 
grep "Name" ../data/gtex/v6p/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct > ../processed_data/tarid/header.gtex.txt
grep "Name" ../processed_data/160527/combined.filter.rpkm > ../processed_data/tarid/header.hcasmc.txt


# combine: 
paste <(cat ../processed_data/tarid/header.gtex.txt) <(cut -f2- ../processed_data/tarid/header.hcasmc.txt) > ../processed_data/tarid/header.txt
paste <(cat ../processed_data/tarid/tarid_rpkm.gtex.txt) <(cut -f2- ../processed_data/tarid/tarid_rpkm.hcasmc.txt) > ../processed_data/tarid/tarid_rpkm.txt
paste <(cat ../processed_data/tarid/tcf21_rpkm.gtex.txt) <(cut -f2- ../processed_data/tarid/tcf21_rpkm.hcasmc.txt) > ../processed_data/tarid/tcf21_rpkm.txt
cat ../processed_data/tarid/header.txt ../processed_data/tarid/tarid_rpkm.txt ../processed_data/tarid/tcf21_rpkm.txt > ../processed_data/tarid/tarid_tcf21_rpkm.txt
less -S ../processed_data/tarid/tarid_tcf21_rpkm.txt