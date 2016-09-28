cat ../data/gtex/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$14}' | sort -k1,1 > ../data/gtex/SAMPID_SMTSD.sorted.txt
cut -f3- ../processed_data/tarid/header.gtex.txt | bash $scripts/tarid/transpose.sh | sort > ../processed_data/tarid/header.gtex.t.sorted.txt
join -1 1 -2 1 -t $'\t' ../processed_data/tarid/header.gtex.t.sorted.txt ../data/gtex/SAMPID_SMTSD.sorted.txt > ../data/gtex/SAMPID_SMTSD.sorted.filtered.txt
cut -f2- ../processed_data/tarid/header.hcasmc.txt | bash $scripts/tarid/transpose.sh > ../data/gtex/tmp.hcasmc.txt
paste ../data/gtex/tmp.hcasmc.txt <(echo "HCASMC" | awk '{for (i=0;i<52;i++) print}') > ../data/gtex/tmp.hcasmc.2.txt
cat ../data/gtex/SAMPID_SMTSD.sorted.filtered.txt ../data/gtex/tmp.hcasmc.2.txt > ../data/gtex/SAMPID_SMTSD.sorted.filtered.with_hcasmc.txt
rm ../data/gtex/tmp.hcasmc.txt ../data/gtex/tmp.hcasmc.2.txt ../data/gtex/SAMPID_SMTSD.sorted.filtered.txt ../data/gtex/SAMPID_SMTSD.sorted.txt
