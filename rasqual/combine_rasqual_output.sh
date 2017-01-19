cat ../processed_data/160527/combined.filter.rpkm | cut -f1 | grep -v Name > ../data/eQTL/rasqual/expressed_genes.txt
n=0
head -n1 ../processed_data/rasqual/output/ENSG00000273493.1_RP11-80H18.4.pval.txt > ../data/eQTL/rasqual/expressedGenes.pval.txt
while read gene_id; do
	n=$((n+1))
	echo $n: $gene_id
	cat ../processed_data/rasqual/output/$gene_id*pval.txt | grep -v 'fid' >> ../data/eQTL/rasqual/expressedGenes.pval.txt
done < ../data/eQTL/rasqual/expressed_genes.txt
rm ../data/eQTL/rasqual/expressed_genes.txt