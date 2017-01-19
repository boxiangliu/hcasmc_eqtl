cat ../data/eQTL/fastqtl/fastqtl_nominal/fastqtl.allpairs.pc3.peer8.padj.txt | awk '{print $1,$2,$4,$5,$8}' > ../processed_data/compare_rasqual_and_fastqtl/fastqtl.tmp.txt
Rscript compare_rasqual_and_fastqtl/gene_id2name.R ../processed_data/compare_rasqual_and_fastqtl/fastqtl.tmp.txt ../processed_data/compare_rasqual_and_fastqtl/fastqtl.txt
cat ../data/eQTL/rasqual/expressedGenes.padj.txt | awk '{print $1,$3"_"$4"_"$5"_"$6,$12,$26,$29}' > ../processed_data/compare_rasqual_and_fastqtl/rasqual.txt
rm ../processed_data/compare_rasqual_and_fastqtl/fastqtl.tmp.txt