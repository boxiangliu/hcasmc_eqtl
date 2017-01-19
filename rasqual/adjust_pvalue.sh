awk '{print $26}' ../data/eQTL/rasqual/expressedGenes.pval.txt > ../data/eQTL/rasqual/expressedGenes.pval.tmp.txt
Rscript $scripts/rasqual/adjust_pvalue.R ../data/eQTL/rasqual/expressedGenes.pval.tmp.txt ../data/eQTL/rasqual/expressedGenes.padj.tmp.txt
paste ../data/eQTL/rasqual/expressedGenes.pval.txt <(cut -f2-3 ../data/eQTL/rasqual/expressedGenes.padj.tmp.txt) > ../data/eQTL/rasqual/expressedGenes.padj.txt
rm ../data/eQTL/rasqual/expressedGenes.pval.tmp.txt ../data/eQTL/rasqual/expressedGenes.padj.tmp.txt
