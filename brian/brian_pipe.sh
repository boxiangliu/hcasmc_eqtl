# Setup: 
mkdir -p ../figures/brian/locuszoom brian ../processed_data/brian

# make locuszoom for AHR: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000106546.8_AHR.pval.txt --pvalcol pval --markercol rsid --refsnp rs10265174 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (AHR)" --prefix ../figures/brian/locuszoom/AHR_eQTL
locuszoom --metal ../processed_data/rasqual/output/ENSG00000106546.8_AHR.pval.txt --pvalcol pval --markercol rsid --refsnp rs10265174 --flank 500KB ymax=5 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (AHR)" --prefix ../figures/brian/locuszoom/AHR_eQTL

# MMP1:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1799750 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1799750 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1711437 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs751547 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1784405 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl

# CYP1B1:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000138061.7_CYP1B1.pval.txt --pvalcol pval --markercol rsid --refsnp rs4670837 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000138061.7_CYP1B1)" --prefix ../figures/brian/locuszoom/ENSG00000138061.7_CYP1B1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000138061.7_CYP1B1.pval.txt --pvalcol pval --markercol rsid --refsnp rs162558 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000138061.7_CYP1B1)" --prefix ../figures/brian/locuszoom/ENSG00000138061.7_CYP1B1_eqtl

