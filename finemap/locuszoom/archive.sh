#### fine-mapping:
# REST:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs17087335 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs17087335 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs66790703 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs66790703 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs56155140 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs56155140 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

# MAP3K7CL: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000156265.11_MAP3K7CL.pval.txt --pvalcol pval --markercol rsid --refsnp rs11911017 --chr chr21 --start 30500000 --end 30700000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000156265.11_MAP3K7CL)" --no-date --prefix ../figures/locuszoom/ENSG00000156265.11_MAP3K7CL_eqtl_chr21:30500000-30700000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11911017 --chr chr21 --start 30500000 --end 30700000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000156265.11_MAP3K7CL)" --no-date --prefix ../figures/locuszoom/ENSG00000156265.11_MAP3K7CL_gwas_chr21:30500000-30700000

# COL4A2: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000134871.13_COL4A2.pval.txt --pvalcol pval --markercol rsid --refsnp rs11838776 --flank 500kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_eqtl
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11838776 --flank 500kB --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_gwas

locuszoom --metal ../processed_data/rasqual/output/ENSG00000134871.13_COL4A2.pval.txt --pvalcol pval --markercol rsid --refsnp rs11838776 --flank 50kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_eqtl
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11838776 --flank 50kB --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_gwas
