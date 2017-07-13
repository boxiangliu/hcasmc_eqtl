howson_fn=../data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt
nikpay_fn=../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt


# TCF21:
grep rs2327429 $howson_fn
grep rs2327429 $nikpay_fn
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'
grep rs2327429 $eqtl_fn
atac_dir='../processed_data//atacseq/rasqual/output/chr6/'
grep 6_134209837_T_C $atac_dir/* | sort -k11nr | head -n1 | awk '{print $11}' 

grep 6:134214525 $howson_fn
grep rs12190287 $nikpay_fn
eqtl_fn='../processed_data/rasqual/output_pval/chr6/ENSG00000118526.6_TCF21.pval.txt'
grep rs12190287 $eqtl_fn
atac_dir='../processed_data//atacseq/rasqual/output/chr6/'
grep 6_134214525_C_G $atac_dir/* | sort -k11nr | head -n1 | awk '{print $11}' 
