# analyzing conform gt result: 

# 10      69083   rs72635988      C       T       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE
# 10      95289   rs74113794      C       A       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE
# 10      95306   rs74651559      C       A       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE
# 10      98506   .       C       T       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE
# 10      111170  rs28376301      C       T       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE
# 10      111471  rs13328780      T       C       NOT_PERFORMED   NOT_PERFORMED   NOT_PERFORMED   REMOVED NOT_IN_REFERENCE

# 10:69083 is a quad-allelic loci and is filtered out in the 1kg p3 dataset 
# 10:95289 has maf of 0.026. But it is not in 1000G phase1 or phase3
# 10:95306 has maf of 0.026. But it is not in 1000G phase1 or phase3
# The last 2 variants seems to be too close to each other. 
# 10:98506 has maf of 0.008621. Not in 1000G phase1 or phase3
# 10:111170 is in recalibrated_variants.pass.norm.split.vcf.gz. 
# It is in 1000G phase1 v3 
# 10      111170  rs28376301      C       T       100     PASS    AC=160;AA=C;AN=2184;VT=SNP;ERATE=0.0204;THETA=0.0012;RSQ=0.9002;LDAF=0.0794;SNPSOURCE=LOWCOV;AVGP
# However, it is not in phase3v5

# 10:111471 has allele frequency of 0.216. But it is not in 1000G phase1 or phase3 