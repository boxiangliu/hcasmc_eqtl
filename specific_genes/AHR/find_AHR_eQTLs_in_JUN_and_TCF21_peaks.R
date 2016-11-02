library(dplyr)
library(data.table)
library(stringr)

eqtl=fread('../processed_data/rasqual/output/ENSG00000106546.8_AHR.pval.txt')


eqtl[which(eqtl$rsid=='rs111472193'),]
eqtl[which(eqtl$rsid=='rs73075165'),]
eqtl[which(eqtl$rsid=='rs62444739'),]
# A TCF21 peak within chr7:17,569,211-17,570,989 contains two SNPs, 
# rs76720739 and rs62444739 were not tested as eQTL because of low MAF. 
# rs62444739 was not a significant (p~0.14) eQTL. 
# rs111472193 was not a significant eQTL (P<0.81). 
# rs73075165 is not a significant eQTL with p=0.9895647. 

eqtl[which(eqtl$rsid=='rs10265174'),]
eqtl[which(eqtl$rsid=='rs12699842'),]
eqtl[which(eqtl$rsid=='rs60457144'),]
# Another TCF21 and JUN peak within chr7:17,410,875-17,412,620 contains rs10265174, rs12699842, rs60457144. 
# rs10265174 is the second most significant eQTL with p-value < 9e-5. 
# rs10265174 also changes the binding affinity of AP-1 and TCF4 (which has similar motif to TCF21). 
# rs12699842 and rs60457144 are not tested due to low MAF. 


eqtl[which(eqtl$rsid=='rs187160691'),]
eqtl[which(eqtl$rsid=='rs149775688'),]
eqtl[which(eqtl$rsid=='rs140860934'),]
eqtl[which(eqtl$rsid=='rs1476080'),]
# Another TCF21 and JUN peak within chr7:17,357,504-17,358,786 contains 4 SNPs.
# rs187160691, rs149775688, and rs140860934 are not tested due to low MAF. 
# rs1476080, which is in JUN/JUND peak but not in TCF21 peak, is nominally significant with p-value 0.033. 

eqtl[which(eqtl$rsid=='rs6951212'),]
# Another TCF21 and JUN peak within chr7:17,344,243-17,345,487 contains several SNPs. 
# rs6951212 is in the center of TCF21 and edge of JUN peak; it is an eQTL for AHR with p<0.03. 
# Other variants are not tested due to low MAF. 


eqtl[which(eqtl$rsid=='rs17137508'),]
eqtl[which(eqtl$rsid=='rs73680651'),]
# Another TCF21 and JUN peak within chr7:17,314,154-17,316,566 contains 2 SNP. 
# However, neither rs17137508 nor rs73680651 are tested due to low MAF. 


eqtl[which(eqtl$rsid=='rs73077634'),]
eqtl[which(eqtl$rsid=='rs6963637'),]
# Another TCF21 peak (with weak JUN/JUND peak nearby) within chr7:17,246,408-17,247,665 contains 2 SNPs. 
# rs73077634 is not a significant eQTL with p < 0.5. 
# rs6963637 (within weak JUN peak and on the edge of TCF21 peak) is also not a significant eQTL with p < 0.3536105



eqtl[which(eqtl$rsid=='rs7804230')]
eqtl[which(eqtl$rsid=='rs137875887')]
eqtl[which(eqtl$rsid=='rs10229472')]
eqtl[which(eqtl$rsid=='rs10259316')]
eqtl[which(eqtl$rsid=='rs7800496')]
# Another JUN/JUND peak in chr7:17,444,757-17,446,288 contains 5 variants. 
# rs7804230 and rs10259316 are not significant eQTL both with p=0.1727379 (possibly due to LD). 
# rs10229472 and rs137875887 are not tested due to low MAF. 



eqtl[which(eqtl$rsid=='rs6966058')]
# Another JUN/JUND peak in chr7:17,447,959-17,448,718 contains 2 variants. 
# rs6966058 is a significant eQTL with p=0.03358077. rs6461315 is not a significant variant with p=0.8475344. 