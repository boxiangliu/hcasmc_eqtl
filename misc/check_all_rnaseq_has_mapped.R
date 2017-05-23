#!/usr/bin/env Rscript 
# bosh liu
# 2016/04/12
# durga

all_samples_filename = '../processed_data/check_all_rnaseq_has_mapped/Stanford_3rd_round_AAA-StudyInfo.txt'
mapped_samples_filename = '../processed_data/rna_wgs_match.tsv'

all_samples = fread(all_samples_filename)
temp = all_samples[,.(SAMP_id,SAMP_name,RULA_Barcode)]
temp = unique(temp)
all_samples = temp

mapped_samples = fread(mapped_samples_filename)
split_names = str_split(mapped_samples$rna,"_")
rna_names = c()
for (split_name in split_names){
	print(split_name)
	rna_names = c(rna_names, split_name[1])
}
mapped_samples$rna_name = rna_names

temp = merge(mapped_samples, all_samples, by.x = 'rna_name', by.y = 'SAMP_name', all = T)
merged = temp
missing = sort(unique(which(is.na(merged), arr.ind=T)[,1]))
print(merged[missing,])

# sample 2913 is not mapped, details see below: 

#              rna_name                                   rna     dna SAMP_id
#  1:            1020301 1020301_26417_4_5_6_GTGAAA_4_5_GTGAAA 1020301      NA
#  2:          1020301.7                                    NA      NA   26417
#  3:       106(1060602)                                    NA      NA   26428
#  4:            1060602                  1060602_26428_CAGATC 1060602      NA
#  5:            2030801 2030801_26423_4_5_6_CTTGTA_4_5_CTTGTA 2030801      NA
#  6:          2030801.2                                    NA      NA   26423
#  7:               2913                                    NA      NA   26442
#  8:       310(3100203)                                    NA      NA   26415
#  9:            3100203 3100203_26415_4_5_6_ACTGAT_4_5_ACTGAT 3100203      NA
# 10:            8072501                                    NA      NA   26426
# 11: 9071501.8000000007                                    NA      NA   26436

