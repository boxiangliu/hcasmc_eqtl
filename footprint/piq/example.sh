piq=/srv/persistent/bliu2/tools/piq-single/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/footprint

Rscript $piq/pwmmatch.exact.r $piq/common.r $piq/pwms/jasparfix.txt 139 $processed_data/piq/motif.matches/
Rscript bam2rdata.r common.r $processed_data/piq/bam/fbs_rep1.RData /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam
Rscript pertf.r common.r $processed_data/piq/motif.matches/ $processed_data/piq/tmp_fbs_rep1/ $processed_data/piq/res/ $processed_data/piq/bam/fbs_rep1.RData 139


Rscript $piq/pwmmatch.exact.r $piq/common.r $piq/pwms/jasparfix.txt 139 $processed_data/piq/motif.matches/
Rscript bam2rdata.r common.r $processed_data/piq/bam/fbs_rep1_masked.RData /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam
Rscript pertf.r common.r $processed_data/piq/motif.matches/ $processed_data/piq/tmp_fbs_rep1_masked/ $processed_data/piq/res_masked/ $processed_data/piq/bam/fbs_rep1_masked.RData 139


Rscript $piq/pwmmatch.exact.r $piq/common.r $piq/pwms/jasparfix.txt 139 $processed_data/piq/motif.matches/
Rscript pairedbam2rdata.r common.r $processed_data/piq/bam/fbs_rep1_masked_pe.RData /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam # slow
Rscript pertf.r common.r $processed_data/piq/motif.matches/ $processed_data/piq/tmp_fbs_rep1_masked_pe/ $processed_data/piq/res_masked_pe/ $processed_data/piq/bam/fbs_rep1_masked_pe.RData 139


Rscript $piq/pwmmatch.exact.r $piq/common.r $piq/pwms/jasparfix.txt 139 $processed_data/piq/motif.matches/
Rscript bam2rdata.r common.r $processed_data/piq/bam/fbs_rep1_masked_nodup.RData /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.nodup.bam
Rscript pertf.r common.r $processed_data/piq/motif.matches/ $processed_data/piq/tmp_fbs_rep1_masked_nodup/ $processed_data/piq/res_masked_nodup/ $processed_data/piq/bam/fbs_rep1_masked_nodup.RData 139
