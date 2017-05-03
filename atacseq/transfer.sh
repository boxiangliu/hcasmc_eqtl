#------- FBS ---------# 
# Copy 2305 from valk:
scp bosh@valkyr:/home/clint/ATAC/150123_NS500418_0078_AH2JNYBGXX/Data/Intensities/BaseCalls/CA2305/atacseq/*FBS*concat*fastq ../data/atacseq/fbs/fastq/
# md5sum checked with no discrepencies. 


# gzip all fastq files:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/2305/fastq
gzip *fastq


# Copy 20805 from Valk: 
scp bosh@valkyr:/home/clint/ATAC/HCASMC/Hiseq/fastq/stimulation/020805.2_L1_TAAGGCGA_L001_R*.fastq.gz ../data/atacseq/fbs/fastq/


#------ SF ----------#
# prepare data for SF samples: 
mkdir -p ../data/atacseq/sf/fastq/
scp bosh@valkyr:/home/clint/ATAC/150123_NS500418_0078_AH2JNYBGXX/Data/Intensities/BaseCalls/*SF*concat*fastq ../data/atacseq/sf/fastq/
cd ../data/atacseq/sf/fastq/
for file in $(ls *fastq);do	echo $file; gzip $file &; done
wait 


