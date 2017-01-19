scp bosh@valkyr:"/home/towerraid/tomq/Sequencing/Sequences/TQ\ 21-23/*" /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/lane1.fastq /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k4me1.fastq
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/lane2.fastq /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k4me3.fastq
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/lane3.fastq /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27me3.fastq
echo -e "The following three are from bosh@valkyr:/home/towerraid/tomq/Sequencing/Sequences/TQ\ 21-23/. They are all serum-fed and from sample 1508.\nh3k4me1 = TQ21\nh3k4me3 = TQ22\nh3k27me3 = TQ23" > /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/README
bgzip /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k4me1.fastq
bgzip /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k4me3.fastq
bgzip /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27me3.fastq

scp bosh@valkyr:/home/clint/ChIPseq/H3K27ac/150513_PINKERTON_0360_BC7489ACXX_L4_ATCACG_*_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/
scp bosh@valkyr:/home/clint/ChIPseq/H3K27ac/150513_PINKERTON_0360_BC7489ACXX_L4_CAGATC_*_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/150513_PINKERTON_0360_BC7489ACXX_L4_ATCACG_1_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27ac_1.fastq.gz
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/150513_PINKERTON_0360_BC7489ACXX_L4_ATCACG_2_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27ac_2.fastq.gz
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/150513_PINKERTON_0360_BC7489ACXX_L4_CAGATC_1_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27ac_IgG_1.fastq.gz
mv /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/150513_PINKERTON_0360_BC7489ACXX_L4_CAGATC_2_pf.fastq.gz /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/h3k27ac_IgG_2.fastq.gz
echo -e "The following 2 samples are from bosh@valkyr:/home/clint/ChIPseq/H3K27ac/; they are all serum-fed and from sample 2108\nh3k27ac = ATCACG\nh3k27ac = CAGATC" >> /srv/persistent/bliu2/HCASMC_eQTL/data/chromatin/fastq/README