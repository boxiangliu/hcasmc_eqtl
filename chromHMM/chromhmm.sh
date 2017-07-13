wd=/srv/persistent/bliu2/tools/ChromHMM/
bam_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/
cell_mark_fn=chromHMM/cell_mark.txt
bin_bam_dir=../data/chromatin/chromHMM/bin_bam/
out_dir=../data/chromatin/chromHMM/out/
mkdir -p $bin_bam_dir $out_dir
state=10

java -mx3200M -jar $wd/ChromHMM.jar \
BinarizeBam $wd/CHROMSIZES/hg19.txt \
$bam_dir $cell_mark_fn \
$bin_bam_dir

java -mx3200M -jar $wd/ChromHMM.jar LearnModel \
-p 15 -printposterior -printstatebyline \
$bin_bam_dir $out_dir $state hg19