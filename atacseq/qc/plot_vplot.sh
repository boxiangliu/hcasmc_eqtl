fn=../data/atacseq/fbs/2989/out/signal/macs2/rep1/CA2989_S7_concat_R1_001.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig
out_dir=../processed_data/atacseq/qc/plot_vplot/
fig_dir=../figures/atacseq/qc/plot_vplot/
mkdir -p $out_dir $fig_dir

extract_sample_name(){
	fn=$1
	sample_name=$(basename $fn)
	sample_name=$(echo $sample_name | cut -d'_' -f1)
	sample_name=$(echo $sample_name | cut -d'-' -f1)
	sample_name=$(echo $sample_name | sed 's/CA//')
	replicate=$(dirname $fn)
	replicate=$(basename $replicate)
	sample_name=${sample_name}_${replicate}
	echo $sample_name
}
export -f extract_sample_name

for fn in $(ls ../data/atacseq/fbs/*/out/signal/macs2/rep*/*.fc.signal.bigwig); do

sample_name=$(extract_sample_name $fn)
echo $sample_name

computeMatrix reference-point -S $fn \
-R atacseq/qc/hg19.refGeneReviewedValidated.tss.chr4.bed \
-o $out_dir/$sample_name.matrix.gz \
--referencePoint TSS \
--binSize 10 --missingDataAsZero -b 2000 -a 2000

plotHeatmap -m $out_dir/$sample_name.matrix.gz \
-out $fig_dir/${sample_name}_heatmap.png \
--xAxisLabel 'Distance (bp)' \
--samplesLabel Insertions \
--zMin 0 -z ATAC 

plotProfile -m $out_dir/$sample_name.matrix.gz \
-out $fig_dir/${sample_name}_profile.png \
--outFileNameData $out_dir/$sample_name.profile.txt \
--samplesLabel Insertions -z ATAC

done
