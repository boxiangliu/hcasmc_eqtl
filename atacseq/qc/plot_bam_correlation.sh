out_dir=../processed_data/atacseq/qc/plot_bam_correlation/
fig_dir=../figures/atacseq/qc/plot_bam_correlation/
mkdir -p $out_dir $fig_dir

bamfiles=$(ls ../data/atacseq/fbs/*/out/align/*/*{1346,1508,1522,200212,2305,2356,2510,2989}*.nodup.bam)
labels="1346 1508-rep1 1508-rep2 1522 200212 2305-rep1 2305-rep2 2305-rep3 2356 2510 2989"
out_fn=$out_dir/bam_coverage.npz
multiBamSummary bins --bamfiles $bamfiles --labels $labels -o $out_fn


fig_fn=$fig_dir/bam_coverage_spearman_correlation.pdf
plotCorrelation --corData $out_fn --corMethod spearman --whatToPlot heatmap --plotFile $fig_fn --plotNumbers --plotHeight 6 --plotWidth 7
