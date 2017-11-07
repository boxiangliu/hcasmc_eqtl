# Setup:
mkdir -p ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/ ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc
ln /srv/persistent/bliu2/HCASMC_eQTL/processed_data/mpra/DHS_expanded/* ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc
ln ../data/atacseq/fbs/*/out/peak/macs2/overlap/*.naive_overlap.narrowPeak.gz ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc


# Preprocess encode data: 
Rscript hcasmc_specific_open_chromatin/process_encode_data.R
Rscript hcasmc_specific_open_chromatin/process_encode_data_without_taking_max.R


# Find shared peaks: 
python hcasmc_specific_open_chromatin/find_shared_peaks.py


# Get signal at shared peaks:
python hcasmc_specific_open_chromatin/get_signal_at_shared_peaks.py
Rscript hcasmc_specific_open_chromatin/merge_signals.R


# Quantile normalization: 
Rscript hcasmc_specific_open_chromatin/quantile_normalization.R


# Define HCASMC-specific peak score:
Rscript hcasmc_specific_open_chromatin/calc_peak_specificity_index.R


# Intersect specificity with sample peaks: 
Rscript hcasmc_specific_open_chromatin/intersect_specificity_and_peaks.R


# Take the minimum specificity index of each peak: 
Rscript hcasmc_specific_open_chromatin/get_min_index.R


# Assign peak specificity score to raw peaks:
Rscript hcasmc_specific_open_chromatin/assign_raw_peak_specificity.R

# Plot the distribution of number of tissues sharing each peak in HCASMC:
Rscript hcasmc_specific_open_chromatin/plot_distribution_of_tissue_sharing.R


# Perform multi-dimensional scaling: 
Rscript hcasmc_specific_open_chromatin/mds_open_chromatin.R 


# Classifying peaks based on TSS, gene body, intergenic:
Rscript hcasmc_specific_open_chromatin/classify_peaks.R


# Find specific open chromatin around TCF21:
Rscript hcasmc_specific_open_chromatin/peak_at_TCF21.R

