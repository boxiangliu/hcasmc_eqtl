# Find HCASMC specific open chromatin
# Boxiang Liu
# 2017-12-21

# Setup:
# Preprocess encode data: 
Rscript hcasmc_specific_open_chromatin2/process_encode_data.R


# Find shared peaks: 
python hcasmc_specific_open_chromatin2/find_shared_peaks.py


# Plot Venn diagram:
Rscript hcasmc_specific_open_chromatin2/plot_venn_diagram.R
