library(data.table)

# Variable: 
out_dir='../processed_data/hcasmc_specific_open_chromatin/raw_peak_specificity/'
if (!dir.exists(out_dir)) dir.create(out_dir)

# Get HCASMC peaks (raw peaks, not summit extend peaks):
peaks=fread('zcat ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc/2305_ppr.IDR0.1.filt.narrowPeak.gz')


# Give meaningfulnames and remove uncessary columms:
setnames(peaks,c('V1','V2','V3','V10'),c('chr','start','end','summit'))
peaks[,c('V4','V5','V6','V7','V8','V9'):=NULL]


# Keep only autosomes: 
peaks=peaks[(chr!='chrX')&(chr!='chrY')&(chr!='chrM'),]


# Calculate summit extend summit boundary:
peaks[,c('start_se','end_se'):=list(start+summit-75,start+summit+75)]


# Get HCASMC specific peaks:
peak_specific=fread('../processed_data/hcasmc_specific_open_chromatin/peak_specificity_filt/HCASMC.bed')


# Assign specificity index to raw peaks: 
peak_ol_specific=merge(peaks,peak_specific,by.x=c('chr','start_se','end_se'),by.y=c('chr','start','end'))
peak_ol_specific=peak_ol_specific[,list(psi=min(psi)),by=c('chr','start','end')]


# Output raw peaks with specificity index:
fwrite(peak_ol_specific,sprintf('%s/HCASMC.bed',out_dir),sep='\t')