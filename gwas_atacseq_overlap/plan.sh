# download Roadmap data
# from /mnt/data/epigenomeRoadmap/peaks/unconsolidated/narrowPeak, find samples with DNase-seq 
# list of adult samples with DNA
# UW.Gastric.ChromatinAccessibility.STL001.DNase.DS20260.narrowPeak.gz 
# UW.Gastric.ChromatinAccessibility.STL003.DNase.DS20748.narrowPeak.gz E094
# UW.Heart.ChromatinAccessibility.STL001.DNase.DS20383.narrowPeak.gz
# UW.Ovary.ChromatinAccessibility.STL002.DNase.DS20827.narrowPeak.gz E097
# UW.Pancreas.ChromatinAccessibility.STL002.DNase.DS20842.narrowPeak.gz 
# UW.Pancreas.ChromatinAccessibility.STL003.DNase.DS20753.narrowPeak.gz E098
# UW.Penis_Foreskin_Melanocyte_Primary_Cells.ChromatinAccessibility.skin02.DNase.DS18668.narrowPeak.gz 
# UW.Penis_Foreskin_Melanocyte_Primary_Cells.ChromatinAccessibility.skin02.DNase.DS19662.narrowPeak.gz
# UW.Psoas_Muscle.ChromatinAccessibility.STL001.DNase.DS20325.narrowPeak.gz E100
# UW.Small_Intestine.ChromatinAccessibility.STL003.DNase.DS20770.narrowPeak.gz E109 

for file in $(ls UW.{Gastric,Heart,Overy,Pancreas,Penis,Psoas,Small}*DNase*); do
	zcat $file | wc -l 
done

# 580277
# 575295
# 419849
# 277424
# 473846
# 403709
# 485361
# 394073
# 356971

# The GWAS loci is located at ../processed_data/160615/compare_hcasmc_and_gtex/GWAS.txt
# The atacseq peak is located at ~/atacseq/2305/out/peaks

