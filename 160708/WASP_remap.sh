# on valk:
cd /home/diskstation/RNAseq/HCASMC/Mapped_Files
ls -d HCASMC_RNASEQ_60lines_mapping_part1/*/ # copy and paste onto WASP_remap.sample_list.txt
ls -d HCASMC_RNASEQ_60lines_mapping_part2/*/ # copy and paste onto WASP_remap.sample_list.txt
sed -e "s/part1/part*/" -e "s/part2/part*/" 160708/WASP_remap.sample_list.txt | sort > 160708/WASP_remap.sample_list.2.txt
cp 160708/WASP_remap.sample_list.2.txt 160708/WASP_remap.sample_list.3.txt
# I manually corrected the correspondence between samples on VALK and those on DURGA

# convert excel to txt:
Rscript $scripts/160708/convert_excel_to_txt.R \
	-excel_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rna_wgs_match.reduced_050616.xlsx' \
	-sheet=5 \
	-subset='Path_to_alignment_on_valk','Path_to_alignment_on_durga' \
	-output_file='/srv/persistent/bliu2/HCASMC_eQTL/scripts/160708/WASP_remap.sample_list.4.txt'


# run STAR: 
bash $scripts/160708/WASP_remap.core.sh /srv/persistent/bliu2/HCASMC_eQTL/scripts/160708/WASP_remap.sample_list.4.txt
