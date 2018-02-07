# Run Metasoft
# Boxiang Liu
# 2018-01-05

mkdir -p ../processed_data/eqtl/metasoft/rs2327429/metasoft_output/
input=../processed_data/eqtl/metasoft/rs2327429/metasoft_input/metasoft_input.txt
output=../processed_data/eqtl/metasoft/rs2327429/metasoft_output/metasoft_output.mcmc.txt
log=../processed_data/eqtl/metasoft/metasoft_output/rs2327429/metasoft_output.mcmc.log

metasoft=/srv/persistent/bliu2/tools/Metasoft/Metasoft.jar
pvalue_table=/srv/persistent/bliu2/tools/Metasoft/HanEskinPvalueTable.txt
echo [start] $(date)
java -jar $metasoft \
	-input $input \
	-output $output \
	-log $log \
	-pvalue_table $pvalue_table \
	-mvalue true \
	-mvalue_method mcmc \
	-mvalue_p_thres 1e0
echo [end] $(date)

