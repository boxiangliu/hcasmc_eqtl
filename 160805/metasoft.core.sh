# run metasoft
input=$1
output=$2
log=$3
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

