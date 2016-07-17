# merge HCASMC read counts:

# variables:
data=/srv/persistent/bliu2/HCASMC_eQTL/data/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
wd=$data/rnaseq2/alignments/


# main: 
tail -n +3 $wd/1020301/report/genes.reads.gct | cut -f1-2 > $processed_data/160715/combined.count

while read sample;do
	sample=${sample///}
	echo $sample
	tail -n +3 $wd/$sample/report/genes.reads.gct | cut -f3 > $processed_data/160715/$sample.tmp
	cp $processed_data/160715/combined.count $processed_data/160715/combined.tmp
	paste -d "\t" $processed_data/160715/combined.tmp $processed_data/160715/$sample.tmp > $processed_data/160715/combined.count
	rm $processed_data/160715/$sample.tmp
done < $wd/sample_list.txt
rm $processed_data/160715/combined.tmp