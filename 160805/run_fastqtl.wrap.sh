vcf=$1
bed=$2
cov=$3
out=$4
args=$5

fastqtl=/usr/bin/fastQTL

n=0
for i in {1..22}; do 
	n=$(($n+1))
	if [[ n -gt 10 ]]; then
		wait
		n=0
	fi 

	echo "$fastqtl --vcf $vcf \
		--bed $bed \
		--out $out \
		--cov $cov \
		$args \
		--region chr$i"

	$fastqtl --vcf $vcf \
		--bed $bed \
		--out $out.chr$i \
		--cov $cov \
		$args \
		--region chr$i > $out.chr$i.log &

done 
wait
cat $out.chr{1..22} > $out
rm $out.chr{1..22}