# Calculate the distance to closest splice donor site
# Boxiang Liu (bliu2@stanford.edu)


# Variables:
SNPsnap_fn=$1
fastQTL_dir=$2
out_dir=$3
SNPsnap_bed=${1/.tab/.dist_to_intron.bed}

mkdir -p $out_dir

SNPsnap_fn='/srv/persistent/bliu2/shared/SNPsnap/kb1000_collection.tab'

# Functions: 
SNPsnap2bed(){
	SNPsnap_fn=$1
	SNPsnap_bed=$2

	echo INFO SNPsnap2bed
	echo INFO - input: $SNPsnap_fn
	echo INFO - output: $SNPsnap_bed

	cut -f1,1 $SNPsnap_fn | \
	grep -v snpID | \
	awk 'BEGIN{FS=":";OFS="\t"}{print "chr"$1,$2-1,$2}' | \
	sort -k1,1 -k2,2n > \
	$SNPsnap_bed
}


fastQTL2bed(){
	fastQTL_fn=$1
	fastQTL_bed=$2

	echo INFO - fastQTL2bed
	echo INFO - input: $fastQTL_fn
	echo INFO - output: $fastQTL_bed

	zcat $fastQTL_fn | \
	cut -f1,1 | \
	awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2-1,$3,$0}' | \
	sort -k1,1 -k2,2n | \
	uniq > \
	$fastQTL_bed
}

find_closest(){
	SNPsnap_bed=$1
	fastQTL_bed=$2
	out_fn=$3

	echo INFO - find_closest
	echo INFO - input: $SNPsnap_bed
	echo INFO - input: $fastQTL_bed
	echo INFO - output: $out_fn

	bedtools closest \
	-a $SNPsnap_bed \
	-b $fastQTL_bed > \
	$out_fn
}

echo INFO - converting SNPsnp to bed format
SNPsnap2bed \
$SNPsnap_fn \
$out_dir/kb1000_collection.bed

for i in {1..22}; do
	echo INFO - chr$i

	echo INFO - subsetting SNPsnap 
	grep -P "^chr$i\t" \
	$out_dir/kb1000_collection.bed > \
	$out_dir/kb1000_collection.chr$i.bed

	echo INFO - converting fastQTL to bed
	fastQTL2bed \
	$fastQTL_dir/chr$i.nominal.txt.gz \
	$out_dir/chr$i.nominal.bed


	echo INFO - finding closest intron
	find_closest \
	$out_dir/kb1000_collection.chr$i.bed \
	$out_dir/chr$i.nominal.bed \
	$out_dir/chr$i.closest.bed

	echo INFO - calculating distance to closest intron
	Rscript sqtl/utils/dist_to_intron.R \
	$out_dir/chr$i.closest.bed \
	$out_dir/chr$i.dist_to_intron.txt \

done 

echo INFO - concatenating distance files 
cat \
$out_dir/chr*.dist_to_intron.txt > \
$out_dir/all.dist_to_intron.txt

sort -k1,1 $SNPsnap_fn > ${SNPsnap_fn/tab/sort.tab}
sort -k1,1 $out_dir/all.dist_to_intron.txt > $out_dir/all.dist_to_intron.sort.txt

echo INFO - joining distance to intron to SNPsnap 
join -1 1 -2 1 \
${SNPsnap_fn/tab/sort.tab} \
$out_dir/all.dist_to_intron.sort.txt \
-t $'\t' > \
$SNPsnap_bed


echo INFO - removing temp directory
# rm -r $out_dir
