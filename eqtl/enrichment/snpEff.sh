in_dir=$1
out_dir=$2
mkdir -p $out_dir
# out_dir=../processed_data/sqtl/fastQTL/snpEff/

# Functions: 
extract_anno(){

	in_fn=$1
	anno=$2
	out_fn=$3

	echo INFO - input: $in_fn
	echo INFO - anno: $anno
	echo INFO - output: $out_fn

	java -Xmx4g \
	-jar /srv/persistent/bliu2/tools/snpEff/SnpSift.jar \
	filter "ANN[ANY].EFFECT has '$anno'" \
	$in_fn | \
	grep -v ^# | \
	awk -v anno=$anno 'BEGIN{OFS="\t"}{print $1,$2,$2,anno}' | \
	bgzip > $out_fn

}

export -f extract_anno


# Annotate using SnpEff:
parallel -j10 \
java -Xmx4g \
-jar /srv/persistent/bliu2/tools/snpEff/snpEff.jar \
hg19 \
$in_dir/chr{}.vcf.gz '|' \
bgzip '>' \
$out_dir/chr{}.snpEff.gz ::: {1..22}

# Concatenate:
zcat $out_dir/chr{1..22}.snpEff.gz | \
bgzip > $out_dir/all.snpEff.gz

# Extract all annotations using SnpSift:
mkdir -p $out_dir/annotation/

parallel -j10 \
extract_anno $out_dir/all.snpEff.gz \
{} \
$out_dir/annotation/{}.bed.gz \
::: downstream_gene_variant exon_variant intron_variant missense_variant \
splice_acceptor_variant splice_donor_variant splice_region_variant \
synonymous_variant upstream_gene_variant 3_prime_UTR_variant 5_prime_UTR_variant

# Concatenate:
zcat $out_dir/annotation/*_variant.bed.gz | bgzip > $out_dir/annotation/all.bed.gz
