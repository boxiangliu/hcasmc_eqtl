mkdir -p ../data/joint3/asvcf_snpEff/
out_dir=../processed_data/sqtl/fastQTL/snpEff/
mkdir -p $out_dir

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
../data/joint3/asvcf/phased_and_imputed.chr{}.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz '|' \
bgzip '>' \
../data/joint3/asvcf_snpEff/chr{}.snpEff.gz ::: {1..22}

# Concatenate:
zcat ../data/joint3/asvcf_snpEff/chr{1..22}.snpEff.gz | \
bgzip > ../data/joint3/asvcf_snpEff/all.snpEff.gz

# Extract all annotations using SnpSift:
parallel -j10 \
extract_anno ../data/joint3/asvcf_snpEff/all.snpEff.gz \
{} \
$out_dir/{}.bed.gz \
::: downstream_gene_variant exon_variant intron_variant missense_variant \
splice_acceptor_variant splice_donor_variant splice_region_variant \
synonymous_variant upstream_gene_variant 3_prime_UTR_variant 5_prime_UTR_variant

# Concatenate:
zcat $out_dir/*_variant.bed.gz | bgzip > $out_dir/all.bed
