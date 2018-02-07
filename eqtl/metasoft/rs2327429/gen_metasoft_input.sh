# Concatenate eQTL files 
# Boxiang Liu
# 2018-01-05

in_dir=../processed_data/eqtl/metasoft/rs2327429/rs2327429/
out_dir=../processed_data/eqtl/metasoft/rs2327429/metasoft_input/
mkdir -p $out_dir

cat $in_dir/*.rs2327429.txt | \
awk 'BEGIN{OFS="\t"}{print $1"_"$2,$4,$5,$6,$7}' | \
sort -k1,1 -k5,5 -V > \
$out_dir/All_Tissues.rs2327429.sorted.txt

cat $out_dir/All_Tissues.rs2327429.sorted.txt | \
python eqtl/metasoft/gen_metasoft_input.py 1 > \
$out_dir/metasoft_input.txt
