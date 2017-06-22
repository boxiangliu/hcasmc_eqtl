in_dir=../processed_data/finemap/finemap/rasqual_data/input/
out_dir=../processed_data/finemap/finemap/rasqual_finemap/
[[ ! -d $out_dir ]] && mkdir -p $out_dir

# Functions: 
tab2space(){
	f=$1
	sed 's/\t/ /g' $f > ${f/ld/tmp}; mv ${f/ld/tmp} $f
}

export -f tab2space

input_fn=$out_dir/data
echo "z;ld;snp;config;log;n-ind" > $input_fn

for f in `ls $in_dir/*.z`;do
	echo $f

	n=52
	echo -e "$f;${f/z/ld};${f/z/snp};${f/z/config};${f/z/log};$n" >> $input_fn
done


parallel -j30 tab2space {} ::: `ls $in_dir/*.ld`


for r in {1..211}; do
/users/bliu2/tools/finemap_v1.1_x86_64/finemap --sss \
--in-files $input_fn --regions $r --n-causal-max 1
done

mkdir -p $out_dir/n_causal_max_1/
mv $in_dir/*.{snp,config} $out_dir/n_causal_max_1/