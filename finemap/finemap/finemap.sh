in_dir=../processed_data/finemap/finemap/howson_data/
out_dir=../processed_data/finemap/finemap/finemap/
[[ ! -d $out_dir ]] && mkdir -p $out_dir

input_fn=$out_dir/data
echo "z;ld;snp;config;log;n-ind" > $input_fn

for f in `ls $in_dir/*.z`;do
	echo $f

	n=`cat ${f/z/n}`
	echo -e "$f;${f/z/ld};${f/z/snp};${f/z/config};${f/z/log};$n" >> $input_fn
	
	sed 's/\t/ /g' ${f/z/ld} > tmp
	mv tmp ${f/z/ld}
done


for r in {1..1699}; do
/users/bliu2/tools/finemap_v1.1_x86_64/finemap --sss \
--in-files $input_fn --regions $r
done

mkdir -p $out_dir/default/
mv $in_dir/*.{snp,config} $out_dir/default/


for r in {1..1699}; do
/users/bliu2/tools/finemap_v1.1_x86_64/finemap --sss \
--in-files $input_fn --regions $r --n-causal-max 1
done

mkdir -p $out_dir/n_causal_max_1/
mv $in_dir/*.{snp,config} $out_dir/n_causal_max_1/


for r in {1..1699}; do
/users/bliu2/tools/finemap_v1.1_x86_64/finemap --sss \
--in-files $input_fn --regions $r --n-causal-max 2
done

mkdir -p $out_dir/n_causal_max_2/
mv $in_dir/*.{snp,config} $out_dir/n_causal_max_2/
