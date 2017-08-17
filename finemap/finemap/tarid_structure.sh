hg19=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa
rm_dir=/users/bliu2/tools/RMDetect/rmdetect_0.0.3/
out_dir=../processed_data/finemap/finemap/tarid_structure/
fig_dir=../figures/finemap/finemap/tarid_structure/

mkdir -p $out_dir $fig_dir

# get TARID sequence: 
samtools faidx $hg19 chr6:134209742-134210350 | \
grep -v '>' | tr -d '\n' > $out_dir/tarid_exon2.seq

# Predict motifs:
$rm_dir/rmdetect.py --both-strands $out_dir/tarid_exon2.seq > $out_dir/tarid_exon2_result.txt
$rm_dir/rmout.py < $out_dir/tarid_exon2_result.txt
$rm_dir/rmout.py --out=$fig_dir/tarid_exon2.pdf --cands=5 fold < $out_dir/tarid_exon2_result.txt


# get TARID alternative sequence:
sed 's/AGGTGTTGG/AGGTGCTGG/' $out_dir/tarid_exon2.seq > $out_dir/tarid_exon2.alt.seq

# Predict motifs: 
$rm_dir/rmdetect.py --both-strands $out_dir/tarid_exon2.alt.seq > $out_dir/tarid_exon2_result.alt.txt
$rm_dir/rmout.py < $out_dir/tarid_exon2_result.alt.txt
$rm_dir/rmout.py --out=$fig_dir/tarid_exon2.alt.pdf --cands=5 fold < $out_dir/tarid_exon2_result.alt.txt


# get TARID reverse kcomplement sequence:
echo 'AGTTGTGAGAAGGAACCCCGAGGACTTCTGCAGCAGAGGAACCGGCTCTCTGTTTCCCCACTCCCAGAGTCGTGGCCGTAGAGGAGGGTGAGCGAGCGCTGAGGAATTTGGTGGACACTAGGAATTTATCTGGGGAATAGAGGGGCGGATCCTTGCAGCCCAGGGAGGGGCCTGCCAGGCCCCAGCAAAGCGCTTAGACCCCTTTCACAACCGGAGGGAAACTCAATGCACAGACCCTGATTTGCAACTTGTAATGTAAATCAACTCAACTGCATCATGTACTTACCAGCCACCTTCTCCCAACTGTGGCCGGGGCCGAGAAGCAGCATGTTCTGTCCATCTGTAAAAGGCCTTCTTTCTCTCTTAAGACGTCACAACTGGTTGTTACTGAGAACTTTAGAAAAACGACTAGATCGTTTGGCTCTTTCTGGGCCTTCCACCTACACTCCATGCCCTCTCTCCTCCCAAGAAAACACGAATTAAAATACAGGGAGTATCAACCGCATCCTGCCAACACCTCATAAGACTCCAAGGCTTTACCGAGATCTTCCACCTTCCCCGGGATAAATGAAGACAAGTGCTTCCTAGAATCTTTATGGTGATGGAAAA' >  $out_dir/tarid_exon2.rc.seq


# Predict motifs: 
$rm_dir/rmdetect.py --both-strands $out_dir/tarid_exon2.rc.seq > $out_dir/tarid_exon2_result.rc.txt
$rm_dir/rmout.py < $out_dir/tarid_exon2_result.rc.txt
$rm_dir/rmout.py --out=$fig_dir/tarid_exon2.rc.pdf --cands=11 fold < $out_dir/tarid_exon2_result.rc.txt


# get TARID alternative sequence:
echo 'AGTTGTGAGAAGGAACCCCGAGGACTTCTGCAGCAGAGGAACCGGCTCTCTGTTTCCCCACTCCCAGAGTCGTGGCCGTAGAGGAGGGTGAGCGAGCGCTGAGGAATTTGGTGGACACTAGGAATTTATCTGGGGAATAGAGGGGCGGATCCTTGCAGCCCAGGGAGGGGCCTGCCAGGCCCCAGCAAAGCGCTTAGACCCCTTTCACAACCGGAGGGAAACTCAATGCACAGACCCTGATTTGCAACTTGTAATGTAAATCAACTCAACTGCATCATGTACTTACCAGCCACCTTCTCCCAACTGTGGCCGGGGCCGAGAAGCAGCATGTTCTGTCCATCTGTAAAAGGCCTTCTTTCTCTCTTAAGACGTCACAACTGGTTGTTACTGAGAACTTTAGAAAAACGACTAGATCGTTTGGCTCTTTCTGGGCCTTCCACCTACACTCCATGCCCTCTCTCCTCCCAAGAAAACACGAATTAAAATACAGGGAGTATCAACCGCATCCTGCCAGCACCTCATAAGACTCCAAGGCTTTACCGAGATCTTCCACCTTCCCCGGGATAAATGAAGACAAGTGCTTCCTAGAATCTTTATGGTGATGGAAAA' > $out_dir/tarid_exon2.rc.alt.seq

# Predict motifs: 
$rm_dir/rmdetect.py --both-strands $out_dir/tarid_exon2.rc.alt.seq > $out_dir/tarid_exon2_result.rc.alt.txt
$rm_dir/rmout.py < $out_dir/tarid_exon2_result.rc.alt.txt
$rm_dir/rmout.py --out=$fig_dir/tarid_exon2.rc.alt.pdf --cands=12 fold < $out_dir/tarid_exon2_result.rc.alt.txt
$rm_dir/rmout.py --out=$fig_dir/tarid_exon2.rc.alt.5.pdf --cands=5 fold < $out_dir/tarid_exon2_result.rc.alt.txt

