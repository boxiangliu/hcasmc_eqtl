#!/bin/bash 
# boxiang liu
# durga
# run count_num_sig_association.py across all matrix eqtl output
processed_data=../processed_data/160530/
scripts=./160530/
for num_geno_pc in 3 4 5; do
	echo $num_geno_pc
	for num_peer_factor in {1..15}; do
		echo $num_peer_factor
		input=$processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc$num_geno_pc.peer$num_peer_factor.2.cis.txt
		
		# count the number of significant associations with FDR:
		output1=$(python $scripts/count_num_sig_association.py \
					$input \
					0.1,0.05,0.01,0.001 \
					fdr)
		echo -e "$num_geno_pc\t$num_peer_factor\t$output1" >> $processed_data/find_optimum_num_PEER_factors_matrixeqtl/num_eqtls_vs_cov.fdr.txt
		
		# count the number of significant associations with p-value
		output2=$(python $scripts/count_num_sig_association.py \
			$input \
			1e-3,1e-5,1e-8 \
			pvalue)
		echo -e "$num_geno_pc\t$num_peer_factor\t$output2" >> $processed_data/find_optimum_num_PEER_factors_matrixeqtl/num_eqtls_vs_cov.pvalue.txt
	done 
done 
