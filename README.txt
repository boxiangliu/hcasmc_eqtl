# TOC 
- 160520 RNA sample outlier 
- 160526 WGS sample contamination
- 160527 PEER correction 
- 160530 eQTL mapping
- 160603 RNA expression clustering 
- 160604 phasing 
- 160614 differential expression
- 160615 GWAS overlap 
- 160618 colocalization
- 160627 sQTL



# 160527

# 160530
## 1 matrix_eQTL
1. prepare covariates
    combine_covariates.R
2. prepare genotype data
    bcftools 
3. prepare expression data
    normalize_rpkm.R
4. prepare genotype location
    gen_snps_loc.R
5. prepare gene location
    gen_gene_loc.R
6. find optimum number of PC and PEER factors: 
    find_optimal_num_PEER_factors_matrix_eQTL.sh 
    run_count_num_sig_association.sh 
    plot_num_eqtl_vs_cov.R
6. run matrix eQTL
    run_matrix_eQTL.R

Conclusion:
3 PCs and 8 PEER factors give the most eQTLs
The output file is here:
    
    $processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer8.2.cis.txt


## 2 fastQTL
1. prepare genotype data
    bcftools
2. prepare expression data
    `gen_bed.R`
3. prepare covariates
    `combine_covariates.R`
4. run fastqtl nominal pass
    `run_fastqtl.nominal.sh`
5. compare fastqtl 10000 with 100000 permutations
    `run_fastqtl.sh`
    `qqplot_fastqtl_pvalue.R`
    `compare_10000_100000_permutations.R`
    `fastqtl_pvalue_corrections.R`
6. find optimal number of PC and PEER factors
    `find_optimal_num_PEER_factors.sh`
    `plot_num_egene_vs_cov.R`
7. compare number of eGenes with GTEx
    `plot_num_egenes_vs_sample_size.R`


## 3 trans-eQTL (using matrix eQTL)
1. run trans eQTL using 3 PCs and 8 PEER factors 
    `run_matrix_eQTL.R`


# 160603 (expressiong clustering)
## 0 reference

[Figure](http://www.broadinstitute.org/collaboration/gtex_analysis/index.php/Temp_-_Transcriptome#4._Tissue_clustering.2FTissue_identity_from_Expression_Profiles) 
[Method](http://www.broadinstitute.org/collaboration/gtex_analysis/index.php/Transcriptome_analysis_(CRG)#Clustering_of_gene_expression.)
    
## 1 gene expression clustering
0. run RNAseQC
    `run_RNAseQC.sh`

1. copy rpkm into correct place
    `prepare_gtex_rpkm.sh`

2. combine rpkm across all tissues, filter for genes >0.1 rpkm in >10 individuals
    `combine_and_filter_rpkm.R`

3. hierarchical clustering
    `hclust.R`
conclusion: HCASMC clustered with fibroblasts...

4. MDS
    `mds.R`
5. 


# 160604_phasing:
## compare the target VCF with reference VCF:
    convert_chrom_names.sh
    subset_vcf_to_first_4_col.sh
    compare_vcf_with_reference.scratch.R

conclusion: 
Many difference exists between the target and the reference VCF, especially in the ALT field. For example: 
22_16858229     C     C     T     A
These variants cannot be phased property, and conform-gt needs to be run. 

## run conform-gt
    make_list_of_non_caucasian_samples.hcasmc.R
    make_list_of_non_EUR_samples.1kg.sh
    run_conform_gt.sh

Result: 
In total, the conformed VCF has 12523334 lines, compared to 15463398 lines in the original VCF file. The difference is 2940064 lines, or 19%. 

## phasing: 
    run_beagle_phasing.sh # with reference
    run_beagle_phasing_without_reference.sh # without reference


# 160614 (differential expression)
1. getting raw counts 
2. run DESeq

## reference: 
[fibroblast marker genes](https://www.researchgate.net/post/Which_genes_are_highly_expressed_in_fibroblast_cells)
[SMC marker genes](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943591/)


# 160615

# 160618 (colocalization)
# TODO: 
## reference 
[CAD prevalance](https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6040a1.htm)

#### sQTL 
# 160627 (sQTL)

# figures

normalize.R
    intron_missingness_rate.pdf


# 160628 (beagle)
beagle_QC.R
    R2_histogram.pdf          
    dosageR2_by_AFbin.pdf
    dosageR2_by_AFbin_size_002.pdf
    dosageR2_by_chrom.pdf
    # TODO: 
2. fine mapping
4. sQTL mapping 


# 160629 (sqtl)
- map sQTLs with fastqtl nominal pass 

- make qqplot and histogram 

- diagnose enrichment of high pvalues
    diagnostic_p_value_dist.R
        gt6_count_vs_adj_min_pval.pdf
        gt6_count_vs_mean_pval.pdf
        gt6_count_vs_min_pval.pdf
        mean_count_vs_adj_min_pval.pdf
        mean_count_vs_mean_pval.pdf
        mean_count_vs_min_pval.pdf
        zero_count_vs_adj_mean_pval.pdf
        zero_count_vs_mean_pval.pdf
        zero_count_vs_min_pval.pdf

- plot number of sig. sqtl vs distance
    + plot_sqtl_vs_distance.R
        * num_sig_sqtl_vs_dist.pdf
        * num_sig_sqtl_within_intron_vs_dist.pdf

- visualize bam files


# 160705 
# diagnose the p-value deflation problem.
- map sQTL on permuted data
    + sqtl.perm.histogram.pdf
    + sqtl.perm.qqplot.pdf 
        * observation: deflation still exists

- map sQTL on permuted data without PEER factors
    + sqtl.perm.histogram.noPEER.pdf
    + sqtl.perm.qqplot.noPEER.pdf
        * observation: deflation disappeared

- map sQTL on permuted data with PEER factors using top 10000 introns
    + sqtl.perm.histogram.top10000.pdf
    + sqtl.perm.qqplot.top10000.pdf
        * observation: deflation became small but still exists

- map sQTL on permuted data with 3,6,9,12,15 PEER factors (top 10000 introns)
    + sqtl.perm.histogram.top10000.peer{3,6,9,12,15}.pdf
    + sqtl.perm.qqplot.top10000.peer{3,6,9,12,15}.pdf
        * observation: deflation lessens as the number of factors decreases 

# 160708
# WASP correction


#### splicing QTLs (sQTLs)
#### 160629
(TBD) Leafcutter sQTL analysis
One of the top sIntrons for the SMAD3 loci is de novo. We visualized the loci by displaying RNAseq coverage summed across all samples on UCSC browser (figures/160629/smad3_loci_denovo_intron_RNAseq_coverage.png). The de novo intron is supported by some, but not much, RNAseq coverage. One sample, 1020301, does not show coverage in support of the de novo intron. One could suspect that this sample possesses a variant that favors against this intron. 



#### HCASMC specific genes 
#### 160715
subsample.list.by.tissue.bl.R subsampled 10 individuals for each tissue and wrote the individuals to one txt file per tissue. 

subsample.R generated subsampled gct files using lists generated by subsample.list.by.tissue.bl.R.

combine_read_counts.sh combined all gct files of HCASMC into one gct file. 

find_housekeeping_genes.R found ~9,435 not differentially expressed across tissues with FDR 99%. The p-value distribution has a peak at 0 and another at 1 (see figures/160715/DE_p_value_distribution.pdf). In addition, p-value is inversely correlated with mean read depth. However, there are some genes with relatively high read depth and insignificant p-values (see figures/160715/mean_read_depth_vs_pvalue.pdf). The plot figures/160715/read_counts_across_tissue_for_DE_genes.pdf shows a couple examples of differentailly expressed genes. MDGA2 (ENSG00000139915.14) is a brain specific gene and FOXP1 (ENSG00000114861.14) is a transcription factor whose expression varies across tissue. 

Eisenberg and Levanon reported 3804 housekeeping genes across 16 human tissues. They also provided a narrow list of 11 highly uniform genes.  The plot figures/160715/read_counts_across_tissue_for_housing_keeping_genes_by_Levanon.pdf shows that some of reported housekeeping genes may not be uniformly expressed. 


To select bona fide invariant housekeeping gene, we employ a equivalence testing procedure in an ANOVA framework. However, since ANOVA assumes samples are normally distribution with equal variances, we needed to transform the data. Among log(x), sqrt(x), x^(2/3) and variance stabilizing transformations, the first and last performed the best in making the distribution symmetric (figures/160715/compare_transformations.pdf). In terms of decompling variance from the mean, the last performed the best(figures/160715/mean_vs_variance.pdf). The figure figures/160715/expression_by_tissue_raw_count_vs_vsn.pdf shows pre and post variance stabilizating transformation. It is difficult to judge whether the variances are more constant, but expression outliers has smaller effects after VSN. In addition, the mean variance dependency has been reduced. The equivalence test is implemented in equivalence.R. The testing procedure follows chapter 7 of the book testing_statistical_hypotheses_of_equivalence_and_noninferiority. However, I did not carry on with this approach since I do not fully understand it. 


An alternative approach is to set a minimum threshold on the expression value and a maximum threshold on the expression variance (as described in Levanon et al). These following filters: 
(i) expression observed in all tissues; 
(ii) low variance over tissues: standard-deviation [log2(RPKM)]<1;
(iii) no log2(rpkm) differed from the averaged by 1 (twofold) or more
identified 3391 housekeeping genes. Two examples of housekeeping and non-housekeeping genes are in figures/160715/{read_counts_for_hk_gene.pdf,read_counts_for_non_hk_gene.pdf}


RUVSeq was used to correct for the unwanted hidden variables. Removing 20 hidden factors were sufficient, as the relative log expression for each sample is distributed with ~0 mean and similar variances (/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/rle_and_pca_after_ruv.*.pdf)


HCASMC specific gene were found using a rank based statistic described as follows. Given a list L1 of N numbers 1:N and a second list L2 of n numbers, the test calculates the probability that the maximum member of an arbitrary list of lengh n is less than the maximum memeber of the list L2. In particular, p=choose(max(L),n)/choose(N,n). P-values using permuted dataset shows that the test is well-calibrated (/srv/persistent/bliu2/HCASMC_eQTL/figures/160715/qqplot_pvalues.pdf)After BH adjustment, 123 genes are statistically significant at FDR 5%. 



# TO READ:
1. WGCNA paper
2. Sazonova (2015) Plos Genetics
3. iPSC eQTL paper
4. leaf cutter

