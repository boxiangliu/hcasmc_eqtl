#!/usr/bin/env Rscript
# durga
# bosh liu
# find putative target genes of 7 GWAS loci

# load GWAS variants:
# TCF21 rs12190287
# * SMAD3 rs17293632 15:67442596 top GWAS variant, disrupts AP1 motif

# load eQTL dataset: 
cis=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer8.2.cis.txt',head=T)
trans=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/trans.txt')

# find genes influenced by GWAS loci:
# CDKN2BAS rs1537373 chr9:22103341
cis[SNP=='chr9_22103341_T_G',.(SNP,gene,`p-value`)]
 #                  SNP               gene    p-value
 # 1: chr9_22103341_T_G  ENSG00000197919.3 0.09053498
# ENSG00000197919.3 IFNA1
# The protein encoded by this gene is produced by macrophages and has antiviral activity. This gene is intronless and the encoded protein is secreted. [provided by RefSeq, Sep 2011]


# SMAD3 rs17293632 chr15_67442596_C_T
cis[SNP=='chr15_67442596_C_T', .(SNP,gene,`p-value`)]
#                SNP               gene     p-value
# 1: chr15_67442596_C_T  ENSG00000103591.8 0.002256956
# 2: chr15_67442596_C_T ENSG00000174444.10 0.006744905
# 3: chr15_67442596_C_T  ENSG00000188501.7 0.028284410
# 4: chr15_67442596_C_T ENSG00000103599.15 0.063346500
# 5: chr15_67442596_C_T  ENSG00000166938.8 0.065635317
# 6: chr15_67442596_C_T ENSG00000166949.11 0.072012422
# ENSG00000103591.8 AAGAB
# The protein encoded by this gene interacts with the gamma-adaptin and alpha-adaptin subunits of complexes involved in clathrin-coated vesicle trafficking. Mutations in this gene are associated with type I punctate palmoplantar keratoderma. Alternatively spliced transcript variants have been found for this gene. [provided by RefSeq, Dec 2012]
# ENSG00000174444.10 RPL4
# Ribosomes, the organelles that catalyze protein synthesis, consist of a small 40S subunit and a large 60S subunit. Together these subunits are composed of 4 RNA species and approximately 80 structurally distinct proteins. This gene encodes a ribosomal protein that is a component of the 60S subunit. The protein belongs to the L4E family of ribosomal proteins. It is located in the cytoplasm. As is typical for genes encoding ribosomal proteins, there are multiple processed pseudogenes of this gene dispersed through the genome. [provided by RefSeq, Jul 2008]
# ENSG00000188501.7 LCTL
# This gene encodes a member of family 1 glycosidases. Glycosidases are enzymes that hydrolyze glycosidic bonds and are classified into families based on primary amino acid sequence. Most members of family 1 have two conserved glutamic acid residues, which are required for enzymatic activity. The mouse ortholog of this protein has been characterized and has a domain structure of an N-terminal signal peptide, glycosidase domain, transmembrane domain, and a short cytoplasmic tail. It lacks one of the conserved glutamic acid residues important for catalysis, and its function remains to be determined (PMID: 12084582). Alternative splicing results in multiple transcript variants. [provided by RefSeq, Jun 2013]
# ENSG00000103599.15 IQCH
# May play a regulatory role in spermatogenesis, IQ motif is a protein motif that binds to cadmodulin
# ENSG00000166938.8 DIS3L
# DIS3L (DIS3 Like Exosome 3'-5' Exoribonuclease) is a Protein Coding gene. Among its related pathways are Deadenylation-dependent mRNA decay. GO annotations related to this gene include RNA binding and 3-5-exoribonuclease activity. An important paralog of this gene is DIS3L2.
# ENSG00000166949.11 SMAD3 
# SMAD3 (SMAD Family Member 3) is a Protein Coding gene. Diseases associated with SMAD3 include smad3-related loeys-dietz syndrome and smad3-related thoracic aortic aneurysms and aortic dissections. Among its related pathways are Signaling by GPCR and Downstream signaling events of B Cell Receptor (BCR). GO annotations related to this gene include transcription factor activity, sequence-specific DNA binding and sequence-specific DNA binding. An important paralog of this gene is SMAD1.
trans[SNP=='chr15_67442596_C_T']


# PDGFD rs2019090 chr11:103668962
cis[SNP=='chr11_103668962_A_T',.(SNP,gene,`p-value`)]
# 1: chr11_103668962_A_T ENSG00000110347.7 0.05730896
# 2: chr11_103668962_A_T ENSG00000260966.1 0.09831423
# * ENSG00000110347.7 MMP12
# Proteins of the matrix metalloproteinase (MMP) family are involved in the breakdown of extracellular matrix in normal physiological processes, such as embryonic development, reproduction, and tissue remodeling, as well as in disease processes, such as arthritis and metastasis. Most MMP's are secreted as inactive proproteins which are activated when cleaved by extracellular proteinases. It is thought that the protein encoded by this gene is cleaved at both ends to yield the active enzyme, but this processing has not been fully described. The enzyme degrades soluble and insoluble elastin. It may play a role in aneurysm formation and studies in mice suggest a role in the development of emphysema. The gene is part of a cluster of MMP genes which localize to chromosome 11q22.3. [provided by RefSeq, Jul 2008]
# ENSG00000260966.1 ENSG00000260966
# ENSG00000260966 is an RNA Gene, and is affiliated with the lncRNA class.
trans[SNP=='chr11_103668962_A_T']

# IL6R rs7549250 chr1_154404336_C_T:
cis[SNP=="chr1_154404336_C_T",.(SNP,gene,`p-value`)]
 # 1: chr1_154404336_C_T  ENSG00000160712.8 0.009881121
 # 2: chr1_154404336_C_T ENSG00000169241.13 0.024570834
 # 3: chr1_154404336_C_T  ENSG00000143624.9 0.040916104
 # 4: chr1_154404336_C_T ENSG00000143621.12 0.067681228
# * ENSG00000160712.8 IL5R
# The interleukin-5 receptor; recognize and respond to cytokines; Through binding to the interleukin-5 receptor, interleukin 5 stimulates B cell growth and increases immunoglobulin secretion. It is also a key mediator in eosinophil activation.
# reference 1: Endogenous interleukin-1 alpha promotes a proliferative and proinflammatory phenotype in human vascular smooth muscle cells.
# referencd 2: Upregulation of fibronectin synthesis by interleukin-1 beta in coronary artery smooth muscle cells is associated with the development of the post-cardiac transplant arteriopathy in piglets.
# ENSG00000169241.13 SLC50A1 
# (Solute Carrier Family 50 (Sugar Efflux Transporter), Member 1) is a Protein Coding gene. Diseases associated with SLC50A1 include cerebral palsy, spastic quadriplegic, 1. GO annotations related to this gene include glucoside transmembrane transporter activity.
# ENSG00000143624.9 INTS3
# INTS3 is a subunit of the Integrator complex, which associates with the C-terminal domain of RNA polymerase II large subunit (POLR2A; MIM 180660) and mediates 3-prime end processing of small nuclear RNAs U1 (RNU1; MIM 180680) and U2 (RNU2; MIM 180690) (Baillat et al., 2005 [PubMed 16239144]). INTS3 is also a subunit of single-stranded DNA (ssDNA)-binding complexes involved in the maintenance of genome stability (Huang et al., 2009) [PubMed 19683501].[supplied by OMIM, Feb 2010]
# ENSG00000143621.12 LF2
# The protein encoded by this gene is a transcription factor required for T-cell expression of the interleukin 2 gene. It also binds RNA and is an essential component for encapsidation and protein priming of hepatitis B viral polymerase. The encoded 45 kDa protein (NF45, ILF2) forms a complex with the 90 kDa interleukin enhancer-binding factor 3 (NF90, ILF3), and this complex has been shown to affect the redistribution of nuclear mRNA to the cytoplasm, to repair DNA breaks by nonhomologous end joining, and to negatively regulate the microRNA processing pathway. Knockdown of NF45 or NF90 protein retards cell growth, possibly by inhibition of mRNA stabilization. Alternative splicing results in multiple transcript variants. Related pseudogenes have been found on chromosomes 3 and 14. [provided by RefSeq, Dec 2014]
trans[SNP=="chr1_154404336_C_T",.(SNP,gene,`p-value`)]


# BMP1 rs73551707
cis[SNP=="chr8_22046423_C_T",.(SNP,gene,`p-value`)]
# 1: chr8_22046423_C_T ENSG00000120910.10 0.06034759 
# 2: chr8_22046423_C_T  ENSG00000168495.8 0.08322178
# 3: chr8_22046423_C_T  ENSG00000104635.9 0.08946980 
# ENSG00000120910.10 PPP3CC
# Calcineurin is a calcium-dependent, calmodulin-stimulated protein phosphatase involved in the downstream regulation of dopaminergic signal transduction. Calcineurin is composed of a regulatory subunit and a catalytic subunit. The protein encoded by this gene represents one of the regulatory subunits that has been found for calcineurin. Three transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Sep 2011]
# ENSG00000168495.8 POLR3D
# POLR3D (Polymerase (RNA) III (DNA Directed) Polypeptide D, 44kDa) is a Protein Coding gene. Among its related pathways are Immune System and Epstein-Barr virus infection. GO annotations related to this gene include chromatin binding and RNA polymerase III activity.
# ENSG00000104635.9 SLC39A14
# Zinc is an essential cofactor for hundreds of enzymes. It is involved in protein, nucleic acid, carbohydrate, and lipid metabolism, as well as in the control of gene transcription, growth, development, and differentiation. SLC39A14 belongs to a subfamily of proteins that show structural characteristics of zinc transporters (Taylor and Nicholson, 2003 [PubMed 12659941]).[supplied by OMIM, Mar 2008]
trans[SNP=="chr8_22046423_C_T",.(SNP,gene,`p-value`)]


# CCDC97 rs2241718 chr19:41829606
cis[SNP=="chr19_41829606_G_A",.(SNP,gene,`p-value`)]
 # 1: chr19_41829606_G_A ENSG00000188493.10 0.01920192
 # 2: chr19_41829606_G_A  ENSG00000167601.7 0.02752588
 # 3: chr19_41829606_G_A ENSG00000142046.10 0.03394394
 # 4: chr19_41829606_G_A ENSG00000105223.14 0.05133411
 # 5: chr19_41829606_G_A  ENSG00000268041.1 0.09700477
# ENSG00000188493.10 C19orf54 
# Chromosome 19 Open Reading Frame 54
# * ENSG00000167601.7 AXL 
# The protein encoded by this gene is a member of the Tyro3-Axl-Mer (TAM) receptor tyrosine kinase subfamily. The encoded protein possesses an extracellular domain which is composed of two immunoglobulin-like motifs at the N-terminal, followed by two fibronectin type-III motifs. It transduces signals from the extracellular matrix into the cytoplasm by binding to the vitamin K-dependent protein growth arrest-specific 6 (Gas6). This gene may be involved in several cellular functions including growth, migration, aggregation and anti-inflammation in multiple cell types. Alternative splicing results in multiple transcript variants of this gene. [provided by RefSeq, Jul 2013]
# ENSG00000142046.10 TMEM91 
# TMEM91 (Transmembrane Protein 91) is a Protein Coding gene. An important paralog of this gene is SYNDIG1, which encodes a protein that belongs to the interferon-induced transmembrane family of proteins. A similar protein in rat is thought to regulate the development of excitatory synapses. [provided by RefSeq, Jul 2013]
# ENSG00000105223.14 PLD3
# This gene encodes a member of the phospholipase D (PLD) family of enzymes that catalyze the hydrolysis of membrane phospholipids.
# ENSG00000268041.1 LOC390937
# LOC390937 (Ets2 Repressor Factor-Like) is a Protein Coding gene.