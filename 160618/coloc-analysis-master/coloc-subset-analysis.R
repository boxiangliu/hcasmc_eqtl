library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(multidplyr)

# Read in the files

colocs <- fread("/srv/scratch/tshimko/results/subset/all_colocs.txt")
colocs <- fread("~/all_colocs.txt")

# Gather all of the hypothesis columns into one key, value pair and keep only the
# hypothesis with the highest posterior probability

colocs2 <- gather(colocs, hypothesis, pp, pp.H0:pp.H4)
colocs2 <- colocs2 %>% group_by(tissue1, tissue2, gene) %>% filter(pp == max(pp))

# Multidplyerified

clust <- create_cluster(cores = 15)
clust <- create_cluster(cores = 3)

colocs2 <- gather(colocs, hypothesis, pp, pp.H0:pp.H4)
system.time({colocs3 <- partition(colocs2, cluster = clust)})
system.time({colocs4 <- filter(colocs3, pp == max(pp))})
system.time({colocs5 <- collect(colocs4)})

# Set a boolean stating whether h4 is the most likely for each tissue combiantion, gene pairing

colocs2$h4 <- ifelse(colocs2$hypothesis == "pp.H4", 1, 0)

# Reduce 

c2 <- colocs2 %>% group_by(gene) %>% summarize(percshared = mean(h4))
c3 <- colocs2 %>% group_by(tissue1, tissue2, gene) %>% summarize(percshared = mean(h4))

c32 <- c3 %>% group_by(tissue1, tissue2) %>% summarise(percshared = mean(percshared)) %>% filter(!is.na(tissue1), !is.na(tissue2))
ggplot(c32, aes(tissue2, tissue1)) +
  geom_tile(aes(fill = percshared)) +
  scale_fill_gradient(low = "yellow", high = "blue") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0)) +
  xlab("Tissue 2") +
  ylab("Tissue 1")



### OVERLAPPING RANGES ###

coloc_snps <- colocs2 %>% filter(h4 == 1) %>%
  separate(snp, into = c("chrom", "pos", "ref", "alt", "junk"), sep = "_") %>%
  select(-junk)

coloc_snps <- colocs2 %>%
  separate(snp, into = c("chrom", "pos", "ref", "alt", "junk"), sep = "_") %>%
  select(-junk)

features <- fread("/srv/scratch/tshimko/roadmap/25state/all_features.bed")
features <- fread("~/all_features.bed")

features <- GRanges(seqnames = features$V1,
                     ranges = IRanges(start = features$V2,
                     end = features$V3),
                     type = features$V4)

shared_snps <- GRanges(seqnames = paste0("chr", coloc_snps$chrom),
                       ranges = IRanges(start = as.numeric(coloc_snps$pos),
                                        width = sapply(coloc_snps$alt, str_length)),
                       gene = coloc_snps$gene,
                       tissue1 = coloc_snps$tissue1,
                       tissue2 = coloc_snps$tissue2,
                       ref = coloc_snps$ref,
                       alt = coloc_snps$alt,
                       snp.pp.H4 = coloc_snps$snp.pp.H4,
                       nsnps = coloc_snps$nsnps,
                       pp = coloc_snps$pp)

shared_snps_max <- GRanges(seqnames = paste0("chr", all_expanded$chrom.y),
                       ranges = IRanges(start = as.numeric(all_expanded$pos.y),
                                        width = sapply(all_expanded$alt.y, str_length)),
                       gene = all_expanded$gene,
                       tissue1 = all_expanded$tissue1,
                       tissue2 = all_expanded$tissue2,
                       ref = all_expanded$ref.y,
                       alt = all_expanded$alt.y,
                       pp = all_expanded$pp)

shared_snps_coloc <- GRanges(seqnames = paste0("chr", all_expanded$chrom.x),
                       ranges = IRanges(start = as.numeric(all_expanded$pos.x),
                                        width = sapply(all_expanded$alt.x, str_length)),
                       gene = all_expanded$gene,
                       tissue1 = all_expanded$tissue1,
                       tissue2 = all_expanded$tissue2,
                       ref = all_expanded$ref.x,
                       alt = all_expanded$alt.x,
                       snp.pp.H4 = all_expanded$snp.pp.H4,
                       nsnps = all_expanded$nsnps,
                       pp = all_expanded$pp)

snps_in_features = subsetByOverlaps(shared_snps, features)
features_hit_by_snps = subsetByOverlaps(features, shared_snps)

coloc_snps_in_features = subsetByOverlaps(shared_snps_coloc, features)
hits <- findOverlaps(shared_snps_coloc, features) %>% as.data.frame(.) %>% group_by(queryHits) %>% summarize(subjectHits = first(subjectHits))
idx <- hits$subjectHits
values <- data.frame(type = as.character(mcols(features)$type[idx]))
mcols(coloc_snps_in_features) <- c(mcols(coloc_snps_in_features), values)

max_snps_in_features = subsetByOverlaps(shared_snps_max, features)
hits <- findOverlaps(shared_snps_coloc, features)
idx <- unique(subjectHits(hits))
values <- data.frame(type = mcols(features)$type[idx])
mcols(max_snps_in_features) <- c(mcols(max_snps_in_features), values)

# Compare to eQTL snps

eqtl <- fread("/srv/scratch/tshimko/analysis/subset_eGenes/all_eQTL.txt") %>%
  as.data.frame(.) %>% dplyr::rename(tissue1 = tissue, snp = snps)
eqtl <- fread("~/all_eQTL.txt") %>%
  as.data.frame(.) %>% dplyr::rename(tissue1 = tissue, snp = snps)

all_shared <- colocs2 %>% filter(h4 == 1)

all_shared_with_snps <- left_join(all_shared, eqtl, by = c("tissue1", "gene"))

all_expanded <- all_shared_with_snps %>%
  separate(snp.x, into = c("chrom.x", "pos.x", "ref.x", "alt.x", "junk.x"), sep = "_") %>%
  separate(snp.y, into = c("chrom.y", "pos.y", "ref.y", "alt.y", "junk.y"), sep = "_") %>%
  select(-junk.x, -junk.y)

all_expanded <- all_expanded %>% filter(!is.na(chrom.x), !is.na(chrom.y))

features <- fread("/srv/scratch/tshimko/roadmap/25state/coloc_features.bed")
hits <- findOverlaps(features, features)
features <- GRanges(seqnames = features$V1,
                    ranges = IRanges(start = features$V2,
                                     end = features$V3),
                    type = features$V4)
values = data.frame(type1 = features$type[queryHits(hits)], type2 = features$type[subjectHits(hits)])

library(gplots)
load("~/calls.Rda")
heatmap(calls)


nhits <- sapply(1:length(shared_snps_coloc), function(x) {
  print(x)
  coloc_snp <- shared_snps_coloc[x,]
  length(findOverlaps(coloc_snp, features))
  })





########### Plotting functions for making the gif plots

# Make a GIF of increasing RPKM cutoff stringency

plot_rpkm_cutoff <- function(cutoff) {
  #   sharing <- colocs2 %>%
  #     filter(median_rpkm_tissue1 > cutoff, median_rpkm_tissue2 > cutoff) %>% 
  #     group_by(tissue1, tissue2) %>%
  #     summarise(percshared = mean(h4))
  
  sharing <- as.data.table(colocs2)
  sharing2 <- sharing[median_rpkm_tissue1 > cutoff]
  sharing3 <- sharing2[median_rpkm_tissue2 > cutoff]
  sharing4 <- sharing3[, mean(h4), by = "tissue1,tissue2"] %>% rename_(percshared = "V1") %>% as.data.frame(.)
  
  #     # Make the matrix to order the variables by heirarchecal clustering
  #     sharing_mat <- sharing4 %>% spread(tissue2, percshared) %>% data.frame(.)
  #     rownames(sharing_mat) <- sharing_mat$tissue1
  #     sharing_mat <- sharing_mat %>% select(-tissue1) %>% as.matrix(.)
  #     
  #     # Cluster on the matrix to set the ordering of the heatmap
  #     distances <- dist(sharing_mat)
  #     clustering <- hclust(distances)
  sharing4$tissue1 <- factor(sharing4$tissue1, levels = unique(sharing4$tissue1)[clustering$order])
  sharing4$tissue2 <- factor(sharing4$tissue2, levels = unique(sharing4$tissue1)[clustering$order])
  
  # Plot the heatmap
  ggplot(sharing4, aes(tissue2, tissue1)) +
    geom_tile(aes(fill = percshared)) +
    scale_fill_gradient(low = "yellow", high = "blue", limits = c(0, .018)) +
    theme(axis.text.x=element_text(angle=-45, hjust = 0)) +
    xlab("Tissue 2") +
    ylab("Tissue 1") + 
    ggtitle(paste0("Percentage of Causal Variant Sharing\nfor all Pairwise Tissue Comparisons\nMedian RPKM Threshold = ", cutoff))
  
  ggsave(filename = paste0("animation4-", sprintf("%02d", cutoff), ".png"), height = 10, width = 12)
}

for (i in seq(0,15)) {
  print(i)
  plot_rpkm_cutoff(i)
}

animate <- function() {
  lapply(seq(0,15), function(i) {
    plot_rpkm_cutoff(i)
  })
}









all_shared <- colocs2
all_shared_with_snps <- left_join(all_shared, eqtl, by = c("tissue1", "gene"))

all_expanded <- all_shared_with_snps %>%
  separate(snp.x, into = c("chrom.x", "pos.x", "ref.x", "alt.x", "junk.x"), sep = "_") %>%
  separate(snp.y, into = c("chrom.y", "pos.y", "ref.y", "alt.y", "junk.y"), sep = "_") %>%
  select(-junk.x, -junk.y) %>% group_by(gene) %>%
  summarize(percshared = mean(h4),
            chrom.x = names(sort(-table(.$chrom.x)))[1],
            pos.x = names(sort(-table(.$chrom.x)))[1],
            chrom.y = names(sort(-table(.$chrom.y)))[1],
            pos.y = names(sort(-table(.$chrom.y)))[1]
            ) %>% arrange(desc(percshared))
  
all_expanded <- all_expanded %>% filter(!is.na(chrom.x), !is.na(chrom.y))