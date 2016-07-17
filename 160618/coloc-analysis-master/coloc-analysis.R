library(data.table)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(stringr)



colocs <- fread("~/MontgomeryRotation/all_colocs.txt")
colocs <- fread("~/all_colocs.txt")

expr1 <- fread("~/MontgomeryRotation/all_expression.txt")
expr2 <- fread("~/MontgomeryRotation/all_expression.txt")

colnames(expr1) <- c('gene', 'expr_tissue1', 'tissue1')
colnames(expr2) <- c('gene', 'expr_tissue2', 'tissue2')

colocs_expr <- left_join(colocs, expr1)
colocs_expr <- left_join(colocs_expr, expr2)

n1 <- fread("~/MontgomeryRotation/n_samples.txt")
n2 <- fread("~/MontgomeryRotation/n_samples.txt")

colnames(n1) <- c("tissue1", "n_tissue1")
colnames(n2) <- c("tissue2", "n_tissue2")

colocs_expr <- left_join(colocs_expr, n1)
colocs_expr <- left_join(colocs_expr, n2)

colocs2 <- gather(colocs_expr, hypothesis, pp, pp.H0:pp.H4)
colocs2 <- gather(colocs, hypothesis, pp, pp.H0:pp.H4)

colocs2 <- colocs2 %>% group_by(tissue1, tissue2, gene) %>% filter(pp == max(pp))

colocs2$h4 <- ifelse(colocs2$hypothesis == "pp.H4", 1, 0)


c2 <- colocs2 %>% group_by(gene) %>% summarize(percshared = mean(h4))
c3 <- colocs2 %>% group_by(tissue1, tissue2, gene) %>% summarize(percshared = mean(h4))


###### Pull out most shared genes to look at overlap with CTCF binding sites

mostshared <- inner_join(c2[c2$percshared >= .9, ], colocs2)

mostshared$chrom <- sapply(mostshared$snp, function (x) paste0("chr", str_split(x, "_")[[1]][1]))
mostshared$start <- sapply(mostshared$snp, function(x) str_split(x, "_")[[1]][2])
mostshared$end <- as.numeric(mostshared$start) + 1
mostshared$name <- mostshared$gene
mostshared$score <- 1000

bed <- mostshared %>% select(chrom, start, end, name, score) %>%
  write.table(., "~/mostshared.bed", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)






genes <- read.table("~/UterusTop100.txt")

c2 <- filter(c2, gene %in% genes$gene)

hist(c2$percshared, breaks = 50, xlab = "Percentage H4 == max for all pairwise comparisons",
     ylab = "Count", main = "H4 Proportion Among All Tissue Combinations\nfor Top 100 Genes in Uterus Sample")



n <- read.table("~/MontgomeryRotation/n_samples.txt", header = FALSE)
colnames(n) <- c("tissue", "n")

top100 <- read.table("~/MontgomeryRotation/Top100_nohead.txt", header = TRUE)

top100 <- left_join(top100, c2)
top100 <- left_join(top100, n)

test <- top100 %>% group_by(tissue) %>% summarise(n = mean(n)) %>% arrange(n) %>% select(tissue)



top100$tissue <- factor(top100$tissue, levels = test$tissue)

ggplot(top100, aes(x = percshared)) +
  geom_bar() + facet_wrap(~ tissue, nrow = 7) +
  geom_text(aes(x = .10, y = 25, label = n)) +
  xlab("Percentage Shared Among All Tissues (Proportion H4+)") + ylab("Count")


# Generate the comparison heatmap

top100 <- read.table("~/MontgomeryRotation/Top100_nohead.txt", header = TRUE)

top100_genes <- left_join(top100, c3)
top100_genes <- top100_genes %>% group_by(tissue1, tissue2) %>% summarise(percshared = mean(percshared)) %>% filter(!is.na(tissue1), !is.na(tissue2))

c32 <- c3 %>% group_by(tissue1, tissue2) %>% summarise(percshared = mean(percshared)) %>% filter(!is.na(tissue1), !is.na(tissue2))
ggplot(c32, aes(tissue2, tissue1)) +
  geom_tile(aes(fill = percshared)) +
  scale_fill_gradient(low = "yellow", high = "blue") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0)) +
  xlab("Tissue 2") +
  ylab("Tissue 1") +
  ggtitle("Percent Colocalization for Top 100 Genes in Each Tissue 1")


ggplot(top100_genes, aes(tissue2, tissue1)) +
  geom_tile(aes(fill = percshared)) +
  scale_fill_gradient(low = "yellow", high = "blue") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0)) +
  xlab("Tissue 2") +
  ylab("Tissue 1") +
  ggtitle("Percent Colocalization for Top 100 Genes in Each Tissue 1")




#Look for dependence on Tajima's D

TD <- read.table("~/out.Tajima.D", header = TRUE)

TD <- TD %>% rename(pos = BIN_START, chrom = CHROM) %>% mutate(chrom = as.character(chrom))

# TD <- TD %>% mutate(BIN_END = BIN_START + 99999)

withpos <- colocs2 %>% separate(snp, into = c("chrom", "pos", "ref", "var", "garbage"))

percsharedwithpos <- withpos %>% group_by(tissue1, tissue2, chrom, pos, gene) %>% summarize(percshared = mean(h4))

round <- function(x) {
  x - (x %% 100000)
}

percsharedwithpos$pos <- round(as.numeric(percsharedwithpos$pos))

test <- inner_join(percsharedwithpos, TD, copy = TRUE) %>% select(chrom, pos, tissue1, tissue2, gene, TajimaD, percshared)

test2 <- test %>% group_by(TajimaD) %>% summarize(ps = mean(percshared))

plot(test2$TajimaD, test2$ps, xlab = "Tajima's D", ylab = "Pecent Shared (H4)", main = "Percent Shared vs. 100kb-binned Tajima's D")

summary(lm(ps ~ TajimaD, data = test2))
hist(test2$TajimaD)

# test <- sqldf("select * from percsharedwithpos f1 left join TD f2 on (f1.chrom == f2.CHROM and f1.pos >= f2.BIN_START and f1.pos <= f2.BIN_END)")


dlist <- list()
for (i in 1:nrow(TD)) {
  print(i)
  pos <- TD$BIN_START[i]:TD$BIN_END[i]
  chrom <- rep(TD$chrom[i], times = length(pos))
  D <- rep(TD$TajimaD[i], times = length(pos))
  out <- as.data.frame(chrom, pos, D)
  dlist <- append(dlist, out)
}

d <- rbind_all(dlist)




ngenes <- colocs2 %>% group_by(tissue1, tissue2) %>% summarize(ngenes = n(), propH4 = mean(h4))

plot(ngenes$ngenes, ngenes$propH4,
     xlab = "Number of genes compared",
     ylab = "Proportion for which H4 is greatest",
     main = "Proportion of Shared Causal Variants by Number of Genes Compared")

ngenesmodel <- lm((propH4) ~ (ngenes), data = ngenes)

summary(ngenesmodel)

lmmodel <- lm(h4 ~ (expr_tissue1) *
              (expr_tissue2) *
              (n_tissue1) *
              (n_tissue2),
            data = colocs2)

lmmodel <- lm(h4 ~ (expr_tissue1) +
                (expr_tissue2) +
                (n_tissue1) +
                (n_tissue2) +
                tissue1 +
                tissue2,
              data = colocs2)

logmodel <- glm(h4 ~ (expr_tissue1) *
                 (expr_tissue2) *
                 (n_tissue1) *
                 (n_tissue2),
               data = colocs2, family = binomial)


summary(lmmodel)
summary(logmodel)

















glm.out = glm(cbind(Menarche, Total-Menarche) ~ Age, family=binomial(logit), data=menarche)


combos <- colocs2 %>% group_by(tissue1, tissue2, gene) %>%
  filter(pp == max(pp)) %>%
  ungroup() %>%
  group_by(tissue1, tissue2) %>%
  summarise(fraction = sum(hypothesis == "pp.H4")/n())

ggplot(combos, aes(tissue2, tissue1)) +
  geom_tile(aes(fill = fraction)) +
  scale_fill_gradient(low = "yellow", high = "blue") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0))









combomat <- spread(combos, tissue2, fraction) %>% as.data.frame(.)
rownames(combomat) <- combomat[,1]
combomat <- combomat %>% select(-tissue1) %>% as.matrix(.)

combomat <- combomat[,-1] %>% as.matrix(.)

combomat[diag(combomat)] <- 1

heatmap(combomat)

d <- dist(combos)
d2 <- as.matrix(d)

d3 <- d2[order(d2[,1]),]



plot <- ggplot(colocs2, aes(x = pp, fill = hypothesis)) +
  geom_bar(position = "dodge") + facet_grid(tissue1~tissue2)

colocs2 <- colocs2 %>% group_by(gene) %>% filter(pp == max(pp)) %>% filter(hypothesis == "pp.H4")











pca1 = prcomp(combomat)



plot <- ggplot(colocs2, aes(x = pp, fill = hypothesis)) +
  geom_bar(position = "dodge") + facet_grid(tissue1~tissue2)



colocs2 <- colocs2 %>% group_by(gene) %>% filter(pp == max(pp)) %>% filter(hypothesis == "pp.H4")
colocs2 <- colocs2 %>% group_by(gene) %>% mutate(all = n()) %>%
  filter(pp == max(pp)) %>% filter(hypothesis == "pp.H4") %>% mutate(subset = n()) %>%
  summarize(mean(all), mean(subset))



c2 <- colocs %>% select(tissue1, tissue2, nsnps, snp.pp.H4, pp.H0, pp.H1, pp.H2, pp.H3, pp.H4) %>%
  group_by(tissue1, tissue2) %>% summarise_each(funs(mean)) %>% data.frame(.)
rownames(c2) <- paste(c2$tissue1, c2$tissue2, sep = "_")
c2 <- select(c2, -tissue1, -tissue2)

pca1 = prcomp(c2)
pca <- princomp(c2)


colocs2$chrom <- sapply(str_split(colocs2$snp, "_"), function(x) x[1])
colocs2$pos <- sapply(str_split(colocs2$snp, "_"), function(x) x[2])

g <- ggplot(colocs2, aes(x = pos, y = pp)) + geom_point() + facet_grid(.~chrom)
ggsave(filename = "~/testplot.pdf", plot = g, width = 500, height = 500, limitsize = FALSE)

colocs3 <- colocs2 %>% group_by(chrom, pos) %>% summarise(count = n())
colocs3 <- colocs2 %>% group_by(chrom) %>% summarise(count = n())

ggplot(colocs3, aes(x = pos)) + geom_bar() + facet_grid(.~chrom)

combos <- colocs2 %>% group_by(tissue1, tissue2) %>%
  filter(pp == max(pp)) %>%
  summarise(fraction = sum(hypothesis == "pp.H4"))

combomat <- spread(combos, tissue2, fraction) %>% as.data.frame(.)
rownames(combomat) <- combomat[,1]
combomat <- combomat %>% select(-tissue1) %>% as.matrix(.)
combomat <- rbind(combomat, NA)
combomat <- rbind(combomat, NA)

g <- graph.adjacency(combomat, weighted = TRUE, mode = "undirected")
g <- simplify(g)
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
set.seed(3952)
layout1 <- layout.fruchterman.reingold(g)
E(g)$weight
plot(g, layout=layout1, edge.width=scale(E(g)$weight))

tkplot(g)

