##### Load and process GSE119262 expression dataset and annotations
# Expression 
gse119262 <- read.delim("~/Desktop/Biomarkers/Data/GSE119262/GSE119262_series_matrix_V2.txt", comment.char="!", header=TRUE, fill=TRUE)
# Header annotation 
gse119262.annot <- read.delim("~/Desktop/Biomarkers/Data/GSE119262/GSE119262_series_matrix_V2.txt", header=FALSE, fill=TRUE, nrows=9)
gse119262.annot <- t(gse119262.annot)
gse119262.annot <- data.frame(gse119262.annot[-1, ])
colnames(gse119262.annot) <- c("Patient", "ID_REF1", "HER2", "ER", "PR", "Response_AKT", "Ki67_Change", "Response_Ki67","ID_REF")
gse119262.annot$Sample <- sapply(strsplit(as.character(gse119262.annot$Patient), ": |\\.|_"), "[", 1)
gse119262.annot$PID <- sapply(strsplit(as.character(gse119262.annot$Patient), ": |\\.|_"), "[", 2)
gse119262.annot$Time <- sapply(strsplit(as.character(gse119262.annot$Patient), ": |\\.|_"), "[", 3)
gse119262.annot$Replicate <- sapply(strsplit(as.character(gse119262.annot$Patient), ": |\\.|_"), "[", 4)

# Array information 
gse119262.genes <- read.delim("~/Desktop/Biomarkers/Data/GSE119262/GPL6104-11576.txt", comment.char="#", header=TRUE)

# Aggregate genes in expression matrix
x <- match(gse119262$ID_REF, gse119262.genes$ID)
agg.gse119262 <- aggregate(as.matrix(gse119262[ 2:ncol(gse119262)]), by=list(gse119262.genes$Symbol[x]), FUN=mean)
gse119262.mat <- as.matrix(agg.gse119262[-c(1:4), 2:ncol(agg.gse119262)]) #not actual genes in first 4 rows
rownames(gse119262.mat) <- agg.gse119262$Group.1[-c(1:4)]

# Aggregate replicates 
temp <- aggregate(t(gse119262.mat), by=list(paste(gse119262.annot$PID, gse119262.annot$Time, sep="_")), FUN=mean)
agg.pre <- as.matrix(temp[grep("Pre", temp$Group.1), 2:ncol(temp)])
agg.post <- as.matrix(temp[grep("Post", temp$Group.1), 2:ncol(temp)])
rownames(agg.pre) <- temp$Group.1[grep("Pre", temp$Group.1)]
rownames(agg.post) <- temp$Group.1[grep("Post", temp$Group.1)]

# Retain paired sampes with both pre and post treatment data
cIN <- match(sapply(strsplit(rownames(agg.post), "_"), "[", 1), sapply(strsplit(rownames(agg.pre), "_"), "[", 1))
agg.pre <- agg.pre[cIN,]

# Response information  
resp.pre <- as.factor(gse119262.annot$Response_Ki67[match(rownames(agg.pre), paste(gse119262.annot$PID, gse119262.annot$Time, sep="_"))]) # 1=Non-responder 2=Responder
resp.post <- as.factor(gse119262.annot$Response_Ki67[match(rownames(agg.post), paste(gse119262.annot$PID, gse119262.annot$Time, sep="_"))]) # 1=Non-responder 2=Responder

agg.pre.H <- gsva(t(agg.pre), Hset, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())
agg.post.H <- gsva(t(agg.post), Hset, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())

agg.pre.C2 <- gsva(t(agg.pre), c2set, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())
agg.post.C2 <- gsva(t(agg.post), c2set, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())

# responders vs. non-responders 
pre.res.mat.H <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.H), ncol=4))
for (i in 1:nrow(agg.pre.H)) {
  res <- t.test(agg.pre.H[i,] ~ resp.pre)
  pre.res.mat.H[i,] <- c(diff(res$estimate), res$stderr, res$statistic, res$p.value)
}
colnames(pre.res.mat.H) <- c("Estimate", "SE", "t-value", "P")
rownames(pre.res.mat.H) <- rownames(agg.pre.H)
pre.res.mat.H$FDR <- p.adjust(pre.res.mat.H$P, method="fdr")
write.table(pre.res.mat.H, file="~/Desktop/Eleni Paper/pre.res.mat.H.tsv", sep="\t")

post.res.mat.H <- data.frame(matrix(data=NA, nrow=nrow(agg.post.H), ncol=4))
for (i in 1:nrow(agg.post.H)) {
  res <- t.test(agg.post.H[i,] ~ resp.post)
  post.res.mat.H[i,] <- c(diff(res$estimate), res$stderr, res$statistic, res$p.value)
}
colnames(post.res.mat.H) <- c("Estimate", "SE", "t-value", "P")
rownames(post.res.mat.H) <- rownames(agg.post.H)
post.res.mat.H$FDR <- p.adjust(post.res.mat.H$P, method="fdr")
write.table(post.res.mat.H, file="~/Desktop/Eleni Paper/post.res.mat.H.tsv", sep="\t")


pre.res.mat.C2 <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.C2), ncol=4))
for (i in 1:nrow(agg.pre.C2)) {
  res <- t.test(agg.pre.C2[i,] ~ resp.pre)
  pre.res.mat.C2[i,] <- c(diff(res$estimate), res$stderr, res$statistic, res$p.value)
}
colnames(pre.res.mat.C2) <- c("Estimate", "SE", "t-value", "P")
rownames(pre.res.mat.C2) <- rownames(agg.pre.C2)
pre.res.mat.C2$FDR <- p.adjust(pre.res.mat.C2$P, method="fdr")
write.table(pre.res.mat.C2, file="~/Desktop/Eleni Paper/pre.res.mat.C2.tsv", sep="\t")

post.res.mat.C2 <- data.frame(matrix(data=NA, nrow=nrow(agg.post.C2), ncol=4))
for (i in 1:nrow(agg.post.C2)) {
  res <- t.test(agg.post.C2[i,] ~ resp.post)
  post.res.mat.C2[i,] <- c(diff(res$estimate), res$stderr, res$statistic, res$p.value)
}
colnames(post.res.mat.C2) <- c("Estimate", "SE", "t-value", "P")
rownames(post.res.mat.C2) <- rownames(agg.post.C2)
post.res.mat.C2$FDR <- p.adjust(post.res.mat.C2$P, method="fdr")
write.table(post.res.mat.C2, file="~/Desktop/Eleni Paper/post.res.mat.C2.tsv", sep="\t")


# pre vs. post 
NR.pre <- grep("Non-Responder", resp.pre)
NR.post <- grep("Non-Responder", resp.post)
R.res.mat.H <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.H), ncol=4))
NR.res.mat.H <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.H), ncol=4))

for (i in 1:nrow(agg.pre.H)) {
  res1 <- t.test(agg.pre.H[i, -NR.pre], agg.post.H[i, -NR.post], paired=F)
  res2 <- t.test(agg.pre.H[i, NR.pre], agg.post.H[i, NR.post], paired=F) 
  
  R.res.mat.H[i,] <- c(diff(res1$estimate), res1$stderr, res1$statistic, res1$p.value)
  NR.res.mat.H[i,] <- c(diff(res2$estimate), res2$stderr, res2$statistic, res2$p.value)
}

colnames(R.res.mat.H) <- c("Estimate", "SE", "t-value", "P")
rownames(R.res.mat.H) <- rownames(agg.post.H)
R.res.mat.H$FDR <- p.adjust(R.res.mat.H$P, method="fdr")
write.table(R.res.mat.H, file="~/Desktop/Eleni Paper/R.res.mat.H.tsv", sep="\t")

colnames(NR.res.mat.H) <- c("Estimate", "SE", "t-value", "P")
rownames(NR.res.mat.H) <- rownames(agg.post.H)
NR.res.mat.H$FDR <- p.adjust(NR.res.mat.H$P, method="fdr")
write.table(NR.res.mat.H, file="~/Desktop/Eleni Paper/NR.res.mat.H.tsv", sep="\t")

R.res.mat.C2 <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.C2), ncol=4))
NR.res.mat.C2 <- data.frame(matrix(data=NA, nrow=nrow(agg.pre.C2), ncol=4))

for (i in 1:nrow(agg.pre.C2)) {
  res1 <- t.test(agg.pre.C2[i, -NR.pre], agg.post.C2[i, -NR.post], paired=F)
  res2 <- t.test(agg.pre.C2[i, NR.pre], agg.post.C2[i, NR.post], paired=F) 
  
  R.res.mat.C2[i,] <- c(diff(res1$estimate), res1$stderr, res1$statistic, res1$p.value)
  NR.res.mat.C2[i,] <- c(diff(res2$estimate), res2$stderr, res2$statistic, res2$p.value)
}

colnames(R.res.mat.C2) <- c("Estimate", "SE", "t-value", "P")
rownames(R.res.mat.C2) <- rownames(agg.post.C2)
R.res.mat.C2$FDR <- p.adjust(R.res.mat.C2$P, method="fdr")
write.table(R.res.mat.C2, file="~/Desktop/Eleni Paper/R.res.mat.C2.tsv", sep="\t")

colnames(NR.res.mat.C2) <- c("Estimate", "SE", "t-value", "P")
rownames(NR.res.mat.C2) <- rownames(agg.post.C2)
NR.res.mat.C2$FDR <- p.adjust(NR.res.mat.C2$P, method="fdr")
write.table(NR.res.mat.C2, file="~/Desktop/Eleni Paper/NR.res.mat.C2.tsv", sep="\t")

# Paired t-tests for all 
res.ttest.C2 <- matrix(data=NA, nrow=nrow(agg.pre.C2), ncol=4)
for (i in 1:nrow(agg.pre.C2)) {
  df <- data.frame("Y" = c( agg.pre.C2[i, -NR.pre], agg.post.C2[i, -NR.post]),
                   "X" = c(rep("Pre-treatment", length(agg.pre.C2[i, -NR.pre])), 
                           rep("Post-treatment", length(agg.post.C2[i, -NR.post]))))
  temp1 <- t.test(df$Y ~ df$X, paired=T)
  
  df <- data.frame("Y" = c( agg.pre.C2[i, NR.pre], agg.post.C2[i, NR.post]),
                   "X" = c(rep("Pre-treatment", length(agg.pre.C2[i, NR.pre])), 
                           rep("Post-treatment", length(agg.post.C2[i, NR.post]))))
  temp2 <- t.test(df$Y ~ df$X, paired=T)
  res.ttest.C2[i, ] <- c(temp1$statistic, temp1$p.value, temp2$statistic, temp2$p.value)
  
}
res.ttest.C2 <- data.frame(res.ttest.C2, row.names = rownames(agg.pre.C2))
colnames(res.ttest.C2) <- c("t_R", "p_R", "t_NR", "p_NR")
res.ttest.C2$FDR_R <- p.adjust(res.ttest.C2$p_R, method = 'fdr')
res.ttest.C2$FDR_NR <- p.adjust(res.ttest.C2$p_NR, method = 'fdr')

write.table(res.ttest.C2, file="~/Desktop/Eleni Paper/res.ttest.C2.tsv", sep="\t")
  
# 

df <- data.frame("Y" = c(agg.pre.C2["REACTOME_CELL_CYCLE", -NR.pre], agg.post.C2["REACTOME_CELL_CYCLE", -NR.post], 
                         agg.pre.C2["REACTOME_CELL_CYCLE", NR.pre], agg.post.C2["REACTOME_CELL_CYCLE", NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["REACTOME_CELL_CYCLE", -NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["REACTOME_CELL_CYCLE", -NR.post])),
                         rep("Pre-treatment", length(agg.pre.C2["REACTOME_CELL_CYCLE", NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["REACTOME_CELL_CYCLE", NR.post]))),
                 "Z" = c(rep("Reponder", length(agg.pre.C2["REACTOME_CELL_CYCLE", -NR.pre]) + length(agg.post.C2["REACTOME_CELL_CYCLE", -NR.post])),
                         rep("Non-Reponder", length(agg.pre.C2["REACTOME_CELL_CYCLE", NR.pre]) + length(agg.post.C2["REACTOME_CELL_CYCLE", NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
df$Z <- factor(df$Z, levels=c("Reponder", "Non-Reponder"))


ggplot(df, aes(x=Z, y=Y, fill=X)) +
  geom_boxplot(width=0.75, show.legend = T) +
  scale_fill_brewer(palette = "Accent") +
  #geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  labs(x = "Time", y="ssGSEA score")

df <- data.frame("Y" = c(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.pre], agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.post], 
                         agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.pre], agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.post])),
                         rep("Pre-treatment", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.post]))),
                 "Z" = c(rep("Reponder", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.pre]) + length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.post])),
                         rep("Non-Reponder", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.pre]) + length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
df$Z <- factor(df$Z, levels=c("Reponder", "Non-Reponder"))

ggplot(df, aes(x=Z, y=Y, fill=X)) +
  geom_boxplot(width=0.75, show.legend = T) +
  scale_fill_brewer(palette = "Accent") +
  #geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  labs(x = "Time", y="ssGSEA score")



df <- data.frame("Y" = c(agg.pre.C2["REACTOME_CELL_CYCLE", -NR.pre], agg.post.C2["REACTOME_CELL_CYCLE", -NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["REACTOME_CELL_CYCLE", -NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["REACTOME_CELL_CYCLE", -NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
t.test(df$Y ~ df$X, paired=T) #p-value = 0.004343
P1 <- ggplot(df, aes(x=X, y=Y, fill=X)) +
  geom_boxplot(show.legend = T, ) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", fill = "Time", y="ssGSEA score")

df <- data.frame("Y" = c( agg.pre.C2["REACTOME_CELL_CYCLE", NR.pre], agg.post.C2["REACTOME_CELL_CYCLE", NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["REACTOME_CELL_CYCLE", NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["REACTOME_CELL_CYCLE", NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
t.test(df$Y ~ df$X, paired=T) #p-value = 0.2049
P2 <- ggplot(df, aes(x=X, y=Y, fill=X)) +
  geom_boxplot(show.legend = T, ) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", fill = "Time", y="ssGSEA score")

df <- data.frame("Y" = c(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.pre], agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", -NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
t.test(df$Y ~ df$X, paired=T) #p-value = 0.007551
P3 <- ggplot(df, aes(x=X, y=Y, fill=X)) +
  geom_boxplot(show.legend = T, ) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", fill = "Time", y="ssGSEA score")

df <- data.frame("Y" = c( agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.pre], agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.post]),
                 "X" = c(rep("Pre-treatment", length(agg.pre.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.pre])), 
                         rep("Post-treatment", length(agg.post.C2["KEGG_OXIDATIVE_PHOSPHORYLATION", NR.post]))))
df$X <- factor(df$X, levels=c("Pre-treatment", "Post-treatment"))
t.test(df$Y ~ df$X, paired=T) #p-value = 0.8329
P4 <- ggplot(df, aes(x=X, y=Y, fill=X)) +
  geom_boxplot(show.legend = T, ) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_jitter(width = 0.1, show.legend = F, aes(color=Y)) +
  scale_color_viridis_c() +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", fill = "Time", y="ssGSEA score")


pdf("~/Desktop/Eleni Paper/Evero_res_Figures.pdf", height=4, width=5)
  print(P1)
  print(P2)
  print(P3)
  print(P4)
dev.off()


save.image("~/Desktop/Eleni Paper/Eleni_Paper.RData")
