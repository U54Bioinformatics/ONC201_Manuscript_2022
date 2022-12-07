setwd("~/Desktop/Eleni Paper/")
list.files()

ONC.exp <- read.table("GSE212369_ONC201_time_course_RNAseq.txt", header=T, as.is = T)
ONC.annot <- read.delim("GSE212369_series_matrix.txt", header = F, as.is = T, fill = T, sep = "\t")

ONC.exp <- ONC.exp[-which(duplicated(ONC.exp$Gene)), ]
rownames(ONC.exp) <- ONC.exp$Gene
ONC.exp <- as.matrix(ONC.exp[, -c(1:2)])

### Run GSVA in ssGSEA with hallmark and C2 gene sets and output to file -------------------------------------------
library(GSEABase)
library(GSVA)
library(BiocParallel)

c2set <- GSEABase::getGmt("~/Desktop/GeneSets/c2.all.v7.0.symbols.gmt")
Hset <- GSEABase::getGmt("~/Desktop/GeneSets/H.all.v7.0.symbols.gmt")

ONC.exp.H <- gsva(ONC.exp, Hset, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())
ONC.exp.C2 <- gsva(ONC.exp, c2set, min.sz=10, verbose=TRUE, kcdf="Gaussian", method="ssgsea", parallel.sz = 5, BPPARAM=MulticoreParam())

colnames(ONC.exp)
ONC.time <- c(rep(0, 6), 
              rep(3, 3),
              rep(6, 3),
              rep(12, 3),
              rep(24, 3))

library(mgcv)
summary(gam(ONC.exp.H[1,] ~ s(ONC.time, k = 4)))
ggplot(data.frame("X"=ONC.time, "Y"=ONC.exp.H[1,]), aes(y=Y, x=X)) + 
  geom_smooth(method = "gam", formula =  y ~ s(x, bs = "re", k = 4)) +
  geom_point(aes(color = Y)) +
  theme_classic(base_size = 18)

res.mat.H <- data.frame(matrix(data=NA, nrow=nrow(ONC.exp.H), ncol=4))
for (i in 1:nrow(ONC.exp.H)) {
  res <- summary(glm(ONC.exp.H[i,] ~ ONC.time))
  res.mat.H[i,] <- res$coefficients[2,]
}
colnames(res.mat.H) <- c("Estimate", "SE", "t-value", "P")
rownames(res.mat.H) <- rownames(ONC.exp.H)
res.mat.H$FDR <- p.adjust(res.mat.H$P, method="fdr")
write.table(res.mat.H, file="~/Desktop/Eleni Paper/res.mat.H.tsv", sep="\t")

res.mat.C2 <- data.frame(matrix(data=NA, nrow=nrow(ONC.exp.C2), ncol=4))
for (i in 1:nrow(ONC.exp.C2)) {
  res <- summary(glm(ONC.exp.C2[i,] ~ ONC.time))
  res.mat.C2[i,] <- res$coefficients[2,]
}
colnames(res.mat.C2) <- c("Estimate", "SE", "t-value", "P")
rownames(res.mat.C2) <- rownames(ONC.exp.C2)
res.mat.C2$FDR <- p.adjust(res.mat.C2$P, method="fdr")
write.table(res.mat.C2, file="~/Desktop/Eleni Paper/res.mat.C2.tsv", sep="\t")

res.mat.exp <- data.frame(matrix(data=NA, nrow=nrow(ONC.exp), ncol=4))
for (i in 1:nrow(ONC.exp)) {
  res <- summary(glm(ONC.exp[i,] ~ ONC.time))
  res.mat.exp[i,] <- res$coefficients[2,]
}
colnames(res.mat.exp) <- c("Estimate", "SE", "t-value", "P")
rownames(res.mat.exp) <- rownames(ONC.exp)
res.mat.exp$FDR <- p.adjust(res.mat.exp$P, method="fdr")
write.table(res.mat.exp, file="~/Desktop/Eleni Paper/res.mat.exp.tsv", sep="\t")

path.list <- c("REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
               "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
               "BIOCARTA_ERK_PATHWAY",
               "BIOCARTA_AKT_PATHWAY",
               "REACTOME_CELL_CYCLE")

pdf("~/Desktop/Eleni Paper/ONC201_Time_Figures.pdf", height=4, width=4.5)
for (i in 1:length(path.list)) {
  x <- match(path.list[i], rownames(ONC.exp.C2))
  P <- ggplot(data.frame("X"=ONC.time, "Y"=ONC.exp.C2[x,]), aes(y=Y, x=X)) + 
        geom_smooth(method = "glm", color="darkred", fill="salmon") +
        geom_point(aes(color = Y), size=4,show.legend = F) +
        theme_classic(base_size = 20) +
        scale_color_viridis_c(alpha=0.75) +
        labs(x="Time (hours)", y="ssGSEA score")
  print(P)
}
dev.off()


ggplot(data.frame("X"=ONC.time, "Y"=ONC.exp.C2["KEGG_OXIDATIVE_PHOSPHORYLATION",]), aes(y=Y, x=X)) + 
  geom_smooth(method = "glm", color="darkred", fill="salmon") +
  geom_point(aes(color = Y), size=4,show.legend = F) +
  theme_classic(base_size = 20) +
  scale_color_viridis_c(alpha=0.75) +
  labs(x="Time (hours)", y="ssGSEA score")
