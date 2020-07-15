library(EnsDb.Mmusculus.v79)
library(DESeq2)
library(tidyverse)

load("../Course_Materials/Robjects/DE.RData")

# annotation

columns(EnsDb.Mmusculus.v79)
keytypes(EnsDb.Mmusculus.v79)

# query

ourCols <- c("SYMBOL", "GENEID", "ENTREZID")
ourKeys <- rownames(resLvV)[1:1000]

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
head(annot)
dim(annot)
annot %>%
  add_count(GENEID) %>%
  dplyr::filter(n>1)

# challenge 1

ourCols <- c("SYMBOL", "GENEID", "ENTREZID", "GENEBIOTYPE")
ourKeys <- rownames(resLvV)

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
dim(annot)
multiples <- annot %>%
  add_count(GENEID) %>%
  dplyr::filter(n>1)
length(unique(multiples$SYMBOL))

load("../Course_Materials/Robjects/Ensembl_annotations.RData")

colnames(ensemblAnnot)

# adding annotation to results

annotLvV <- as.data.frame(resLvV) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC = log2FoldChange, FDR = padj)

write_tsv(annotLvV, "../Course_Materials/data/VirginVsLactating_Results_Annotated.txt")

annotLvV %>%
  arrange(FDR) %>%
  head(10)

# Visualisation

ddsShrink <- lfcShrink(ddsObj, coef = "Status_lactate_vs_virgin")

shrinkLvV <- as.data.frame(ddsShrink) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC = log2FoldChange, FDR = padj)

hist(shrinkLvV$pvalue)

# MA plots

plotMA(ddsShrink, alpha = 0.05)

cutoff <- sort(shrinkLvV$pvalue)[10]
shrinkLvV <- shrinkLvV %>%
  mutate(TopGeneLabel = ifelse(pvalue<=cutoff, Symbol, ""))

ggplot(shrinkLvV, aes(x = log2(baseMean), y = logFC)) +
  geom_point(aes(colour = FDR < 0.05), shape = 20, size = 0.5) +
  geom_text(aes(label = TopGeneLabel)) +
  labs(x = "mean of normalised counts", y = "log fold change")

# challenge 2

shrinkLvV <- shrinkLvV %>%
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(shrinkLvV, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = FDR < 0.05), size = 1)

# heatmaps

library(ComplexHeatmap)
library(circlize)

sigGenes <- as.data.frame(shrinkLvV) %>%
  top_n(150, wt = -FDR) %>%
  pull("GeneID")

plotDat <- vst(ddsObj)[sigGenes,] %>%
  assay()
z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

myPalette <- c("red3","ivory", "blue3")
myRamp <- colorRamp2(c(-2,0,2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        cluster_columns = FALSE)

# cutting the tree

hcDat <- hclust(dist(z.mat))
cutGroups <- cutree(hcDat, h = 4)

ha1 <- HeatmapAnnotation(df = colData(ddsObj)[,c("CellType","Status")])

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        split = cutGroups, 
        rect_gp = gpar(col = "darkgrey", lwd=0.5),
        top_annotation = ha1)

