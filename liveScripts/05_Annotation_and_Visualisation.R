library(EnsDb.Mmusculus.v79)
library(DESeq2)
library(tidyverse)

# load in our data

load("RNAseq/Robjects/DE.RData")

resLvV

# Query the Database

columns(EnsDb.Mmusculus.v79)

keytypes(EnsDb.Mmusculus.v79)

# let's set up the query

ourCols <- c("SYMBOL", "GENEID", "ENTREZID")
ourKeys <- rownames(resLvV)
head(ourKeys)
ourKeys <- rownames(resLvV)[1:1000]

# run the query
annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               key="GENEID")

head(annot)

length(unique(annot$ENTREZID))

dim(annot)
length(unique(annot$GENEID))

annot %>% 
  add_count(GENEID) %>% 
  filter(n>1) %>% 
  head()

# Challenge 1

columns(EnsDb.Mmusculus.v79)
ourCols <- c("SYMBOL", "GENEID", "ENTREZID", "GENEBIOTYPE")
ourKeys <- rownames(resLvV)

# run the query
annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               key="GENEID")
dim(annot)

annot %>% 
  add_count(GENEID) %>% 
  filter(n>1) %>% 
  distinct(GENEID) %>% 
  count()

## One we made ealier

load("RNAseq/Robjects/Ensembl_annotations.RData")

# annotate our results table

annotLvV <- resLvV %>% 
  as.data.frame() %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, by="GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)
head(annotLvV)  

write_tsv(annotLvV, "RNAseq/results/LactateVsVirgin_Results_Annotated.txt")

# Visualising our data

# P-value histogram

hist(resLvV$pvalue)

# Shrinking the fold changes

ddsShrink <- lfcShrink(ddsObj, coef = "Status_lactate_vs_virgin")
ddsShrink

shrinkLvV <- ddsShrink %>% 
  as.data.frame() %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, by = "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)
head(shrinkLvV)

# MA plots

par(mfrow=c(1, 2))
plotMA(resLvV, alpha = 0.05)
plotMA(ddsShrink, alpha = 0.05)

# ggplot2

ggplot(shrinkLvV, mapping = aes(x = log2(baseMean), y = logFC)) +
  geom_point(aes(colour = FDR < 0.05), shape = 20, size = 0.05) +
  geom_text(data=~top_n(.x, 10, wt=-FDR), aes(label = Symbol)) +
  labs(x = "Mean of normalised counts", y= "log2(Fold Change)")

# Challenge - volcano plot

shrinkLvV %>% 
  mutate(`-log10(pvalue)` = -log10(pvalue)) %>% 
  ggplot(aes(x = logFC, y = `-log10(pvalue)`)) +
    geom_point(aes(colour = FDR < 0.05), size = 1)

# Heatmap with ComplexHeatmap

library(ComplexHeatmap)
library(circlize)

sigGenes <- shrinkLvV %>% 
  top_n(150, wt = -FDR) %>% 
  pull("GeneID")

# filter and normalise the counts
plotDat <- vst(ddsObj)[sigGenes, ] %>% 
  assay()

# z scaling

z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE))

# colour palette
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)
myRamp

# generate the heatmap
Heatmap(z.mat,
        name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        cluster_columns = FALSE)

# create annotation

ha1 <- HeatmapAnnotation(df = colData(ddsObj)[,c("CellType", "Status")])


Heatmap(z.mat,
        name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        top_annotation = ha1,
        split = 7,
        rect_gp = gpar(col="darkgrey", lwd=0.5))





























