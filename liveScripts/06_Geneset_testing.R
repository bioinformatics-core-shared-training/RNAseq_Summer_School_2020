library(clusterProfiler)
library(tidyverse)

search_kegg_organism('mmu', by='kegg_code')

load("RNAseq/Robjects/Annotated_Results_LvV.RData")

sigGenes <- shrinkLvV %>% 
  drop_na(Entrez, FDR) %>% 
  filter(FDR < 0.05 & abs(logFC) > 1) %>% 
  pull(Entrez)

sigGenes

kk <- enrichKEGG(gene = sigGenes, organism = 'mmu')
head(kk, n=10)

# visualise the pathway

browseKEGG(kk, 'mmu03320')

# use pathview
library(pathview)

logFC <- shrinkLvV[match(sigGenes, shrinkLvV$Entrez), "logFC"]
names(logFC) <- sigGenes

pathview(gene.data = logFC,
         pathway.id = "mmu03320",
         species = 'mmu',
         limit = list(gene = 5, cpd = 1))

# GSEA 

library(fgsea)

rankData <- shrinkLvV %>% 
  drop_na(Entrez) %>% 
  pull(logFC)

names(rankData) <- shrinkLvV %>% 
  drop_na(Entrez) %>% 
  pull(Entrez)

head(rankData)

# load in the pathways

load("RNAseq/Robjects/mouse_H_v5.RData")

str(Mm.H)

# conduct the GSEA

fgseaRes <- fgsea(Mm.H,
                  rankData,
                  minSize = 15,
                  maxSize = 500)

head(fgseaRes)

fgseaRes %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, -padj)

# plot the enrichment score
plotEnrichment(Mm.H[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], rankData)


