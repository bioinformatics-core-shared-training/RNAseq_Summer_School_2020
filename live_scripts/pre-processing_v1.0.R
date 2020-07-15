
# load library
library( DESeq2 )
library( tidyverse )

# read sample metadat
sampleinfo <- read_tsv( 'data/SampleInfo.txt' )
View(sampleinfo)
dim(sampleinfo)

# reads counts data
seqdata <- read_tsv( 'data/GSE60450_Lactation.featureCounts', comment = '#')
View(seqdata)
dim(seqdata)

# dplyr
# select : selects the columns
# filter : filter the rows
# reanme : raname the columns

# reformat
class(seqdata)

# %>% - MAC : COMMAND SHIFT M
# 
countdata <- seqdata %>% 
  column_to_rownames( 'Geneid') %>%  # column to row name
  rename_all( str_remove, '.bam') %>% # rename columns
  select( sampleinfo$Sample) %>% # select columns
  as.matrix()

head(countdata)

dim(countdata)

# filtering

keep <- rowSums( countdata) > 5

# see how many genes retained
table( keep, useNA='always')

# select rows
countdata <- countdata[ keep, ]
dim(countdata)

# box plots for read count distribution
head(countdata)
# get range
range( countdata)

# raw read counts box plot
boxplot( countdata)

# data trasformation
log2(0 + 1)

# log2 transformation
logcounts <- log2( countdata + 1 )
range(logcounts)
boxplot(logcounts)

# color box-plot based on groups
statusCol <- match( sampleinfo$Status, c( 'virgin', 'pregnant', 'lactate')) + 1 

boxplot( logcounts,
         xlab='',
         ylab='log2(counts + 1)',
         las=2,
         col=statusCol
)


# rlog transformation

rlogcounts <- rlog( countdata)

# box-plot from rlog counts
# color box-plot based on groups
statusCol <- match( sampleinfo$Status, c( 'virgin', 'pregnant', 'lactate')) + 1 

boxplot( rlogcounts,
         xlab='',
         ylab='log2(counts + 1)',
         las=2,
         col=statusCol
)
abline( h=median(rlogcounts), col='blue')

# PCA analysis
dim(countdata)

# Raw data : mean variance plot
plot(rowMeans(countdata), rowSds(countdata), main='Raw data', xlim=c(0,10000), ylim=c(0,5000))

# log2 trans. mean  var relationship
plot(rowMeans(logcounts), rowSds(logcounts), main='Log2 data')

# rlog data : mean-var plot
plot(rowMeans(rlogcounts), rowSds(rlogcounts), main='rlog data')

# PCA plot
pcDat <- prcomp( t(rlogcounts) )


library( ggfortify)

autoplot( pcDat)

barplot((pcDat$sdev^2 / sum(pcDat$sdev^2)) * 100)

# colour the PCA
autoplot( pcDat,
          data=sampleinfo,
          colour='CellType',
          shape='Status',
          size=5
)


library( ggrepel)

autoplot( pcDat,
          data=sampleinfo,
          colour='CellType',
          shape='Status',
          size=5
) +
  geom_text_repel( aes(x=PC1, y=PC2, label=Sample), box.padding = 1 )


# correct labels

sampleinfo <- sampleinfo %>% 
  mutate(CellType=ifelse(Sample=="MCL1.DG", "basal", CellType)) %>% 
  mutate(CellType=ifelse(Sample=="MCL1.LA", "luminal", CellType))

write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")
