#module load R/4.0.0
library(biomaRt)
library(AnnotationHub)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyr)
library(GenomicFeatures)
library(Repitools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
prom1 <- promoters(genes(txdb), upstream = 500, downstream = 100)
gc<- gcContentCalc(prom1, organism=Hsapiens, verbose=TRUE)
gc <- as.data.frame(gc)
gc$gene_id <- as.vector(prom1@elementMetadata$gene_id)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

hg_gene_ids  <- as.vector(gc$gene_id)

foo <- getBM(attributes=c('ensembl_gene_id',
                           'percentage_gene_gc_content',
                           'hgnc_symbol',
                           'start_position','end_position'),
              filters = 'ensembl_gene_id',
              values = genes,
              mart = ensembl)
foo$gene_length=foo$end_position - foo$start_position
gc$gene_id <- as.numeric(gc$gene_id)
gc_content <- dplyr::full_join(foo,gc,by=c("entrezgene_id"="gene_id"))
write.table('GC_content_promoter.tsv', gc_content, sep='\t')
