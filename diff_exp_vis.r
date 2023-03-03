##########################################################################################
##########################################################################################
##########################################################################################

# Autor: Nikita Vaulin, BI, Skoltech. vaulin@ro.ru
# Date: 26.2.2023
# Based on E. Khrameeva 'Omics Data Analysis' skoltech course materials
# Created for the BI practicals project on differential expression

# For the first time you may need to run manually 'ah <- AnnotationHub()' line (№ 54)
# to agree with database downloading
# Also ensure that packages installed

##########################################################################################
##########################################################################################
##########################################################################################

BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("pathview")

library(dplyr)
library(tidyverse)
library(ensembldb)
library(AnnotationHub)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(pathview)
library(enrichplot)


### # set constants
padj.cutoff <- 0.05
lfc.cutoff <- 0.6


### upload file from bash pipes
#diffgenes = read.table("stdin", header=TRUE, dec='.', sep="\t", row.names=1 )

### when running manually in RSturio, upload by name
diffgenes = read.table("result.txt", header=TRUE, dec='.', sep="\t", row.names=1)


##########################################################################################
##########################################################################################
#                                    PREPARE THE DATA                                    #              
##########################################################################################
##########################################################################################

### getting Sac Cer genes annotation from web
ah <- AnnotationHub()
yeast_ens <- query(ah, c("Saccharomyces cerevisiae", "EnsDb"))
yeast_ens <- yeast_ens[[tail(yeast_ens$ah_id, 1)]]


### creating a table to compare transcripts codes in data with real genes names
genedb <- genes(yeast_ens, return.type = "data.frame")
genedb <- dplyr::select(genedb, gene_id, symbol)
txdb <- transcripts(yeast_ens, return.type = "data.frame")
txdb <- txdb %>% select(tx_id, gene_id)
tx2gene <- right_join(txdb, genedb, by='gene_id')

### fixing genes names in our data
diffgenes <- diffgenes %>% mutate(gene_id = str_remove(id, 'gene-'))

###  attaching real genes names
diffgenes.an <- left_join(diffgenes, genedb, by='gene_id')

rm(list = c('diffgenes', 'genedb', 'txdb', 'tx2gene'))
gc()

##########################################################################################
##########################################################################################
#                                   BUILD SOME PLOTS                                     #
##########################################################################################
##########################################################################################

dir.create('plots', mode = "0777")
setwd('plots/')

### build volcano
volcano <-  EnhancedVolcano(diffgenes.an,
                  lab = diffgenes.an$symbol,
                  x = 'log2FoldChange', y = 'padj',
                  title = NULL, subtitle = NULL,
                  pCutoff = padj.cutoff, FCcutoff = lfc.cutoff)

ggsave('Volcano_plot.png', volcano, height = 7, width = 10)

rm(list = c('volcano'))


### dotplots to find most momular GO terms in significant genes
sig_genes <- diffgenes.an %>%dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

ego_bp <- enrichGO(gene = sig_genes$gene_id, universe = diffgenes.an$gene_id,
                keyType = "ENSEMBL", OrgDb = org.Sc.sgd.db,
                pAdjustMethod = "BH", pvalueCutoff = padj.cutoff,
                ont = "BP")
ego_mf <- enrichGO(gene = sig_genes$gene_id, universe = diffgenes.an$gene_id,
                   keyType = "ENSEMBL", OrgDb = org.Sc.sgd.db,
                   pAdjustMethod = "BH", pvalueCutoff = padj.cutoff,
                   ont = "MF")
ego_cc <- enrichGO(gene = sig_genes$gene_id, universe = diffgenes.an$gene_id,
                   keyType = "ENSEMBL", OrgDb = org.Sc.sgd.db,
                   pAdjustMethod = "BH", pvalueCutoff = padj.cutoff,
                   ont = "CC")

dot_bp <- dotplot(ego_bp, showCategory = 10)
dot_mf <- dotplot(ego_mf, showCategory = 10)
dot_cc <- dotplot(ego_cc, showCategory = 10)

ggsave('dotplot_BP.png', dot_bp, height = 9, width = 9)
ggsave('dotplot_MF.png', dot_mf, height = 9, width = 9)
ggsave('dotplot_CC.png', dot_cc, height = 9, width = 9)

rm(list = c('ego_bp','ego_mf', 'ego_cc'))
rm(list = c('dot_bp','dot_mf', 'dot_cc'))
gc()

##########################################################################################
##########################################################################################
#                                  VISUALISE KEGG PATHWAYS                               #
##########################################################################################
##########################################################################################

# get the KEGG data and ensembl id's


entrez = bitr(sig_genes$gene_id,
              fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Sc.sgd.db")
sig_entrez <- left_join(sig_genes, entrez, by=c('gene_id'='ENSEMBL'))


ekegg <- enrichKEGG(gene = sig_entrez$gene_id, universe=diffgenes.an$gene_id, 
                    organism = 'sce', 
                    pvalueCutoff = 1, # padj.cutoff, 
                    qvalueCutoff = 1000) #lfc.cutoff)


write.table(ekegg[,c('ID', 'p.adjust', 'qvalue')], file = 
              "KEGG_pathways_pvalues.txt", sep = "\t", row.names = F)
# extract the significant ensembl id's
fc <- sig_entrez$log2FoldChange
names(fc) <- sig_entrez$ENTREZID
fc.sorted <- sort(fc, decreasing = TRUE)

dir.create('kegg_pathways', mode = "0777")
setwd('kegg_pathways/')

# prints the pathways as 'sce00010.pathway.png' (and downloads some files from web)
for(pathway in unique(c(ekegg$ID, 'sce00010'))){
    pathview(gene.data = fc.sorted, species = "sce",
           pathway.id = pathway, high = "green", low = 'red')
    }


diffgenes.an %>%dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff) %>%  nrow
diffgenes.an %>%dplyr::filter(abs(log2FoldChange) >= lfc.cutoff) %>%  nrow
