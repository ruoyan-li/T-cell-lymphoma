setwd('/Users/ruoyanli/Desktop/ALL/TCL_temp')

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
data(geneList, package="DOSE")

df <- read.table('malignant_TLBL_vs_DN_panfetal_pseudobulk_DEGs_allup.xls',
                 header = T, sep='\t')
df %>% head
de1 <- df$Gene[df$padj < 0.01]

########## KEGG
de1_ID <- bitr(de1, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
ekegg = enrichKEGG(de1_ID, pvalueCutoff=0.05, use_internal_data=T)
dotplot(ekegg, showCategory=20)

library(ggplot2)
library(forcats)
library(DOSE)

ggplot(ekegg, showCategory = 20,
       aes(GeneRatio, fct_reorder(Description, GeneRatio))) + 
  #geom_segment(aes(xend=0, yend = Description, color=p.adjust)) + 
  geom_point(aes(color=p.adjust, size = Count)) + 
  scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) + 
  scale_size_continuous(range=c(2, 10)) +
  #coord_cartesian(xlim = c(0, 0.08))+
  theme_dose(12) + 
  xlab("GeneRatio") + 
  ggtitle("KEGG enrichment")

gene_id <- ekegg$geneID[ekegg$Description=='T cell receptor signaling pathway']
gene_id <- strsplit(gene_id,'/')[[1]]
bitr(gene_id, "ENTREZID", "SYMBOL", OrgDb = org.Hs.eg.db)

#### pathway view 
library(pathview)
de1_ID <- bitr(de1, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

hsa04660 <- pathview(gene.data = de1_ID,
                     pathway.id = "hsa04660", ## TCR signaling
                     species = "hsa", kegg.native = FALSE)
