##GSEA TCGA
# step 1; differential analysis based on limma packages between low and high score groups
library(limma)
#load expression matrix
tcga_exp <- read.table("tcga_tumor.exp.txt",head=T,sep="\t",
                             stringsAsFactors = F,row.names = 1)

#load sample groups
tcga_samples <- all_riskdat
tcga_samples$High_score <- ifelse(tcga_samples$groups == 'high-score',0,1)
tcga_samples$Low_score <- ifelse(tcga_samples$groups == 'low-score',0,1)
used_samples <- intersect(row.names(tcga_samples),colnames(tcga_exp))
  
tcga_groups <- tcga_samples[used_samples, c("High_score","Low_score")]
tcga_limma_exp <- tcga_exp[,colnames(tcga_exp) %in% row.names(tcga_groups)]

# limma
design1 <- tcga_groups
fit <- lmFit(tcga_limma_exp,design1)

contrast.matrix <- makeContrasts( Low_score - High_score,levels=design1) 
fit2 <- contrasts.fit(fit,contrast.matrix)

fit2 <- eBayes(fit2)

all_diff_tcga <- topTable(fit2, adjust.method = 'fdr', coef = 1, p.value = 1, lfc = log(1,2), number = 60000, sort.by = 'logFC')

#GSEA
library(clusterProfiler)
library(enrichplot)
all_diff_tcga_sort <-all_diff_tcga[with(all_diff_tcga,order(-logFC)),]
write.table(all_diff_tcga_sort,"tcga_gene_rank.txt",sep="\t",quote=F)

geneList <- all_diff_tcga_sort$logFC
names(geneList) = row.names(all_diff_tcga_sort)

#GO
set.seed(1234)
gogmt <- read.gmt("c5.go.bp.v7.2.symbols.gmt")
GO <- GSEA(geneList,TERM2GENE = gogmt,
           pvalueCutoff = 0.05,nPerm = 20000,seed = TRUE)
GOdf <- data.frame(GO)
write.table(GOdf,"supplementaty table GSEA result GO in TCGAdataset.txt",sep="\t",
            row.names = F,quote = F)

# High risks

pdf("Fig.TCGA_GSEA_GO_Highrisk.pdf",height = 6, width = 9)
gseaplot2(GO,geneSetID = c("GO_KERATINIZATION","GO_KERATINOCYTE_DIFFERENTIATION",
                           "GO_CORNIFICATION","GO_LYMPHOCYTE_MIGRATION",
                           "GO_RESPONSE_TO_INTERLEUKIN_1"),
          pvalue_table = F,subplots = 1:2)
dev.off()

# Low risks

pdf("Fig.TCGA_GSEA_GO_Lowrisk.pdf",height = 6, width = 9)
gseaplot2(GO,geneSetID = c("GO_B_CELL_ACTIVATION","GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                           "GO_B_CELL_MEDIATED_IMMUNITY","GO_LYMPHOCYTE_MEDIATED_IMMUNITY",
                           "GO_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY"),
          pvalue_table = F,subplots = 1:2)
dev.off()

# merge all

pdf("Fig.TCGA_GSEA_GO_ALL.pdf",height = 9, width = 9)
gseaplot2(GO,geneSetID = c("GO_KERATINIZATION","GO_KERATINOCYTE_DIFFERENTIATION",
                           "GO_CORNIFICATION","GO_LYMPHOCYTE_MIGRATION",
                           "GO_RESPONSE_TO_INTERLEUKIN_1","GO_B_CELL_ACTIVATION","GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                           "GO_B_CELL_MEDIATED_IMMUNITY","GO_LYMPHOCYTE_MEDIATED_IMMUNITY",
                           "GO_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY"),
          pvalue_table = F,subplots = 1:2)
dev.off()

#KEGG
set.seed(1234)
kegg_gmt <- read.gmt("c2.cp.kegg.v7.2.symbols.gmt")
KEGG <- GSEA(geneList,TERM2GENE = kegg_gmt,
           pvalueCutoff = 0.05,nPerm = 20000,seed = TRUE)
KEGGdf <- data.frame(KEGG)
write.table(KEGGdf,"supplementaty table GSEA result TCGAdataset_KEGG.txt",sep="\t",
            row.names = F,quote = F)


pdf("Fig.TCGA_GSEA_KEGG.pdf",height = 9, width = 9)
gseaplot2(KEGG,geneSetID = c("KEGG_PROTEASOME","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_RIBOSOME","KEGG_OLFACTORY_TRANSDUCTION"),pvalue_table = F,subplots = 1:2)
dev.off()

#hallmark gene sets

set.seed(1234)
hallmark_gmt <- read.gmt("h.all.v7.4.symbols.gmt")
hallmark <- GSEA(geneList,TERM2GENE = hallmark_gmt,
           pvalueCutoff = 0.05,nPerm = 20000,seed = TRUE)
hallmarkdf <- data.frame(hallmark)
write.table(hallmarkdf,"supplementaty table GSEA result TCGAdataset_hallmark.txt",sep="\t",
            row.names = F,quote = F)

pdf("Fig.TCGA_GSEA_HALLMARK.pdf",height = 9, width = 9)
gseaplot2(hallmark,geneSetID = c("HALLMARK_KRAS_SIGNALING_DN",
"HALLMARK_INFLAMMATORY_RESPONSE",
"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
"HALLMARK_APOPTOSIS",
"HALLMARK_MYOGENESIS",
"HALLMARK_PROTEIN_SECRETION",
"HALLMARK_FATTY_ACID_METABOLISM"),pvalue_table = F,subplots = 1:2)
dev.off()
