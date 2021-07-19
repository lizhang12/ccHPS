#load expression data
tcga_expression <- read.table("tcga_tumor.exp.txt",head=T,sep="\t",
                              stringsAsFactors = F,row.names = 1)
# low- and high-index groups
#sample_info <- read.table("tcga_all_sample_groups.txt",head=T,sep="\t",
#                          stringsAsFactors = F,row.names = 1)

sample_info <- all_riskdat
high_samples <- row.names(sample_info[which(sample_info$groups == 'high-score'),])
low_samples <- row.names(sample_info[which(sample_info$groups == 'low-score'),])

# load immune genes
immune_file <- read.table("LM22.txt", head=T,sep="\t",
		stringsAsFactors = F,row.names = 1)

immune_genes <- row.names(immune_file)
tcga_high_exp <- tcga_expression[row.names(tcga_expression)%in% immune_genes,
                                names(tcga_expression) %in% high_samples]
tcga_low_exp <- tcga_expression[row.names(tcga_expression)%in% immune_genes,
                                names(tcga_expression) %in% low_samples]
write.table(tcga_high_exp,"tcga_cibersort_high.txt",row.names = T,sep="\t",quote=F)
write.table(tcga_low_exp,"tcga_cibersort_low.txt",row.names = T,sep="\t",quote=F)


source("CIBERSORT.R")
library(ggpubr)

##### high-ccHPS ##########
ciber_high <- CIBERSORT("LM22.txt",
                        'tcga_cibersort_high.txt',perm = 100, QN = TRUE)
write.table(ciber_high, "ciber_high.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

# convert matrix to ggplot2 input dataframe.
ciber_high.plot <- as.data.frame(as.table(as.matrix(ciber_high[,c(1:22)])))
names(ciber_high.plot) <- c("PATIENT_ID","CellType","Composition")

# Barplot of cell componment of each sample using hclust
high.sample.index <- hclust(dist(ciber_high[,c(1:22)]), method = "ward.D")$order
high.sample.order <- row.names(ciber_high)[high.sample.index]

highP<- ggbarplot(
  ciber_high.plot,
  x = "PATIENT_ID",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType",
  order = high.sample.order
) +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.line=element_blank(),
    legend.position = "bottom",
    
  )+
  labs(title = "High-ccHPS")

ggsave(
  filename = "Fig5a_high_cibersort.pdf",highP,device = "pdf",
  height = 5, width = 10
)

##### low-ccHPS ##########

ciber_low <- CIBERSORT("LM22.txt",
                       'tcga_cibersort_low.txt', perm = 100, QN = TRUE)

write.table(ciber_low, "ciber_low.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

# convert matrix to ggplot2 input dataframe.
ciber_low.plot <- as.data.frame(as.table(as.matrix(ciber_low[,c(1:22)])))
names(ciber_low.plot) <- c("PATIENT_ID","CellType","Composition")

# Barplot of cell componment of each sample using hclust
low.sample.index <- hclust(dist(ciber_low[,c(1:22)]), method = "ward.D")$order
low.sample.order <- row.names(ciber_low)[low.sample.index]

lowP<- ggbarplot(
  ciber_low.plot,
  x = "PATIENT_ID",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType",
  order = low.sample.order
) +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.line=element_blank(),
    legend.position = "bottom",
    
  )+
  labs(title = "Low-ccHPS")

ggsave(
  filename = "Fig5a_low_cibersort.pdf",lowP,device = "pdf",
  height = 5, width = 10
)
highP
lowP

# box plot group by Group
ciber_high_conv <- read.table("ciber_high.results.output.txt",head=T,sep="\t",
                             row.names = 1,stringsAsFactors = F)
ciber_high_sig <- subset(ciber_high_conv, ciber_high_conv$P.value < 0.05)
ciber_high_sig.plot <- as.data.frame(as.table(as.matrix(ciber_high_sig[,c(1:22)])))
names(ciber_high_sig.plot) <- c("PATIENT_ID","CellType","Composition")


ciber_low_conv <- read.table("ciber_low.results.output.txt",head=T,sep="\t",
                             row.names = 1,stringsAsFactors = F)
ciber_low_sig <- subset(ciber_low_conv, ciber_low_conv$P.value < 0.05)
ciber_low_sig.plot <- as.data.frame(as.table(as.matrix(ciber_low_sig[,c(1:22)])))
names(ciber_low_sig.plot) <- c("PATIENT_ID","CellType","Composition")


ciber_high_sig.plot$Group <- "High-score"
ciber_low_sig.plot$Group <- "Low-score"
merge_sig.plot <- rbind(ciber_high_sig.plot,ciber_low_sig.plot)
write.table(merge.plot,"ciber_boxplot_input_sig.txt",sep="\t",quote=F,row.names = F)

library(ggthemes)
symnum.args <- list(cutpoints = c(0,0.01, 1), symbols = c("%", "ns"))

p1tmp<- ggboxplot(
  merge_sig.plot,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Group",
  xlab = "",
  ylab = "Cell composition",
  #main = "TME Cell composition group by Group",
  outlier.shape = NA,
  ylim = c(0,0.8)
  #palette = c( "#ED553B","#646970")
) +
  stat_compare_means(aes(group = Group),
    label = "p.signif",
    method = "t.test",
    #ref.group = ".all.",
    
    hide.ns = T
  ) +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))

 p1<- p1tmp+ theme(axis.text.x =element_text(size=10),
        axis.title.y =element_text(size=10))
 ggsave(
  filename = "Fig_compare_cibersort_sig_chagecolor.pdf",p1,device = "pdf",
  height = 6, width = 10
)
