# load the expression matrix of high-score samples
# gene in the row and samples in the column.
mcp_high_input <- read.table("tcga_expression_high_risk.txt",head=T,sep="\t",
                             stringsAsFactors = F)
mcp_estimate_high <- MCPcounter::MCPcounter.estimate(mcp_high_input,
                                                     featuresType = c('HUGO_symbols'))
high_mcp <- t(as.data.frame(mcp_estimate_high))
high_mcp.plot <- as.data.frame(as.table(as.matrix(high_mcp)))
names(high_mcp.plot) <- c("PATIENT_ID","CellType","Abundance")
high_mcp.plot$Group <- "High-score"

# low-score samples
mcp_low_input <- read.table("tcga_expression_low_risk.txt",head=T,sep="\t",
                             stringsAsFactors = F)
mcp_estimate_low <- MCPcounter::MCPcounter.estimate(mcp_low_input,
                                                     featuresType = c('HUGO_symbols'))
low_mcp <- t(as.data.frame(mcp_estimate_low))
low_mcp.plot <- as.data.frame(as.table(as.matrix(low_mcp)))
names(low_mcp.plot) <- c("PATIENT_ID","CellType","Abundance")
low_mcp.plot$Group <- "Low-score"

# box plot group by Group
merge_mcp.plot <- rbind(high_mcp.plot,low_mcp.plot)

write.table(merge_mcp.plot,"MCPcount_boxplot.txt",sep="\t",quote=F,row.names = F)
library(ggthemes)
symnum.args <- list(cutpoints = c(0,0.01, 1), symbols = c("%", "ns"))

ptemp<- ggboxplot(
  merge_mcp.plot,
  x = "CellType",
  y = "Abundance",
  color = "black",
  fill = "Group",
  xlab = "",
  ylab = "Cell Abundance",
  #main = "TME Cell composition group by Group",
  outlier.shape = NA
  #ylim = c(0,2),
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

 p2<- ptemp + theme(axis.text.x =element_text(size=10),
        axis.title.y =element_text(size=10))
 ggsave(
  filename = "Fig_compare_MCPcounter.pdf",p2,device = "pdf",
  height = 6, width = 10
)
