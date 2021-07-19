# Rscript ccHPS.R  expression_file.txt

args<-commandArgs(T)
expDat <- read.table(args[1], head=T, sep="\t", 
		stringAsFactors = F, row.names =1)

# the expression of nine genes in ccHPS model

gene.min <- c("EFNA1","IER3", "ISG20","KLF7", "LDHC",
		"P4HA2", "PGM1",  "RBPJ",  "STC1")
samples_riskdat <- dplyr::select(expDat , gene.min)

# coefficient
coeff <- c(0.16972057, 0.13021757, -0.34117323, 0.13001858,
	-0.34334006, 0.40544016,  0.30079628, -0.17812931,
	 0.05886149)
score_df <- data.frame(gene.min, coeff)
names(score_df) <- c("name","coefficient")

#score_df

#    name coefficient
# 1 EFNA1  0.16972057
# 2  IER3  0.13021757
# 3 ISG20 -0.34117323
# 4  KLF7  0.13001858
# 5  LDHC -0.34334006
# 6 P4HA2  0.40544016
# 7  PGM1  0.30079628
# 8  RBPJ -0.17812931
# 9  STC1  0.05886149

samples_riskdat$score <- apply(samples_riskdat[,c(gene.min)], 1, function(x){
  x %*% score_df[c(gene.min),"coefficient" ]
})

write.table(samples_riskdat,"riskdat_ccHPS.tsv",sep="\t",quote = F)
