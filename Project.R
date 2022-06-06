rm(list=ls())

lusc <- read.table(file = 'LUSC1', sep = '\t', header = TRUE, row.names = 1)
lusc[2] <- read.table(file = 'LUSC2', sep = '\t', header = TRUE, row.names = 1) 
lusc[3] <- read.table(file = 'LUSC3', sep = '\t', header = TRUE, row.names = 1)

luad <- read.table(file = 'LUAD1', sep = '\t', header = TRUE, row.names = 1)
luad[2] <- read.table(file = 'LUAD2', sep = '\t', header = TRUE, row.names = 1) 
luad[3] <- read.table(file = 'LUAD3', sep = '\t', header = TRUE, row.names = 1)

#merge into the countTable
countTable <- merge(lusc, luad, by = "row.names", all = TRUE)
rownames(countTable) <- countTable$Row.names
countTable$Row.names <- NULL
countTable <- data.matrix(countTable, rownames.force = NA)
countTable <- countTable[-1:-5,]
head(countTable)

# calling DE genes with edgeR
library("edgeR")
# remove .* numbers
ensembl <- gsub("\\..*", "", rownames(countTable))
rownames(countTable) <- ensembl
head(countTable)

group <- as.factor(rep(c("LUSC","LUAD"), c(3,3))) 
#countTable[is.na(countTable)] = 0
y <- DGEList(counts=countTable)
y$samples$group <- group
y

#ensembl <- c("TCGA-21-1079-01A-01R-0692-07", "TCGA-22-4599-01A-01R-1443-07", "TCGA-37-4132-01A-01R-1100-07", "TCGA-95-7944-01A-11R-2187-07", "TCGA-78-8648-01A-11R-2403-07", "TCGA-50-6595-01A-12R-1858-07")
ensembl <- c("0692-07", "1443-07", "1100-07", "2187-07", "2403-07", "1858-07")
barplot(y$samples$lib.size*1e-6, names=ensembl, ylab="Library size (million)", las=2, col=rep(c("blue", "light blue"), each=3))
legend("top", legend = c("LUSC", "LUAD"), fill=c("blue", "light blue"), bty = "n", inset = c(-0.05,-0.1), xpd=TRUE, horiz=T)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

y <- calcNormFactors(y, method = "TMM")
y

logcpm <- cpm(y, log=TRUE)

boxplot(logcpm, las=2, names = ensembl, col=rep(c("blue", "light blue"), each=3), outline=FALSE , notch = TRUE)
legend("top", legend = c("LUSC", "LUAD"), fill=c("blue", "light blue"), bty = "n", inset = c(-0.05,-0.1), xpd=TRUE, horiz=T)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

plotMDS(logcpm, labels=group, col=rep(c("blue", "light blue"), each=3))
plotMDS(y, col=rep(c("blue", "light blue"), each=3))   

library(statmod)
y <- estimateDisp(y, design)
plotBCV(y, main="BCV Plot")

y$common.dispersion

fit <- glmFit(y, design)

qlf.1vs2 <- glmLRT(fit, contrast=c(-1,1))
results <- as.data.frame(topTags(qlf.1vs2,n = Inf))
library("org.Hs.eg.db")
results$symbol <- mapIds(org.Hs.eg.db,keys=rownames(results), keytype="ENSEMBL", column="SYMBOL", multiVals="first") 
plotMD(qlf.1vs2)
#text(results$logCPM,results$logFC,labels = rownames(results),col="black",cex=0.5,pos=3)
summary(decideTests(qlf.1vs2, p.value=0.01, lfc=2))
qlf.1vs2$table$symbol <- mapIds(org.Hs.eg.db,keys=rownames(qlf.1vs2$table), keytype="ENSEMBL", column="SYMBOL", multiVals="first") 
deg.1vs2 <- topTags(qlf.1vs2, n=Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)$table
up.genes.1vs2 <- deg.1vs2[deg.1vs2$logFC > 0,]
down.genes.1vs2 <- deg.1vs2[deg.1vs2$logFC < 0,]

#dummy volcano plot:
volcanoData <- cbind(results$logFC, -log10(results$FDR))
colnames(volcanoData) <- c("log2FC", "-log10(FDR)")
plot(volcanoData, pch=19)
abline(v=0, col="black", lty=3, lwd=1.0)

res <- cbind(results$symbol,results$logFC,results$PValue,results$FDR)
colnames(res) <- c("Gene", "log2FoldChange", "pvalue", "padj")
write.table(res, "res.txt", append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
res <- read.table("res.txt", header=TRUE)

library(EnhancedVolcano)
library(ggplot2)

EnhancedVolcano(res,
                subtitle = "Volcano plot",
                lab = res$Gene,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'LUSC versus LUAD',
              selectLab = c('FLNC','CALML3','CTSE','SLC30A8','SLC38A11','SLC47A2','SLCO1A2'),
              xlab = bquote(~Log[2]~ 'Fold Change'),
              ylab = bquote(-~Log[10]~ 'FDR'),
               pCutoff = 10e-2,
               FCcutoff = 2.5,
              pointSize = 2.0,
              labSize = 6.0,
              legendLabels=c('NS','Log2 (FC)','FDR','FDR & Log2 (FC)'),
              legendPosition = 'right',
              legendLabSize = 12,
              legendIconSize = 4.0,
              drawConnectors = TRUE,
              widthConnectors = 0.75)

#Results Median Expression Ratio of log2 normalized counts + 1:
lusc <- read.table(file = 'TCGA-LUSC.htseq_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
ensembl <- gsub("\\..*", "", rownames(lusc))
rownames(lusc) <- ensembl
lusc$symbol <- mapIds(org.Hs.eg.db,keys=rownames(lusc), keytype="ENSEMBL", column="SYMBOL", multiVals="first") 
rownames(lusc) <- make.names(lusc$symbol, unique=TRUE)
lusc$symbol <- NULL

luad <- read.table(file = 'TCGA-LUAD.htseq_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
ensembl <- gsub("\\..*", "", rownames(luad))
rownames(luad) <- ensembl
luad$symbol <- mapIds(org.Hs.eg.db,keys=rownames(luad), keytype="ENSEMBL", column="SYMBOL", multiVals="first") 
rownames(luad) <- make.names(luad$symbol, unique=TRUE)
luad$symbol <- NULL

a <- c(median(as.numeric(lusc["FLNC",])),
median(as.numeric(lusc["CALML3",])),
median(as.numeric(lusc["CTSE",])),
median(as.numeric(lusc["SLC30A8",])),
median(as.numeric(lusc["SLC38A11",])),
median(as.numeric(lusc["SLC47A2",])),
median(as.numeric(lusc["SLCO1A2",])))

b <- c(median(as.numeric(luad["FLNC",])),
median(as.numeric(luad["CALML3",])),
median(as.numeric(luad["CTSE",])),
median(as.numeric(luad["SLC30A8",])),
median(as.numeric(luad["SLC38A11",])),
median(as.numeric(luad["SLC47A2",])),
median(as.numeric(luad["SLCO1A2",])))

a - b 

lmts <- range(as.numeric(lusc["SLC47A2",]),as.numeric(luad["SLC47A2",]))
par(mfrow = c(1, 2))
boxplot(as.numeric(lusc["SLC47A2",]), las=2, names = c("LUSC FLNC"), col=rep(c("blue")), notch = TRUE, ylim=lmts)
boxplot(as.numeric(luad["SLC47A2",]), las=2, names = c("LUAD FLNC"), col=rep(c("light blue")), notch = TRUE, ylim=lmts)
legend("top", legend = c("LUSC", "LUAD"), fill=c("blue", "light blue"), bty = "n", inset = c(-0.05,-0.1), xpd=TRUE, horiz=T)

wilcox.test(as.numeric(lusc["SLCO1A2",]), as.numeric(luad["SLCO1A2",]), alternative = "two.sided")

lmts <- range(as.numeric(lusc["CTSE",]),as.numeric(luad["CTSE",]))
par(mfrow = c(1, 2)) 
boxplot(as.numeric(lusc["CTSE",]), las=2, names = c("LUSC FLNC"), col=rep(c("blue")), notch = TRUE, ylim=lmts)
boxplot(as.numeric(luad["CTSE",]), las=2, names = c("LUAD FLNC"), col=rep(c("light blue")), notch = TRUE, ylim=lmts)
legend("top", legend = c("LUSC", "LUAD"), fill=c("blue", "light blue"), bty = "n", inset = c(-0.05,-0.1), xpd=TRUE, horiz=T)

wilcox.test(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",]), alternative = "two.sided")

#Dichotomization based on median:
luscSLC47A2hn <- ifelse(as.numeric(lusc["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))), "High", "Normal")
table(luscSLC47A2hn)
luadSLC47A2hn <- ifelse(as.numeric(luad["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))), "High", "Normal")
table(luadSLC47A2hn)

luscCTSEhn <- ifelse(as.numeric(lusc["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",]))), "High", "Normal")
table(luscCTSEhn)
luadCTSEhn <- ifelse(as.numeric(luad["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",]))), "High", "Normal")
table(luadCTSEhn)

observed <- matrix(c(367, 183, 263, 287), nrow=2, ncol=2, byrow = T)
colnames(observed) <- c("High", "Normal")
rownames(observed) <- c("luscSLC47A2", "luscCTSE")
barplot(observed, beside=T, legend=T, args.legend = list(x = "top"))
chi_test <- chisq.test(observed)
chi_test

observed <- matrix(c(123, 462, 416, 169), nrow=2, ncol=2, byrow = T)
colnames(observed) <- c("High", "Normal")
rownames(observed) <- c("luadSLC47A2", "luadCTSE")
barplot(observed, beside=T, legend=T, args.legend = list(x = "top"))
chi_test <- chisq.test(observed)
chi_test

# SLC47A2 normal CTSE normal
luscnn <- as.numeric(lusc["SLC47A2",]) <= mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(lusc["CTSE",]) <= mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luscnn == TRUE))

luadnn <- as.numeric(luad["SLC47A2",]) <= mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(luad["CTSE",]) <= mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luadnn == TRUE))

# SLC47A2 high CTSE normal
luschn <- as.numeric(lusc["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(lusc["CTSE",]) <= mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luschn == TRUE))

luadhn <- as.numeric(luad["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(luad["CTSE",]) <= mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luadhn == TRUE))

# SLC47A2 normal CTSE high
luscnh <- as.numeric(lusc["SLC47A2",]) <= mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(lusc["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luscnh == TRUE))

luadnh <- as.numeric(luad["SLC47A2",]) <= mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(luad["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luadnh == TRUE))

# SLC47A2 high CTSE high
luschh <- as.numeric(lusc["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(lusc["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luschh == TRUE))

luadhh <- as.numeric(luad["SLC47A2",]) > mean(c(as.numeric(lusc["SLC47A2",]), as.numeric(luad["SLC47A2",]))) & as.numeric(luad["CTSE",]) > mean(c(as.numeric(lusc["CTSE",]), as.numeric(luad["CTSE",])))
length(which(luadhh == TRUE))

# Kaplan meier & Log rank test by tumor type
library(survival)
library(tidyverse)
library(survminer)

survivalLUSC <- read.table(file = 'TCGA-LUSC.survival.tsv', sep = '\t', header = TRUE)
survivalLUAD <- read.table(file = 'TCGA-LUAD.survival.tsv', sep = '\t', header = TRUE)
survivalLUSC['type'] = 'LUSC'
survivalLUAD['type'] = 'LUAD'
write.table(survivalLUSC, "survival.txt", append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(survivalLUAD, "survival.txt", append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)
survival <- read.table("survival.txt", header=TRUE, row.names = 1)

fit_KM <- survfit(Surv(OS.time, OS) ~ type, data = survival)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ type, data = survival)

hazard_ratio <- (fit_log_rank$obs[1] / fit_log_rank$exp[1]) / (fit_log_rank$obs[2] / fit_log_rank$exp[2])
hazard_ratio

# HR = 0.8313035 < 1 indicating that the risk of deaths for LUAD patients is 0.8313035 times the risk of LUSC ones.

# Kaplan meier & Log rank test by sex
library("readxl")
colnames(survivalLUSC)[1] <- "submitter_id.samples"
phenotypeLUSC <- read_excel("TCGA-LUSC.GDC_phenotype.xlsx")
phenotypeLUSC <- merge(survivalLUSC, phenotypeLUSC, by = "submitter_id.samples")
colnames(phenotypeLUSC)[which(names(phenotypeLUSC) == "gender.demographic")] <- "sex"

fit_KM <- survfit(Surv(OS.time, OS) ~ sex, data = phenotypeLUSC)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ sex, data = phenotypeLUSC)
# p= 0.1

colnames(survivalLUAD)[1] <- "submitter_id.samples"
phenotypeLUAD <- read_excel("TCGA-LUAD.GDC_phenotype.xlsx")
phenotypeLUAD <- merge(survivalLUAD, phenotypeLUAD, by = "submitter_id.samples")
colnames(phenotypeLUAD)[which(names(phenotypeLUAD) == "gender.demographic")] <- "sex"

fit_KM <- survfit(Surv(OS.time, OS) ~ sex, data = phenotypeLUAD)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ sex, data = phenotypeLUAD)
# p= 0.8 

# Kaplan meier & Log rank test by age_at_initial_diagnosis
age_at_initial_diagnosis <- ifelse(as.integer(phenotypeLUSC$age_at_initial_pathologic_diagnosis) < 71, 
                                   ifelse(as.integer(phenotypeLUSC$age_at_initial_pathologic_diagnosis) < 56, "39-55", "56-70"),
                                   ifelse(as.integer(phenotypeLUSC$age_at_initial_pathologic_diagnosis) < 81, "71-80", "81-85")
                                  )
fit_KM <- survfit(Surv(OS.time, OS) ~ age_at_initial_diagnosis, data = phenotypeLUSC)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ age_at_initial_diagnosis, data = phenotypeLUSC)
# Chisq= 10.6  on 3 degrees of freedom, p= 0.01 

age_at_initial_diagnosis <- ifelse(as.integer(phenotypeLUAD$age_at_initial_pathologic_diagnosis) < 71, 
                                   ifelse(as.integer(phenotypeLUAD$age_at_initial_pathologic_diagnosis) < 56, "33-55", "56-70"),
                                   ifelse(as.integer(phenotypeLUAD$age_at_initial_pathologic_diagnosis) < 81, "71-80", "81-88")
)
fit_KM <- survfit(Surv(OS.time, OS) ~ age_at_initial_diagnosis, data = phenotypeLUAD)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ age_at_initial_diagnosis, data = phenotypeLUAD)
# Chisq= 11.7  on 3 degrees of freedom, p= 0.009

# Kaplan meier & Log rank test by new_tumor_event_after_initial_treatment
fit_KM <- survfit(Surv(OS.time, OS) ~ new_tumor_event_after_initial_treatment, data = phenotypeLUSC)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ new_tumor_event_after_initial_treatment, data = phenotypeLUSC)
# Chisq= 46.3  on 1 degrees of freedom, p= 1e-11 

hazard_ratio <- (fit_log_rank$obs[1] / fit_log_rank$exp[1]) / (fit_log_rank$obs[2] / fit_log_rank$exp[2])
hazard_ratio

# HR = 0.4339004 < 1 indicating that the risk of deaths in patients with no new tumor a.i.t. is ≈ 0.43 times the risk of patients with a new tumor a.i.t..  

fit_KM <- survfit(Surv(OS.time, OS) ~ new_tumor_event_after_initial_treatment, data = phenotypeLUAD)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ new_tumor_event_after_initial_treatment, data = phenotypeLUAD)
# Chisq= 71.8  on 1 degrees of freedom, p= <2e-16 

hazard_ratio <- (fit_log_rank$obs[1] / fit_log_rank$exp[1]) / (fit_log_rank$obs[2] / fit_log_rank$exp[2])
hazard_ratio

# HR = 0.3156665 < 1 indicating that the risk of deaths in patients with no new tumor a.i.t. is ≈ 0.43 times the risk of patients with a new tumor a.i.t..  


# Kaplan meier & Log rank test by SLC
slcDE <- as.data.frame(luscSLC47A2hn)
slcDE$submitter_id.samples <- str_replace_all(colnames(lusc["SLC47A2",]), '[.]', '-')

survivalLUSC <- merge(survivalLUSC,slcDE, by="submitter_id.samples")
colnames(survivalLUSC)[which(names(survivalLUSC) == "luscSLC47A2hn")] <- "SLC47A2"

fit_KM <- survfit(Surv(OS.time, OS) ~ SLC47A2, data = survivalLUSC)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ SLC47A2, data = survivalLUSC)
# p= 0.6

slcDE <- as.data.frame(luadSLC47A2hn)
slcDE$submitter_id.samples <- str_replace_all(colnames(luad["SLC47A2",]), '[.]', '-')

survivalLUAD <- merge(survivalLUAD,slcDE, by="submitter_id.samples")
colnames(survivalLUAD)[which(names(survivalLUAD) == "luadSLC47A2hn")] <- "SLC47A2"

fit_KM <- survfit(Surv(OS.time, OS) ~ SLC47A2, data = survivalLUAD)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ SLC47A2, data = survivalLUAD)
# p = 0.08 

hazard_ratio <- (fit_log_rank$obs[1] / fit_log_rank$exp[1]) / (fit_log_rank$obs[2] / fit_log_rank$exp[2])
hazard_ratio

# HR = 0.7351113 < 1 indicating that the risk of deaths for LUAD patients is 0.8313035 times the risk of LUSC ones.

# Kaplan meier & Log rank test by CTSE
slcDE <- as.data.frame(luscCTSEhn)
slcDE$submitter_id.samples <- str_replace_all(colnames(lusc["CTSE",]), '[.]', '-')

survivalLUSC <- merge(survivalLUSC,slcDE, by="submitter_id.samples")
colnames(survivalLUSC)[which(names(survivalLUSC) == "luscCTSEhn")] <- "CTSE"

fit_KM <- survfit(Surv(OS.time, OS) ~ CTSE, data = survivalLUSC)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ CTSE, data = survivalLUSC)
# Chisq= 4.4  on 1 degrees of freedom, p= 0.04 

hazard_ratio <- (fit_log_rank$obs[1] / fit_log_rank$exp[1]) / (fit_log_rank$obs[2] / fit_log_rank$exp[2])
# 1.30785

slcDE <- as.data.frame(luadCTSEhn)
slcDE$submitter_id.samples <- str_replace_all(colnames(luad["CTSE",]), '[.]', '-')

survivalLUAD <- merge(survivalLUAD,slcDE, by="submitter_id.samples")
colnames(survivalLUAD)[which(names(survivalLUAD) == "luadCTSEhn")] <- "CTSE"

fit_KM <- survfit(Surv(OS.time, OS) ~ CTSE, data = survivalLUAD)
ggsurvplot(fit_KM,surv.median.line = "hv",risk.table = T, ggtheme = theme_minimal())

fit_log_rank <- survdiff(Surv(OS.time, OS) ~ CTSE, data = survivalLUAD)
# p = 0.4 

# Pairplot
phenotypeLUSCt <- as.data.frame(phenotypeLUSC$age_at_initial_pathologic_diagnosis, phenotypeLUSC$submitter_id.samples)
colnames(phenotypeLUSCt)[which(names(phenotypeLUSCt) == "phenotypeLUSC$age_at_initial_pathologic_diagnosis")] <- "age"
phenotypeLUSCt$new_tumor_a.i.t. <- phenotypeLUSC$new_tumor_event_after_initial_treatment
phenotypeLUSCt$submitter_id.samples <- rownames(phenotypeLUSCt)
survivalLUSC <- merge(survivalLUSC,phenotypeLUSCt, by="submitter_id.samples")
RSEM <- as.data.frame(t(lusc["CTSE",]))
RSEM$submitter_id.samples <- str_replace_all(rownames(RSEM), '[.]', '-')
survivalLUSC$CTSE <- NULL
survivalLUSC <- merge(survivalLUSC,RSEM, by="submitter_id.samples")
RSEM <- as.data.frame(t(lusc["SLC47A2",]))
RSEM$submitter_id.samples <- str_replace_all(rownames(RSEM), '[.]', '-')
survivalLUSC$SLC47A2 <- NULL
survivalLUSC <- merge(survivalLUSC,RSEM, by="submitter_id.samples")

survivalLUSC <- na.omit(survivalLUSC)

library(GGally)
ggpairs(survivalLUSC, c("CTSE", "SLC47A2", "age", "new_tumor_a.i.t."))
ggpairs(survivalLUSC, c("CTSE", "SLC47A2", "age", "new_tumor_a.i.t."), upper = list(continuous = "density", combo = "box_no_facet"), lower = list(continuous = "points", combo = "dot_no_facet"))

# Cox model
survivalLUSC$submitter_id.samples <- NULL
survivalLUSC$X_PATIENT <- NULL
survivalLUSC$type <- NULL
survivalLUSC$SLC47A2 <- NULL

fit_LUSCCTSE <- coxph(Surv(OS.time, OS) ~ ., data=survivalLUSC)
summary(fit_LUSCCTSE)
cox.zph(fit_LUSCCTSE)
ggcoxzph(cox.zph(fit_LUSCCTSE))

# Logistic regression
library(bestglm)
library(pROC)
library(mlbench)
phenotypeLUSC <- read_excel("TCGA-LUSC.GDC_phenotype.xlsx")
RSEM <- as.data.frame(t(lusc["CTSE",]))
RSEM$submitter_id.samples <- str_replace_all(rownames(RSEM), '[.]', '-')
logisticLUSC <- merge(phenotypeLUSC,RSEM, by="submitter_id.samples")
RSEM <- as.data.frame(t(lusc["SLC47A2",]))
RSEM$submitter_id.samples <- str_replace_all(rownames(RSEM), '[.]', '-')
logisticLUSC <- merge(logisticLUSC,RSEM, by="submitter_id.samples")
logisticLUSC$dead <- ifelse(logisticLUSC$vital_status.demographic=="Dead",1,0)
logisticLUSC$stage <- ifelse(logisticLUSC$tumor_stage.diagnoses=="stage ii" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iia" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iib" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iii" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iiia" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iiib" |
                               logisticLUSC$tumor_stage.diagnoses=="stage iv",1,0)
# Other features/columns containing clinical variables have too many NA values

logisticLUSC_full <- as.data.frame(logisticLUSC$stage)
colnames(logisticLUSC_full)[which(names(logisticLUSC_full) == "logisticLUSC$stage")] <- "stage"
logisticLUSC_full$CTSE <- logisticLUSC$CTSE
logisticLUSC_full$SLC47A2 <- logisticLUSC$SLC47A2
logisticLUSC_full$age <- logisticLUSC$age_at_initial_pathologic_diagnosis
logisticLUSC_full$dead <- logisticLUSC$dead
logisticLUSC_full$new_tumor <- logisticLUSC$new_tumor_event_after_initial_treatment
logisticLUSC_full$tissue_source_site <- logisticLUSC$name.tissue_source_site
logisticLUSC_full$stage_complete <- as.factor(logisticLUSC$tumor_stage.diagnoses)
logisticLUSC_full <- na.omit(logisticLUSC_full)

# dead ~ CTSE + age_at_initial_pathologic_diagnosis + new_tumor_event_after_initial_treatment
fit_glm <- glm(dead ~ CTSE + age + new_tumor, data = logisticLUSC_full, family = binomial)
summary(fit_glm)

ROC_curve <- roc(response = logisticLUSC_full$dead, predictor = fit_glm$fitted.values,
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")

y <- logisticLUSC_full$dead
y_hat <- ifelse(fit_glm$fitted.values <= 0.5, 0, 1)
confmat <- table(y, y_hat)
sum(diag(confmat))/sum(confmat)

# AS EXPECTED (confirming our survival analysis' results) THE MODEL DOES NOT IMPROVE ADDING SLC47A2:
fit_glm <- glm(dead ~ CTSE + SLC47A2 + age + new_tumor, data = logisticLUSC_full, family = binomial)
summary(fit_glm)

ROC_curve <- roc(response = logisticLUSC_full$dead, predictor = fit_glm$fitted.values,
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")

y <- logisticLUSC_full$dead
y_hat <- ifelse(fit_glm$fitted.values <= 0.5, 0, 1)
confmat <- table(y, y_hat)
sum(diag(confmat))/sum(confmat)

# stage ~ CTSE + age
fit_glm <- glm(stage ~ CTSE + age + new_tumor, data = logisticLUSC_full, family = binomial)
summary(fit_glm)

ROC_curve <- roc(response = logisticLUSC_full$stage, predictor = fit_glm$fitted.values,
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")

y <- logisticLUSC_full$stage
y_hat <- ifelse(fit_glm$fitted.values <= 0.5, 0, 1)
confmat <- table(y, y_hat)
sum(diag(confmat))/sum(confmat)

# stage ~ CTSE + SLC47A2 + age --> interestingly SLC47A2 is useful for predicting stage!
fit_glm <- glm(stage ~ CTSE + SLC47A2 + age + new_tumor, data = logisticLUSC_full, family = binomial)
summary(fit_glm)

ROC_curve <- roc(response = logisticLUSC_full$stage, predictor = fit_glm$fitted.values,
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")

y <- logisticLUSC_full$stage
y_hat <- ifelse(fit_glm$fitted.values <= 0.5, 0, 1)
confmat <- table(y, y_hat)
sum(diag(confmat))/sum(confmat)

# Logistic mixed effects models with random intercept
library(lme4)
library(insight)

# death
fit_glm <- glmer(dead ~ CTSE + new_tumor + (1 | tissue_source_site), data = logisticLUSC_full, family = binomial)
summary(fit_glm)

# Fixed effects
fixef(fit_glm)

alpha <- 0.05
se <- sqrt(diag(vcov(fit_glm))) #standard errors
# table of estimates with 95% CI using errors we obtained above
CI_betas <- cbind(Est = fixef(fit_glm), 
                  Lower = fixef(fit_glm) - qnorm(1-alpha/2) * se, 
                  Upper = fixef(fit_glm) + qnorm(1-alpha/2) * se)
round(CI_betas,3)

CI_OR <- exp(CI_betas)
round(CI_OR,3)

# Having a new tumor a i t leads to a 7.738 times the risk of death

# Random effect
ranef(fit_glm)

dotplot(ranef(fit_glm))

# Variance Partitioning Coefficient
sigma2_lat <- pi^2/3
sigma2_lat

# Variances of the school intercept
print(VarCorr(fit_glm), comp = c("Variance", "Std.Dev."))

sigma2_b<- as.numeric(get_variance_intercept(fit_glm))

VPC <- sigma2_b/(sigma2_b+sigma2_lat)
VPC
#This VPC value means that 33% of variation in the response is attributed to the classification by tissue source site

logit_p.hat <- predict(fit_glm, logisticLUSC_full)
gl <- binomial(link=logit)
p.hat<-gl$linkinv(logit_p.hat)
p_threshold = sum(logisticLUSC_full$dead)/dim(logisticLUSC_full)[1]
Y.hat <- ifelse(p.hat<p_threshold, 0, 1) 
confusion.matrix <- table(Predicted = Y.hat, Observed = logisticLUSC_full$dead)
confusion.matrix 

# stage ~ new_tumor
fit_glm <- glmer(stage ~ CTSE + SLC47A2 + age + new_tumor + (1 | tissue_source_site), data = logisticLUSC_full, family = binomial)
summary(fit_glm)

# Fixed effects
fixef(fit_glm)

se <- sqrt(diag(vcov(fit_glm))) #standard errors
# table of estimates with 95% CI using errors we obtained above
CI_betas <- cbind(Est = fixef(fit_glm), 
                  Lower = fixef(fit_glm) - qnorm(1-alpha/2) * se, 
                  Upper = fixef(fit_glm) + qnorm(1-alpha/2) * se)
round(CI_betas,3)

CI_OR <- exp(CI_betas)
round(CI_OR,3)

# Having a new tumor a i t leads to a 7.738 times the risk of death

# Random effect
ranef(fit_glm)

dotplot(ranef(fit_glm))

# Variances of the school intercept
print(VarCorr(fit_glm), comp = c("Variance", "Std.Dev."))

sigma2_b<- as.numeric(get_variance_intercept(fit_glm))

VPC <- sigma2_b/(sigma2_b+sigma2_lat)
VPC
#This VPC value means that 1.4% of variation in the response is attributed to the classification by tissue source site

logit_p.hat <- predict(fit_glm, logisticLUSC_full)
gl <- binomial(link=logit)
p.hat<-gl$linkinv(logit_p.hat)
p_threshold = sum(logisticLUSC_full$stage)/dim(logisticLUSC_full)[1]
Y.hat <- ifelse(p.hat<p_threshold, 0, 1) 
confusion.matrix <- table(Predicted = Y.hat, Observed = logisticLUSC_full$stage)
confusion.matrix 

### violin plots ###

library(vioplot)

#LUSC vs LUAD for CALML3
#I tried to add the mean point but I don't know why the points(...) do not show anything
par(mfrow = c(1,2))

vioplot(as.numeric(lusc["CALML3",]), horizontal = T, col = "blue", rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red")
points(mean(as.numeric(lusc["CALML3",])), pch = 19, col = "black", cex = 2)

vioplot(as.numeric(luad["CALML3",]), horizontal = T, col = "lightblue", rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red")
points(mean(as.numeric(luad["CALML3",])), pch = 19, col = "black", cex = 2)

legend("bottom", legend = c("LUAD", "LUSC", "median", "mean"), fill = c("blue", "lightblue", 0, 0), border = NA, bty = "n", inset = c(-0.3,-0.3), xpd = T, horiz = T, pch = c(NA, NA, 16, 19), col = c("red","black"))

par(mfrow = c(1,1))

#or, to make it  a little prettier;
par(mfrow = c(1,2))

vioplot(as.numeric(lusc["CALML3",]), horizontal = T, col = "dodgerblue", rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red", main = "LUAD")
points(mean(as.numeric(lusc["CALML3",])), pch = 23, col = "black", cex = 2)
legend("bottomright", legend = "median",border = NA, bty = "n", inset = c(-0.2,-0.3), xpd = T, horiz = T, pch = 16, col = "red")

vioplot(as.numeric(luad["CALML3",]), horizontal = T, col = "lightblue1", rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red", main = "LUSC")
points(mean(as.numeric(luad["CALML3",])), pch = 23, col = "black", cex = 2)
legend("bottomleft", legend = "mean", border = NA, bty = "n", inset = c(-0.3,-0.3), xpd = T, horiz = T, pch = 19, col = "black")

par(mfrow = c(1,1))

# we could also compare the profile of different genes (if we choose to study more than one) I randomly chose these two just to give an idea
vio_luad_df <- data.frame(as.numeric(luad["CTSE",]), as.numeric(luad["SLC7A2",]))
vio_lusc_df <- data.frame(as.numeric(lusc["CTSE",]), as.numeric(lusc["SLC7A2",]))

par(mfrow = c(1,2))

vioplot( vio_luad_df, col = c("dodgerblue", "dodgerblue4"), names=c("CTSE", "SLC7A2"), , rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red", main = "LUAD genes CTSE vs SLC7A2")
vioplot( vio_lusc_df, col = c("lightblue1", "lightblue4"), names=c("CTSE", "SLC7A2"), , rectCol = "white", lineCol = "white", pchMed = 16, colMed = "red", main = "LUSC genes CTSE vs SLC7A2")
legend("bottomleft", legend = "mean", border = NA, bty = "n", inset = c(-0.3,-0.3), xpd = T, horiz = T, pch = 19, col = "black")

par(mfrow = c(1,1))

#META-ANALYSIS
library(metafor)
library(metadat)
HR <- c(1.30785, 1.41649)
yi = log(HR)
vi<-c(0.01,0.009)
dat <- data.frame(yi=yi, sei=sei)
trial <- 1:2 
author <- c("Our project", "Yang Y, Wang M, Liu B")
year <- c(2022, 2018)
dat <- cbind.data.frame(trial=trial, author=author, year=year, dat)
fit_FEM <- rma(yi, vi=vi, data=dat, method="FE")
fit_FEM
forest(fit_FEM, slab = paste(dat$author, dat$year, sep = ", "), xlim = c(-3, 3))
text(-3, 10, "Study", pos = 4)
text(3, 10, "Standardized Mean Difference [95% CI]", pos = 2)
