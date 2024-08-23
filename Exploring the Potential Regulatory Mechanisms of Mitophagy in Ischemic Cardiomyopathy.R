#####SJZZK-451-3
rm(list=ls())

setwd('/data/nas1/dailihui/re-project/01.SJZZK-451-3/')

mycol <- c('#02b1e6', '#E81D22','#F9BC15','#20e81d')
##
####### GEO data#########
library(tidyr)
library(tibble)
library(dplyr)
library(data.table)
library(GEOquery)
library(org.Hs.eg.db)
library(clusterProfiler)

matrix<-read.delim2('GSE5406_series_matrix.txt',sep = '\t',check.names = F,header = T)
anno <-data.table::fread("GSE5406_family.soft",skip ="ID",header = T)
# write.csv(anno,"anno.csv")
View(anno)
anno$`Gene Symbol`

#filter repeat value
probe2symbol <- anno %>%
  dplyr::select("ID","Gene Symbol") %>% dplyr::rename(probeset = "ID",symbol="Gene Symbol") %>%
  filter(symbol!= "") %>%
  tidyr::separate_rows( `symbol`,sep="///")


rownames(matrix)=matrix[,1]
matrix=matrix[,-1]
exprSet <- matrix %>% as.data.frame() %>%
  rownames_to_column(var="probeset") %>%
  inner_join(probe2symbol,by="probeset") %>%
  dplyr::select(-probeset) %>%
  dplyr::select(symbol,everything())%>%
  filter(symbol != "NA") %>%
  distinct(symbol,.keep_all = T) %>%
  dplyr::select(-rowMean) %>%
  column_to_rownames(var = "symbol")


# exprSet <- read.csv("GSE36961_expr.csv",check.names = F)
duplicated(exprSet$symbol)
exprSet1 = aggregate(.~exprSet$symbol,max,data = exprSet)

# exprSet1 <- read.csv('GSE165004-expr.csv',check.names = F)
b <- read.csv("D:\\human_mRNAs.csv",check.names = F)
out = exprSet1[exprSet1$symbol %in% b$Genename,]
rownames(out) = out[,1]
out = out[,-c(1,2)]
write.csv(out,"GSE5406_mRNA.csv")#output



###########1.Mitochondrial autophagy score calculation#####

# dir.create('1.mt')

library(GSVA)

matrix  <- read.csv('GSE116250-mRNA.csv',check.names = F,row.names = 1)
b <- read.csv('GSE116250-sample.csv',check.names = F)
id = b$Title
matrix = matrix[,id]
write.csv(matrix,'GSE116250-mRNA-ICM.csv')
geneSet <- read.delim2('PMID35198564-MRGs.txt',sep = '\t',check.names = F,header = T)

exprSet <- as.matrix(matrix)

gsva_matrix <- gsva(exprSet,geneSet,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix = t(gsva_matrix)
write.csv(gsva_matrix,'1.mt//Mitophagy.csv')

##boxplot
dbox1 <- read.csv('1.mt/Mitophagy.csv',check.names = F,row.names = 1)
library(ggpubr)
dbox1$group <- c(rep('ICM',13),rep('NC',14))
table(dbox1$group)
my_comparition <- list(c("ICM","NC"))


pdf('1.vlin.pdf',width = 4.5,height = 4.5)
p <- dbox1 %>%
  ggviolin(x= "group",y = colnames(dbox1)[1],
           fill = "group",combine = T,
           palette = mycol,
           ylab = "",xlab = "",
           add = "boxplot",add.params = list(fill="white"))

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition)
dev.off()

###############1.WGCNA#######


library(WGCNA)
library(tidyr)
enableWGCNAThreads() #multi-core

#BiocManager::install("matrixStats")
#load("dataExpr.Rdata")
dataExpr <- read.csv("GSE116250-mRNA-ICM.csv", check.names = F, row.names = 1) %>%
  t() %>%
  as.data.frame()
dataExpr <- apply(dataExpr, 2, function(x) { log2(x + 1)})
sample <- read.csv("Mitophagy评分.csv", row.names = 1, check.names = F)
id <- rownames(sample)
head(sample)
dataExpr <- dataExpr[id, ]


# t(expression)，save dataExpr.Rdata
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste(
      "Removing genes:",
      paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")
    ))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste(
      "Removing samples:",
      paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")
    ))
  }
  # Remove the offending genes and samples from the data:
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes <- ncol(dataExpr)
nSamples <- nrow(dataExpr)
dim(dataExpr)

# fliter again
dataExpr <- dataExpr[, colMeans(dataExpr) > 1]
save(dataExpr, file = "dataExpr.Rdata")
sampleTree <- hclust(dist(dataExpr), method = "ward.D2")
pdf(file = "Fig1.Sample_cluster.pdf", width = 10, height = 5.5)
par(mar = c(0, 4, 2, 0), xpd = F)
plot(sampleTree, main = "Sample clustering to detect outliers", cex.axis = 1, cex.main = 1, sub = "", xlab = "")
abline(h = 230, col = "red")
dev.off()

############## 
### Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 230, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
# dataExpr = dataExpr[keepSamples,]

# Taking disease and health as traits，number from 0-9
dataTraits <- sample
dataTraits <- dataTraits[rownames(dataExpr), ]
# rownames(dataTraits) <- rownames(dataExpr)
head(dataTraits)
save(dataTraits, file = "dataTraits.Rdata")
sampleTree2 <- hclust(dist(dataExpr), method = "ward.D2")
traitColors <- numbers2colors(dataTraits, signed = FALSE)

pdf("Fig2.cluster.dendrogram.pdf", width = 15, height = 7)
# par(pin=c(14, 8), mai=c(10,500,5,5), mar = c(10,500,5,5))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(dataTraits), cex.dendroLabels = 1,
                    marAll = c(1, 7, 3, 1),
                    # cex.colorLabels = 0.6,
                    main = "Sample dendrogram and trait heatmap"
)
dev.off()
# save(dataExpr,file='dataExpr.Rdata')
## power 
powers <- seq(from = 1, to = 20, by = 1) # Power exponent range 1:20
sft <- pickSoftThreshold(dataExpr,
                         RsquaredCut = 0.85, # networkType = "signed",
                         powerVector = powers, verbose = 5
)
sft$powerEstimate
# save(sft,file='sft.Rdata')
pdf("Fig3.Soft.Threshold.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
cex1 <- 0.85
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.857, col = "red") 



plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()


softPower <- sft$powerEstimate # the best power value
#softPower = 9

adjacency <- adjacency(dataExpr, power = softPower)
#softPower
###
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
# save(dissTOM,file='dissTOM.Rdata')
# load('dissTOM.Rdata')
### 
geneTree <- hclust(as.dist(dissTOM), method = "average")
save(geneTree, file = "geneTree.Rdata")
# load('geneTree.Rdata')
pdf(file = "Gene.cluster.pdf", width = 6, height = 5)
plot(geneTree,
     xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04
)
dev.off()


##
minModuleSize <- 100 # module gene number
dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
save(dynamicMods, file = "dynamicMods.Rdata")
# load('dynamicMods.Rdata')
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

pdf(file = "Fig4.DynamicTree.pdf", width = 6, height = 5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors"
)
dev.off()

###############1. WGCNA#########
library(WGCNA)


load("1.mt//dataExpr.Rdata")
load("1.mt//dynamicMods.Rdata")
dynamicColors <- labels2colors(dynamicMods)
MEList <- moduleEigengenes(dataExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# save(MEList, file = "MEList.RDdata")

# moduleColors = MEList$validColors
# table(moduleColors)
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

pdf("1.mt//Fig5.merge.cluster.pdf", width = 7, height = 6)
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "", sub = ""
)
MEDissThres <- 0.2 # change the hight
abline(h = MEDissThres, col = "red")
dev.off()

### rbind the similar module
merge <- mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
load("1.mt/geneTree.Rdata")
pdf(file = "1.mt//Fig6.DynamicTree.pdf", width = 6, height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors"
)
dev.off()
moduleColors <- mergedColors
# save(moduleColors,file = 'moduleColors.Rdata')

table(moduleColors)
colorOrder <- c("grey", standardColors(30))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

### heatmap
nGenes <- ncol(dataExpr)
nSamples <- nrow(dataExpr)
load("1.mt/dataTraits.Rdata")
moduleTraitCor <- cor(MEs, dataTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")",
                    sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
write.table(moduleTraitPvalue, "1.mt/moduleTraitPvalue.txt", sep = "\t", quote = F, row.names = T)
write.table(moduleTraitCor, "1.mt/moduleTraitCor.txt", sep = "\t", quote = F, row.names = T)

pdf(file = "1.mt/Fig7.Module_trait.pdf", width = 9, height = 9)
par(mar = c(10, 10, 3, 3))


labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(dataTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(500),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  main = paste("Module-trait relationships")
)
dev.off()

## cor
dir.create("1.mt/across_modules")
for (i in 1:ncol(dataTraits)) {
  which.trait <- colnames(dataTraits)[i]
  
  moduleTraitCor[, which.trait]
  
  y <- dataTraits[, which.trait]
  
  GS <- as.numeric(cor(y, dataExpr, use = "p"))
  
  GeneSignificance <- abs(GS)
  
  ModuleSignificance <- tapply(
    
    GeneSignificance,
    moduleColors, mean,
    na.rm = T
  )
  pdf(file = paste0("1.mt/across_modules/", which.trait, "-", "Gene significance across modules.pdf"), width = 12, height = 5)
  plotModuleSignificance(GeneSignificance, moduleColors)
  dev.off()
}


### MM GS
# modNames = substring(names(MEs), 3)
# modNames <- c('brown','greenyellow','yellow','cyan','grey60','pink')
modNames <- unique(moduleColors)
geneModuleMembership <- as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) <- paste("p.", colnames(MMPvalue), sep = "")

traitNames <- colnames(dataTraits)
geneTraitSignificance <- as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")



pdf("1.mt/Fig8.Gene.trait_MEturquoise.pdf", width = 7, height = 7)
# par(mfrow=c(6,6))

for (i in colnames(geneTraitSignificance)) {
  moduleGenes <- moduleColors == "turquoise"
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, "MEturquoise"]),
                     abs(geneTraitSignificance[moduleGenes, i]),
                     xlab = ("Module Membership in module"),
                     ylab = "Gene significance",
                     main = paste0("Trait ", substr(i, 4, nchar(i)), "\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = "turquoise"
  )
  # abline(v=0.8,h=0.5,col="yellow")
}
dev.off()

pdf("1.mt/Fig8.Gene.trait_MEblue.pdf", width = 7, height = 7)
# par(mfrow=c(6,6))

for (i in colnames(geneTraitSignificance)) {
  moduleGenes <- moduleColors == "blue"
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, "MEblue"]),
                     abs(geneTraitSignificance[moduleGenes, i]),
                     xlab = ("Module Membership in module"),
                     ylab = "Gene significance",
                     main = paste0("Trait ", substr(i, 4, nchar(i)), "\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = "blue"
  )
  # abline(v=0.8,h=0.5,col="yellow")
}
dev.off()


### output
probes <- colnames(dataExpr)
geneInfo0 <- data.frame(
  probes = probes,
  moduleColor = moduleColors
)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, geneTraitSignificance[, Tra],
    GSPvalue[, Tra]
  )
  names(geneInfo0) <- c(
    oldNames, names(geneTraitSignificance)[Tra],
    names(GSPvalue)[Tra]
  )
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, geneModuleMembership[, mod],
    MMPvalue[, mod]
  )
  names(geneInfo0) <- c(
    oldNames, names(geneModuleMembership)[mod],
    names(MMPvalue)[mod]
  )
}
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder, ]
write.table(geneInfo, file = "1.mt/GS_MM.csv", sep = ",", row.names = F)

####



###############2. GSE116250-DEGS########
library(limma)
# dir.create('2.modulegene')
a <- read.csv('GSE116250-mRNA-ICM.csv',check.names = F,row.names = 1)
b <- read.csv('GSE116250-sample.csv',check.names = F,row.names = 1)

id  = b$Title
table(b$disease)
aa = a[,id]
aa = apply(aa,2,function(x){log2(x+1)})

group_list=c(rep('ICM',13),rep('Control',14))
group_list <- factor(group_list,levels = c("ICM","Control"))

pvalue <-0.05
logFoldChange <- 1
dat <- aa
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design

# DEGs contrast matrix
contrast.matrix <- makeContrasts(ICM-Control, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#Retain 4 digits
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)

write.csv(cbind(Symbol=rownames(allDiff),allDiff),file="2.modulegene//limmaOut.csv")

##DEGs

diffSig = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.csv(diffSig,"2.modulegene/DEGs.csv")

diffUp = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC>logFoldChange)),]
write.csv(diffUp,"2.modulegene/up-DEGs.csv")

diffDown = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.csv(diffDown,"2.modulegene/down-DEGs.csv")



##########TOP10########
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
#test
df <- read.csv('2.modulegene/limmaOut.csv',check.names = F,row.names = 1)

df$Symbol = rownames(df)

colnames(df)
head(df)
##p＜0.05，|log2FC|＞1
pvalue = 0.05
log2FC = 1
#up down type
df$group <- case_when(
  df$logFC > log2FC& df$adj.P.Val < pvalue ~ "up",
  df$logFC < -log2FC & df$adj.P.Val < pvalue ~ "down",
  TRUE ~ 'none'
)
head(df)
table(df$group)

#turn to the factor
df$'-log10(padj)' <- -log10(df$adj.P.Val) 
df$group <- factor(df$group, levels = c("up","down","none"))

p <- ggplot(data = df,
            aes(x = logFC, y = -log10(adj.P.Val), color = group)) + 
  geom_point(size = 2.2) 
p

#coordinate axis
p1 <- p +
  scale_x_continuous(limits = c(-11.5, 11.5), breaks = seq(-12, 12, by = 4))+
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 11), breaks = seq(0, 11, by = 10))
p1

#color
mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
p2 <- p1 +
  scale_colour_manual(name = "GSE116250 DEGs", values = alpha(mycol, 0.7)) +
  mytheme
p2

#Add threshold line：
p3 <- p2 +
  geom_hline(yintercept = c(-log10(pvalue)),size = 0.7,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed") 
p3




DEG<-df
rownames(DEG)
up_DEG <- DEG[which(DEG$group == "up"),]
up_DEG <- up_DEG[order(-up_DEG$logFC , up_DEG$P.Value),]
down_DEG <- DEG[which(DEG$group == "down"),]
down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$P.Value ),]
data_repel <- rbind(up_DEG[1:10,], down_DEG[1:10,])



p5 <- p3 +
  geom_text_repel(data = data_repel,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)
p5
# ggsave('./volcano.pdf',width=7,height =6)
# ggsave('/volcano.png',width=7,height =6)
ggsave('2.modulegene//volcano.pdf',width=7,height =6)


###
############TOP20 heatmap######
library(pheatmap)
library(ComplexHeatmap)


c <- read.csv("2.modulegene//DEGs.csv",check.names = F,row.names=1)
c = c[1:20,]
gene = rownames(c)
expr = dat[gene,]
expr = na.omit(expr)
# expr = apply(expr,2,function(x){log2(x+1)})##TCGA need log,GEO do not log

gene=expr
densityHeatmap(gene)

ha1 = HeatmapAnnotation(group = c(rep("ICM", 13), rep("Control", 14)))
ha2 = HeatmapAnnotation(foo = anno_points(rnorm(27)))

densityHeatmap(gene, top_annotation = ha1, bottom_annotation = ha2)

pdf('2.modulegene/TOP20_heatmap.pdf',width = 9,height = 7)

densityHeatmap(gene, top_annotation = ha1, bottom_annotation = ha2) %v%
  HeatmapAnnotation(foo = anno_barplot(1:27)) %v%
  Heatmap(gene,show_column_names = F,
          
          name = "expression", height = unit(9, "cm"),)

dev.off()
###

################3. module venn#####
library(ggVennDiagram)
emrg <- read.csv('1.mt/modulegene.csv',check.names = F)
emrg <- read.csv('1.mt/GS_MM.csv',check.names = F)
emrg = emrg[emrg$moduleColor == 'turquoise',]

deg <- read.csv('2.modulegene/DEGs.csv',check.names = F,row.names = 1)


pdf("2.venn.pdf", width = 6, height = 6)

venn::venn(list(`EMRGs`= emrg$probes,`DEGs`= rownames(deg)),
           zcolor= c("#3F60AA","#CC340C",'darkgreen'),box=T,sncs = 2,ilcs = 1.5)
dev.off()

x1 = intersect(emrg$probes,rownames(deg))
xx = deg[x1,]
write.csv(xx,'2.modulegene/DE-WGCNA.csv')


x = list(`EMRGs`= emrg$probes,`DEGs`= rownames(deg))

pdf('2.modulegene/venn.pdf',width = 6, height = 4)
ggVennDiagram(x) + scale_fill_gradient(low="#0099b499",high= "#ed000099")
dev.off()
###



################3.gene enrichment analysis#######
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

# dir.create('3.modulefunc')
########GO enrichment################
# transform ID to ENTREZID

deg.aging.related <- read.csv('2.modulegene/DE-WGCNA.csv',check.names = F,row.names = 1)

deg.aging.related = rownames(deg.aging.related)
target_gene_id<- bitr(deg.aging.related,fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
write.csv(target_gene_id,file='gene_entriz.csv')

result_go_aging_genes <- enrichGO(target_gene_id$ENTREZID,
                                  OrgDb = "org.Hs.eg.db",
                                  keyType = "ENTREZID",
                                  pvalueCutoff = 1,qvalueCutoff = 1,
                                  ont = "ALL",readable = T)

write.csv(result_go_aging_genes,"3.modulefunc/enrichGO.csv")


##function
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

##terms
display_number = c(10, 10, 10)
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$ENTREZID,
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
# write.csv(ego_result_MF,"GO_MF.csv")
ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$ENTREZID,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$ENTREZID,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]

go_enrich_df <- data.frame(
  ID = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
  GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type = factor(c(
    rep("biological process", display_number[3]), rep("cellular component", display_number[2]),
    rep("molecular function", display_number[1])
  ), levels = c("molecular function", "cellular component", "biological process"))
)

write.csv(go_enrich_df,"go_enrich_df.csv")
# go_enrich_df <- read.csv("go_enrich_df.csv",check.names = F)

## numbers as data on x axis

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

labels <- (sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names
))

# names(labels) <- rev(1:nrow(go_enrich_df))
labels <- as.factor(rev(go_enrich_df$Description))



## colors for bar // green, blue, orange
CPCOLS <- c('#02b1e6', '#E81D22','#F9BC15')

p <- ggplot(data = go_enrich_df, aes(x = number, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_bw() +
  scale_x_discrete(labels = labels) +
  xlab("") +
  theme(axis.text = element_text(face = "bold", color = "gray10")) +
  labs(title = "The Most Enriched GO Terms")

p

pdf("2.go_barplot.pdf", width = 8, height = 6)
p
dev.off()


kk = result_go_aging_genes
pdf(file="2.GO_barplot_p.pdf",width = 8,height = 7)

barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY",font.size = 12) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=80))+
  scale_colour_gradient(low="blue",high="#cc340c")+
  facet_grid(ONTOLOGY~., scale='free')

dev.off()


pdf(file="2.GO_bubble_p.pdf",width = 8,height = 8)
enrichplot::dotplot(kk,showCategory = 10,split="ONTOLOGY",font.size = 12) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=80))+
  scale_colour_gradient(low="midnightblue",high="#cc340c")+
  facet_grid(ONTOLOGY~.,scale='free')
dev.off()


###P-value color

kegg <- read.csv("enrichGO.csv",row.names = 1)
kegg = kegg[kegg$ONTOLOGY =="BP",]
#order by p.value
kegg <- kegg[order(kegg$pvalue),]
#top15 KEGG pathway
kegg <- kegg[1:15,]
#gene numbers of per-pathway
top10 <- data.frame(kegg$Description,kegg$Count ,kegg$pvalue)
colnames(top10) <- c("Description","count","P-value")
#fill=padj fill
p <- ggplot(data=top10,aes(x=Description,y=count,fill=`P-value`))
#coord_flip()
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
#ylim(0,30) 
p3 <- p2 + ylim(0,4.2) + scale_fill_gradient(low="firebrick1",high="royalblue")
p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="",title="GO BP Terms")

#output
pdf("GO_barplot_BP.pdf",width=8,height = 6)
print(p4)
dev.off()
##################3. GO chord################
library(GOplot)
DEGs <- read.csv('2.modulegene/DE-WGCNA.csv',row.names = 1,check.names = F)

DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]
entrezIDs <- mget(rownames(DEGs), org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(DEGs,entrezID=entrezIDs)
out = out[,c(1,7)]
out <- cbind(symbol = rownames(out),out)
write.table(out,file="GSEA_id.xls",sep="\t",quote=F,row.names = F)


id.fc <- read.table("GSEA_id.xls",sep="\t",header = T)
head(id.fc)
ego <- enrichGO(gene = id.fc$entrezID,
                #mouse
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                #human
                #OrgDb = org.Hs.eg.db,
                #non-model organism
                #OrgDb = maize.db,
                ont = "BP", 
                # pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1) 
dim(ego)






go <- data.frame(Category = "CC",
                 ID = ego$ID,
                 Term = ego$Description, 
                 Genes = gsub("/", ", ", ego$geneID), 
                 adj_pval = ego$pvalue)

id.fc <- read.table("GSEA_id.xls",sep="\t",header = T)
head(id.fc)

genelist <- data.frame(ID = id.fc$entrezID, logFC = id.fc$logFC)             

circ <- circle_dat(go, genelist)
head(circ)

#transform ENTREZ ID to gene symbol
id.gsym <- bitr(circ$genes, 
                fromType = "ENTREZID",
                toType = "SYMBOL",
                OrgDb = "org.Hs.eg.db") 


rownames(id.gsym) <- id.gsym$ENTREZID
circ.gsym <- circ
circ.gsym$genes <- id.gsym[circ$genes,]$SYMBOL
head(circ.gsym)


n = 10 
chord <- chord_dat(circ, genelist, go$Term[1:n])
head(chord)
#transform ENTREZ ID to gene symbol
id.gsym <- bitr(row.names(chord), 
                fromType = "ENTREZID", 
                toType = "SYMBOL", 
                OrgDb = "org.Hs.eg.db") 

rownames(id.gsym) <- id.gsym$ENTREZID
head(id.gsym)

chord.gsym <- chord
row.names(chord.gsym) <- id.gsym[row.names(chord),]$SYMBOL
head(chord.gsym)
circ.gsym=circ.gsym[1:10,]
GOChord(chord.gsym, 
        space = 0.02, 
        gene.order = 'logFC', 
        lfc.col = c('#cc340c', 'white', 'darkblue'), #up down color
        gene.space = 0.25, #distance
        gene.size = 6, #size 
        border.size = 0, 
        process.label = 8,#term size
        #ribbon.col=colorRampPalette(c("blue", "red"))(length(EC$process)),#GO term color
        ribbon.col=brewer.pal(length(circ.gsym$term), "Paired"),#GO term color
) 

ggsave("3.modulefunc/GOChord-BP.pdf", width = 14, height = 14)

###
########KEGG enrichment#################

result_KEGG <- enrichKEGG(target_gene_id$ENTREZID,organism = "hsa",
                          use_internal_data = T,
                          pvalueCutoff =1,qvalueCutoff =1)

write.csv(result_KEGG,"3.modulefunc//enrichKEGG.csv")

###transform entrize to Symbol

enrichOutput<- read.csv("3.modulefunc//enrichKEGG.csv",check.names = F)
biotype <- read.csv("gene_entriz.csv",check.names = F)

enrichOutput$geneSymbol <- unlist(lapply(enrichOutput$geneID, function(v)
  
  paste(biotype$SYMBOL[match(strsplit(as.character(v), '/', fixed=TRUE)[[1]],
                             biotype$ENTREZID)], collapse = '/')))
write.csv(enrichOutput,"3.modulefunc//enrichKEGG.csv")

result_KEGG <- as.data.frame(result_KEGG)
result_KEGG = result_KEGG[1:20,]
result_KEGG$number <- factor(rev(1:nrow(result_KEGG)))

KEGG_enrich_df <- data.frame(
  ID = c(result_KEGG$ID),
  Description = c(result_KEGG$Description),
  GeneNumber = c(result_KEGG$Count),
  type = factor(c(
    rep("KEGG Pathway"))
  ), levels = c("KEGG Pathway"))


labels <- sapply(
  levels(KEGG_enrich_df$Description)[as.numeric(KEGG_enrich_df$Description)],
  shorten_names
)

# names(labels) <- rev(1:nrow(go_enrich_df))
labels <- as.factor(rev(KEGG_enrich_df$Description))
# DATA <- KEGG_enrich_df[1:20,]
## colors for bar // green, blue, orange
CPCOLS <- c("#f2cb48")

p <- ggplot(data = KEGG_enrich_df, aes(x = ID, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.3) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_bw() +
  scale_x_discrete(labels = labels) +
  xlab("") +
  theme(axis.text = element_text(face = "bold", color = "gray10")) +
  labs(title = "The Most Enriched KEGG Pathways")

p

pdf("2. DE-MRGs/KEGG_barplot.pdf", width = 8, height = 6)
p
dev.off()



k = result_KEGG
pdf(file="3.KEGG_barplot_p.pdf",width = 8,height = 6)
# png(file="KEGG_barplot.png",width = 900,height = 600)
barplot(k,showCategory = 12,font.size = 15,fill = pvalue)
dev.off()


pdf("3.KEGG_bubble_p.pdf",width = 8,height = 6)
# png(file="KEGG_bubble.png",width = 900,height = 600)
enrichplot::dotplot(k,showCategory = 12,font.size = 15)+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=80))+ 
  scale_colour_gradient(low="midnightblue",high="#cc340c")
dev.off()



kegg <- read.csv("enrichKEGG.csv",row.names = 1)


kegg <- kegg[order(kegg$pvalue),]
kegg <- kegg[1:6,]
top10 <- data.frame(kegg$Description,kegg$Count ,kegg$pvalue)
colnames(top10) <- c("Description","count","P-value")
p <- ggplot(data=top10,aes(x=Description,y=count,fill=`P-value`))
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
p3 <- p2 + ylim(0,4.2) + scale_fill_gradient(low="#ffcc25",high="royalblue")
p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="",title="KEGG Pathway")

#output
pdf("KEGG_barplot_p.pdf",width=9,height = 7)
print(p4)
dev.off()



#GO dotplot##
goplot <- as.data.frame(GO@result)
goBP <- subset(goplot,subset = (ONTOLOGY == "BP") & goplot$pvalue < 0.05) %>% head(10)
goBP <- goBP[order(goBP$Count, decreasing = F), ]
goCC <- subset(goplot,subset = (ONTOLOGY == "CC") & goplot$pvalue < 0.05) %>% head(10)
goCC <- goCC[order(goCC$Count, decreasing = F), ]
goMF <- subset(goplot,subset = (ONTOLOGY == "MF") & goplot$pvalue < 0.05) %>% head(10)
goMF <- goMF[order(goMF$Count, decreasing = F), ]
go.df <- rbind(goMF,goCC,goBP)
go.df1 <- go.df %>% separate(GeneRatio, into = c("enrichGene","totalGene"), sep = "/")
go.df1$Generatio <- as.numeric(go.df1$enrichGene)/as.numeric(go.df1$totalGene)
##levels ordering factor
x=go.df1$Generatio
y=factor(go.df1$Description,levels = go.df1$Description)
##
pdf(file = '03_DEGEnrich/DEGsEnrichGO.bubble.pdf', width = 10,height = 8)
go_bubb = ggplot(go.df1, aes(x,y))+
  geom_point(aes(size=Count,color=-1*log(p.adjust)))+
  scale_color_gradient(low = "SpringGreen", high = "DeepPink")+
  facet_grid(ONTOLOGY~., scale='free')+
  #facet_wrap(~ ONTOLOGY, nrow = 1)
  #colors = c(rep('#FF6347', 10), rep('#FFA54F', 10), rep('#40E0D0', 10))
  labs(color=expression(-log(p.adjust)),
       size="Gene Number",
       x="EnrichmentScore",
       y="Go terms",
       title="Go enrichment")+
  scale_y_discrete(labels = function(x) str_wrap(x,width = 50))+ 
  theme(axis.title = element_text(size = 13,face = "bold"), 
        axis.text = element_text(size = 12,face = "bold"),
        #axis.text.y = element_text(colour = colors),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 13,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
print(go_bubb)
dev.off()

#KEGG##
KEGGoutT30 <- subset(KEGGout, KEGGout$pvalue < 0.05) %>% head(15)
KEGGoutT30 <- KEGGoutT30[order(KEGGoutT30$Count),]
x=KEGGoutT30$Count
y=factor(KEGGoutT30$Description,levels = KEGGoutT30$Description)
pdf(file = '03_DEGEnrich/DEGsEnrichKEGG.bar.pdf',width = 9,height = 7)
kegg_bar <- ggplot(KEGGoutT30, aes(y=y, x=x, fill=-log(p.adjust)))+
  geom_bar(stat = "identity",width=0.9)+####column width
  scale_y_discrete(labels = function(x) str_wrap(x,width = 50))+
  scale_fill_gradient(low = "#20e81d",high ="#8015f9" )+
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme_bw()+
  theme(axis.title = element_text(size = 13,face = "bold"), 
        axis.text = element_text(size = 12,face = "bold"), 
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 13,face = "bold"), 
        legend.text = element_text(size = 12,face = "bold"), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
print(kegg_bar)  
dev.off()





library(GOplot)
diffSig <- read.csv("2.AML-EMRGs.csv",check.names = F)
kkid_file <- read.csv('3.enrichKEGG.csv',check.names = F)
# kkid_file = kkid_file[1:20,]
kkid <- data.frame(Category = 'All',
                   ID = kkid_file$ID,
                   Term = kkid_file$Description,
                   Genes = gsub('/',',',kkid_file$geneID),
                   adj_pval = kkid_file$pvalue)

allentrez <- mapIds(org.Hs.eg.db,keys = diffSig$x,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='first')

diffSig1 <- data.frame(ID = allentrez, logFC = diffSig$log2FoldChange) 

cirkegg <- circle_dat(kkid,diffSig1) 
table(cirkegg$ID)

pdf('3.KEGG.gossip.pdf',width = 10,height = 7)
GOCircle(cirkegg,rad1=2.5,rad2=3.5,label.size=3.5,nsub=10,zsc.col=c('#E81D22', 'white','#02b1e6'))  #nsub是富集出来的通路数目
dev.off()
###
#############Treemap########
library(treemap)
library(RColorBrewer)
ExampleWOS<-data.frame(group=c("USA","USA","China","France","USA","Ireland","USA","France","UK","USA"),
                       subgroup=c("UNIVERSITY OF CALIFORNIA SYSTEM","HARVARD UNIVERSITY","CHINESE ACADEMY OF SCIENCES",
                                  "INSTITUT NATIONAL DE LA SANTE ET DE LA RECHERCHE MEDICALE INSERM",
                                  "UNIVERSITY OF NORTH CAROLINA","UNIVERSITY COLLEGE CORK",
                                  "UNIVERSITY OF CALIFORNIA SAN DIEGO","INRAE","UNIVERSITY OF LONDON",
                                  "UNIVERSITY OF TEXAS SYSTEM"),
                       
                       value=c(4.770,3.083,1.821,1.758,1.578,1.495,1.495,1.466,1.389,1.374)
)
treemap(ExampleWOS,index=c("group","subgroup"),vSize="value",type="index",palette = "Set1")
result_KEGG <- read.csv('3.modulefunc/enrichKEGG.csv',check.names = F)
result_KEGG = result_KEGG[1:15,]
result_KEGG$pvalue
mycol <- colorRampPalette(brewer.pal(8,"Set1"))(15)

pdf('3.modulefunc/KEGG_treemap.pdf',width = 8,height = 5)
treemap(result_KEGG,index=c("Description",'Description'),
        vSize="Count",
        title = 'KEGG Enrich Treemap',
        vColor="pvalue",
        type="manual",
        border.col='#63B8FF',
        # palette=c("#FFFFFF00", "#1C86EE00"))
        palette = 'Set2')

dev.off()
##
#################4. LASSO回归分析########



# a <- read.csv('GSE116250-mRNA-ICM.csv',check.names = F,row.names = 1)
# gene = x1
# aa = a[gene,]
# # a = na.omit(a)
# data = t(aa)
# data = data.frame(data)
# data$Type <- c(rep('1',13),rep('0',14))
# write.csv(data,'4.hubgene//input.csv')

# a <- read.csv('GSE5406_mRNA.csv',check.names = F,row.names = 1)
# b <- read.csv('GSE5406-sample.csv',check.names = F,row.names = 1)
# id = rownames(b)
# a = a[,id]
# aa = a[gene,]
# aa = na.omit(aa)
# data = t(aa)
# data = data.frame(data)
# data$Type <- c(rep('1',194),rep('0',16))


library(e1071)
library(glmnet)
library(pROC)
library(caret)


rtt <- read.csv("4.hubgene/input.csv",check.names = F,row.names = 1)
# roc_g <- c('PPDPF','DPEP2','LTBP1','SOCS2','C9orf16','FZD7','CDKN1A')

rtt = rtt[,c(roc_g,'Type')]


roc <- data.frame()
roc1 <- data.frame()
#j=2,i=41
#j=9,i=5

for (j in 1:100) {
  
  set.seed(90)
  index <- caret::createDataPartition(rtt[,"Type"], p =0.7)
  lasso_train <- rtt[index$Resample1,]
  
  lasso_test<- rtt[-index$Resample1,]
  # lasso_test <- read.csv("GSE196582_input.csv",check.names = F,row.names = 1)
  # lasso_test$Type <- c(rep('1',6),rep('0',6))
  # lasso_test=rtt
  
  nx <- as.matrix(lasso_test[,c(1:(ncol(lasso_test)-1))])
  
  x=as.matrix(lasso_train[,c(1:(ncol(lasso_train)-1))])
  y=data.matrix(lasso_train$Type)
  
  
  for (i in 1:10) {
    tryCatch({
      
      set.seed(9)
      fit <- glmnet(x, y, family = "binomial", type.measure="class")
      
      pdf(file = "4.hubgene/lasso_lambda.pdf",width = 9,height = 5)
      par(mfrow=c(1,2))
      plot(fit, xvar = "lambda", label = TRUE)
      #dev.off()
      cvfit <- cv.glmnet(x, y, family="binomial",type.measure = 'class')
      plot(cvfit)
      abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
      dev.off()
      coef <- coef(fit, s = cvfit$lambda.min)
      index <- which(as.matrix(coef)!= 0)
      actCoef <- coef[index]
      lassoGene=row.names(coef)[index]
      lassoGene
      write.csv(lassoGene,"4.hubgene/lassoGene.csv")
      #validation
      out <- predict(cvfit,newx=nx,s='lambda.min',type='class')
      group <- data.frame(Sample=rownames(lasso_test),Type=lasso_test$Type)
      rownames(group) <- group$Sample
      group.test <- group[rownames(out),]
      out <- cbind(data.frame(out,check.names = F),group.test$Type)
      colnames(out) <- c('predict','Type')
      #out$predict <- sapply(out$predict,function(x){as.numeric(x)})
      if(length(unique(out$predict)) == 1) {next}
      out$predict <- factor(out$predict,labels = c(0,1))
      out$Type <- factor(out$Type,labels = c(0,1))
      
      #train
      out1 <- predict(cvfit,s='lambda.min',newx=x,type='class')
      
      group1 <- data.frame(Sample=rownames(lasso_train),Type=lasso_train$Type)
      rownames(group1) <- group1$Sample
      group1.test <- group1[rownames(out1),]
      out1 <- cbind(data.frame(out1,check.names = F),group1.test$Type)
      colnames(out1) <- c('predict','Type')
      #out1$predict <- sapply(out1$predict,function(x){as.numeric(x)})
      if(length(unique(out1$predict)) == 1) {next}
      out1$predict <- factor(out1$predict,labels = c(0,1))
      out1$Type <- factor(out1$Type,labels = c(0,1))
      
      pdf('4.hubgene/lasso.ROC.pdf',width = 12,height = 6)
      par(mfrow=c(1,2))
      
      p1 <- plot.roc(as.numeric(out1[,2]),as.numeric(out1[,1]),ylim=c(0,1),xlim=c(1,0),
                     smooth=F, #smooth curve
                     main="Training Set", auc.polygon=T,
                     auc.polygon.col="cornflowerblue",
                     #print.thres="best", #sensitivity+ specificity
                     col='MidnightBlue',
                     lwd=2, 
                     legacy.axes=T,print.auc=T)
      
      p <- plot.roc(as.numeric(out[,2]),as.numeric(out[,1]),ylim=c(0,1),xlim=c(1,0),
                    smooth=F,
                    main="Validation Set", auc.polygon=T,
                    auc.polygon.col="cornflowerblue",
                    #print.thres="best", 
                    col='MidnightBlue',
                    lwd=2, 
                    legacy.axes=T,print.auc=T)
      print(paste(j,i,as.numeric(p$auc),as.numeric(p1$auc),as.numeric(p$auc)+as.numeric(p1$auc),sep = '\t'))
      roc <- rbind(roc,data.frame(LG = paste(lassoGene[-1],sep = '',collapse = ','),
                                  AUC.test= as.numeric(p$auc),
                                  AUC.train= as.numeric(p1$auc),
                                  seed1=j,seed2=i,sum=as.numeric(p$auc)+as.numeric(p1$auc)))
      # roc1 <- rbind(roc1,data.frame(LG = paste(lassoGene[-1],sep = '',collapse = ','),
      #                             AUC.test= as.numeric(p$auc),
      #                             AUC.train= as.numeric(p1$auc),
      #                             seed1=j,seed2=i))
      
      dev.off()
    },error=function(e){})
    
  }
  
}


pdf('lasso.ROC.pdf',width = 4.5,height = 4.5)


p1 <- plot.roc(as.numeric(out1[,2]),as.numeric(out1[,1]),ylim=c(0,1),xlim=c(1,0),
               smooth=F, 
               main="GSE21321 ROC Curve", auc.polygon=T,
               auc.polygon.col="cornflowerblue",
               #print.thres="best",
               col='MidnightBlue',
               lwd=3, 
               legacy.axes=T,print.auc=T)
dev.off()

write.csv(roc,"4.hubgene/roc.csv")

###
#################4. SVM-RFE#######

library(tidyverse)
library(glmnet)
source('msvmRFE.R')  
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)


input<- read.csv('4.hubgene//svm-input.csv',check.names = F,row.names = 1)
# (k-fold crossValidation）

#fix(input)
set.seed(99)

svmRFE(input, k = 5, halve.above = 100) #seed
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #features

top.features = WriteFeatures(results, input, save=F) 
head(top.features)

#output
write.csv(top.features,"4.hubgene//svm.csv")

featsweep = lapply(1:99, FeatSweep.wrap, results, input) #19个
featsweep


no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

# dev.new(width=4, height=4, bg='white')
pdf("4.hubgene//svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()

pdf("4.hubgene//svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 
dev.off()

svm_gene <- top.features[1:which.min(errors), "FeatureName"]
write.csv(svm_gene,'4.hubgene//svm_gene.csv')

##############4. hubgene-venn######
a <- read.csv('4.hubgene//svm_gene.csv',check.names = F)
b <- read.csv('4.hubgene/lassoGene.csv',check.names = F)

pdf("4.hubgene//hubgene-venn.pdf", width = 6, height = 6)

venn::venn(list(`SVM-RFE`= a$x,`LASSO`= b$x), 
           zcolor= c('blue','red','darkgreen'),box=T,sncs = 1.5,ilcs = 1.5)
dev.off()
x1 = intersect(a$x,b$x)
write.csv(x1,'4.hubgene//hubgene.csv')


x = list(`SVM-RFE`= a$x,`LASSO`= b$x)

pdf('4.hubgene/hubgene-venn.pdf',width = 6, height = 4)
ggVennDiagram(x) + scale_fill_gradient(low="#0099b499",high= "#ed000099")
dev.off()
####
##############4. ROC######

library(pROC)
data <- read.csv('4.hubgene/input.csv',check.names = F,row.names = 1)
data1 = data[,x1]
hubexp = data1
hubexp$Type = data$Type
rocdata <- data.frame(Sample = rownames(hubexp),
                      exp = hubexp$PPDPF,
                      Type = hubexp$Type)
#mycol <- brewer.pal(10,'Set3')
#install.packages("RColorBrewer")
library(RColorBrewer)
mycol <- brewer.pal(9,'Set1')

pdf("4. hubgene/hubgene-ROC.pdf",width = 6,height = 6)

x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, 
              ci = F,
              main="ROC Curve",
              #print.thres="best",
              col=mycol[3],
              lwd=3, 
              # print.auc = T,
              legacy.axes=T)

j=1
auc.test <- paste0(colnames(hubexp)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(hubexp[,2:(ncol(hubexp)-1)])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp[,i],
                        Type = hubexp$Type)
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, 
                ci = F,
                main="ROC Curve",
                #print.thres="best", 
                col=mycol[j],
                lwd=3, 
                # print.auc = T,
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',
                       format(as.numeric(x$auc),digits=2)))
}
legend(0.55,0.5, auc.test,lwd=3,bty="n",col=mycol,cex = 1.2)
dev.off()




data <- read.csv('4.hubgene/input.csv',check.names = F,row.names = 1)

hub  = x1

for( i in hub){
  
  # i = "RFC4"
  hubexp = data.frame(exp = data[,i],Type=data$Type)
  
  hubexp <- na.omit(hubexp)#delete namit
  
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp$exp,
                        Type = hubexp$Type) 

  mycol <- brewer.pal(8,'Set1')
  
  pdf(file = paste0('4.hubgene/GSE116250-ROC/',i,'-ROC.pdf'),width = 4.5,height = 4.5)
  
  x <-  plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                 smooth=F,
                 ci = T,
                 main=paste0('GSE116250 ',i," ROC Curve"),
                 print.thres="best", 
                 col=mycol[1],
                 lwd=3, 
                 print.auc = T,
                 legacy.axes=T)
  
  # auc.test <- paste0(gene,' AUC : ',format(as.numeric(x$auc),digits=2))
  # 
  # legend(0.5,0.5, x,lwd=2,bty="n",col=mycol[1],cex = 1.2)
  # 
  dev.off()
  
}

library(ggpubr)

data <- read.csv('4.hubgene//input.csv',check.names = F,row.names = 1)

# data[,1:303] = apply(data[,1:303],2,function(x){log2(x+1)})
dbox1 = data[,hub]
dbox1 = apply(dbox1,2,function(x){log2(x+1)})
dbox1 = data.frame(dbox1)
dbox1$level <- ifelse(data$Type == 0,'Control','ICM')


table(dbox1$level)



######

library(ggplot2)
library(reshape2)
library(plyr)
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))

source('Function_for_violin_plot.R')

dbox1 = gather(dbox1,gene,expression,1:(ncol(dbox1)-1))

Data_summary <- summarySE(dbox1, measurevar="expression", groupvars=c("level","gene"))
head(Data_summary)



if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}


ggplot(dbox1,aes(x= gene,y= expression,fill= level))+
  geom_split_violin(trim= F,color="white",scale = "area") + 
  
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ 
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.1, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c( "#0099b499","#ed000099"))+ 
  labs(y=("Normalized Expression"),x=NULL,title = "") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = level),
                     label = "p.signif",size = 5,
                     method = "wilcox.test",
                     label.y = max(dbox1$expression),
                     hide.ns = F)

ggsave(filename = "9. hubgene-venn//GSE116250-violin.pdf",height = 4,width = 6)
###

##############4. ROC######

library(pROC)
data <- read.csv('4.hubgene/GSE5406_input.csv',check.names = F,row.names = 1)


hubexp = data
hubexp$Type = data$Type
rocdata <- data.frame(Sample = rownames(hubexp),
                      exp = hubexp$CX3CL1,
                      Type = hubexp$Type)
#mycol <- brewer.pal(10,'Set3')
#install.packages("RColorBrewer")
library(RColorBrewer)
mycol <- brewer.pal(9,'Set1')

pdf("4.hubgene/GSE112943-ROC.pdf",width = 6,height = 6)

x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              ci = F,
              main="GSE112943 ROC Curve",
              #print.thres="best", 
              col=mycol[3],
              lwd=3, 
              # print.auc = T,
              legacy.axes=T)

j=1
auc.test <- paste0(colnames(hubexp)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(hubexp[,2:(ncol(hubexp)-1)])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp[,i],
                        Type = hubexp$Type)
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, 
                ci = F,
                main="GSE112943 ROC Curve",
                #print.thres="best", 
                col=mycol[j],
                lwd=3, 
                # print.auc = T,
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',
                       format(as.numeric(x$auc),digits=2)))
}
legend(0.55,0.5, auc.test,lwd=3,bty="n",col=mycol,cex = 1.2)
dev.off()

# write.csv(auc.test,'4.hubgene//testing-ROC.csv')



data <- read.csv('4.hubgene/GSE5406_input.csv',check.names = F,row.names = 1)

hub = x1

for( i in hub){
  
  # i = "RFC4"
  hubexp = data.frame(exp = data[,i],Type=data$Type)
  
  hubexp <- na.omit(hubexp)
  
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp$exp,
                        Type = hubexp$Type) 
  #mycol <- brewer.pal(10,'Set3')
  #install.packages("RColorBrewer")
  # dir.create('4.hubgene/GSE5406-单基因ROC')
  
  mycol <- brewer.pal(8,'Set1')
  
  pdf(file = paste0('4.hubgene/GSE5406-ROC/',i,'-ROC.pdf'),width = 4.5,height = 4.5)
  
  x <-  plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                 smooth=F, 
                 ci = T,
                 main=paste0('GSE5406 ',i," ROC Curve"),
                 print.thres="best", 
                 col=mycol[2],
                 lwd=3, 
                 print.auc = T,
                 legacy.axes=T)
  
  # auc.test <- paste0(gene,' AUC : ',format(as.numeric(x$auc),digits=2))
  # 
  # legend(0.5,0.5, x,lwd=2,bty="n",col=mycol[1],cex = 1.2)
  # 
  dev.off()
  
}
###boxplot
library(ggpubr)

data <- read.csv('4.hubgene//GSE5406_input.csv',check.names = F,row.names = 1)


dbox1 = data[,c(hub)]
dbox1$level <- ifelse(data$Type == 0,'Control','ICM')


table(dbox1$level)



dbox1 = gather(dbox1,gene,expression,1:(ncol(dbox1)-1))
##caculate mean and error
Data_summary <- summarySE(dbox1, measurevar="expression", groupvars=c("level","gene"))
head(Data_summary)



if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}



ggplot(dbox1,aes(x= gene,y= expression,fill= level))+
  geom_split_violin(trim= F,color="white",scale = "area") + 
  
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ 
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.1, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#0099b499","#ed000099"))+ 
  labs(y=("Normalized Expression"),x=NULL,title = "") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = level),size = 5,
                     label = "p.signif",
                     method = "wilcox.test",
                     label.y = max(dbox1$expression),
                     hide.ns = F)

ggsave(filename = "9. hubgene-expression//GSE5406-violin.pdf",height = 4,width = 6)
###


#################5. hubgene clinical application########
library("foreign")
# library(Matrix)
library(rms)
library(regplot)
# dir.create('5.hubgene-clinical')

dat <- read.csv('4.hubgene//input.csv',check.names = F,row.names = 1)
hub = x1
dat = dat[,c('Type',hub)]

train_data = dat[,c(1:5)]
f2 <- lrm(formula = Type ~ ., 
          data = train_data, x = T, y = T) 

regplot(f2,
        # observation=nomogram.data[2, ], 
        points = TRUE, #If FALSE the regression scores of each βx contribution are shown. Otherwise contributions are represented by a 0-100 "points" scale.
        title = "Nomogram",
        odds = TRUE,
        showP = TRUE,
        droplines=TRUE,
        # observation = train_data[1,],
        interval="confidence") 

##ROC############

# dt$y = ifelse(dt$y == 'Control',0,1)
mydata = dat
mydata = data.frame(mydata)
y <- mydata$Type

PPDPF <- mydata$PPDPF
DPEP2 <- mydata$DPEP2
LTBP1 <- mydata$LTBP1
SOCS2 <- mydata$SOCS2


dt <- data.frame(y,PPDPF,DPEP2,LTBP1,SOCS2)
# dt$y = as.numeric(dt$y)
str(dt)

f.glm <- glm(y~PPDPF+DPEP2+LTBP1+SOCS2,data=dt,family = binomial(link = "logit"))


P1 <- predict(f.glm,type = 'response')  

pdf(file = "4.hubgene/Calibration curve.pdf",width = 5.5,height = 5.5)
val.prob(P1,y)     ##P1,probit;y,Outcome variable
dev.off()



#################6. CIBERSORT#######
library("e1071")
library("preprocessCore")
library("GSVA")
library("estimate")
library("ggpubr")
library("tidyr")
library("ggplot2")
library("pheatmap")
source("CIBERSORT.R")

# dir.create('6.immune-analysis')
#CIBERSORT
set.seed(1)
results=CIBERSORT("LM22.txt", 'GSE116250-mRNA-ICM.txt', perm=100, QN=F)

#p<0.05
output <- read.table('CIBERSORT-Results.txt',sep='\t',check.names = F,header = T)
# output <- output[output$`P-value`<0.05,]
output = output[,-c(26,24,25)]
write.csv(output,"6.immune-analysis//CIBERSORT-filter.csv")



data <- read.csv("CIBERSORT-filter.csv",check.names=F,row.names=1)
data = data[,c(1:22)]
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf('2-Cluster-CIBERSORT-proportion.pdf',height=7,width=14)
#png('infi_heatmap.png',height=800,width=1200)
par(las=1,mar=c(4,4,4,10))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F,tick = F)
par(xpd=T);text(100,0,'TCGA CHOL (N = 31)',adj=0.9,cex=1.7);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),
       col=col,pch=15,bty="n",cex=1.2)

dev.off()

#################6.cell proportion#######
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot","ggpubr","bslib","ggthemes")
# BiocManager::install("ggthemes")
library(ggthemes)
library(dplyr)
library(tidyr)
library(tibble)
# Read in results
####
cibersort_raw <- read.csv("6.immune-analysis/mcp_counter.csv",row.names = 1,check.names = F)


head(cibersort_raw)
# View(cibersort_raw)
dd1 <- cibersort_raw %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:9,
               names_to = "CellType",
               values_to = "Composition")
# dd1[,3] = apply(dd1[,3],2,function(x){log2(x+1)})
View(dd1)
plot.info <- dd1[,c(1,3,2)]

mycol = colorRampPalette(brewer.pal(8,"Set1"))(34)
mycol = colorRampPalette(brewer.pal(8,"Set3"))(8)


ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",palette = mycol,
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "") +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))
ggsave(filename = "6.immune-analysis/boxplot.pdf",width = 6,height = 4)




ggbarplot(
  plot.info,
  x = "sample",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType",
  palette = 'Set3'
  
) +
  theme_base() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 1
    ),
    legend.position = "bottom"
  )
ggsave(filename = "6-ssGSEA-barplot.pdf",width = 14,height = 8)


##################6. Differential immune cells#####
dbox1 <- read.csv('6.immune-analysis/mcp_counter.csv',check.names = F,row.names = 1)
dbox1$level = c(rep('ICM',13),rep('Control',14))

table(dbox1$level)

# my_comparition <- list(c("LN","control"))

# mycol <- c("#E81D22","#02b1e6")
dbox1 = gather(dbox1,gene,expression,1:(ncol(dbox1)-1))

pdf('6.immune-analysis//CIBERSORT-boxplot.pdf',width = 7,height = 4)

p <- ggboxplot(dbox1,x='gene',y='expression',fill ='level',
               ylab = 'Proportion ',
               xlab ='',palette = mycol,
               # add ='jitter',
               size =0.4)+
  rotate_x_text(45)+
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        axis.line = element_line(size=0.3))+
  stat_compare_means(size = 4.5,aes(group=level),
                     label = 'p.signif',label.x = 1.3,method = 'wilcox.test',label.y = 29)

p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
  font('x.text',face = 'bold')+font('y.text',face = 'bold')+
  font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

dev.off()


Data_summary <- summarySE(dbox1, measurevar="expression", groupvars=c("level","gene"))
head(Data_summary)



if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}



ggplot(dbox1,aes(x= gene,y= expression,fill= level))+
  geom_split_violin(trim= F,color="white",scale = "area") + 
  
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ 
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.1, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c( "#8DD4B5","#FBABFB"))+ 
  labs(y=("Normalized Expression"),x=NULL,title = "") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = level),
                     label = "p.signif",
                     method = "wilcox.test",
                     label.y = max(dbox1$expression),
                     hide.ns = F)

ggsave(filename = "4.hubgene/GSE5406-violin.pdf",height = 4,width = 6)


##################6. Differential immune cells########

dbox1 <- read.csv('6.immune-analysis/mcp_counter.csv',check.names = F,row.names = 1)
cell = colnames(dbox1)
data = dbox1
data$level = c(rep('ICM',13),rep('Control',14))
pl <- list()
for(i in cell){
  
  # i = 'T cell'
print(i)
pl[[i]] <- ggboxplot(data, x = "level", y = i, color  = "level", 
                     size = 0.5,palette = "npg",notch = TRUE,
                     add = "dotplot",add.params=list(size=1)) +
  labs(y = paste0(i," content"), x = "")+
  rotate_x_text(angle =0 )+
  stat_compare_means(aes(group = level),label = "p.signif",size =5,label.x = 1.5)
}

pdf('6.immune-analysis/boxplot.pdf',width = 10,height = 6)
ggarrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]],
          pl[[5]],pl[[6]],pl[[7]],pl[[8]],
         
          print = TRUE,nrow = 2,ncol = 4)

dev.off()





###################6. Correlation between key genes and immune cells######


b<-read.csv("6.immune-analysis/mcp_counter.csv",check.names = F,row.names =1)
a <- read.csv('GSE116250-mRNA-ICM.csv',check.names = F,row.names = 1)
tcga_gsva = b
# View(tcga_gsva)
### read expression
gene = hub

aa = a[gene,]
tcga_expr = aa

tcga_expr = as.matrix(tcga_expr)


### read gene
# genelist <- read.csv("m5C_expr_DE.csv",check.names = F)
genelist <-gene

gene <- genelist
immuscore <- function(gene){
  y <- tcga_expr[gene,]
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
# immuscore("")
### Batch calculation of the results of the correlation between genelist and immune infiltration
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
# p = 0.01
# cor = 0.6
# data1 = data[(data$p.value < p & (data$cor > cor |  data$cor < -(cor))),]
write.csv(data, "6.immune-analysis//cor.csv", quote = F, row.names = F)

data$pstar <- ifelse(data$p.value < 0.05&abs(data$cor) >0,
                     ifelse(data$p.value < 0.01&abs(data$cor) >0,"**","*"),
                     "")
data$pstar[1:20]

ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "royalblue",mid = "white",high = "#cc340c")+
  geom_text(aes(label=pstar),col ="midnightblue",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))

ggsave("9.immune-analysis/cor.pdf", width = 6, height = 4)


###基因相关性圈图

blue <- "#338CF4"
red <- "#ff0c0a"
orange <- "#ff0c0a"
corCnaExpr = data
# corCnaExpr <- read.csv("file:///E:/22年6月项目/58. YQXX0251/5. 免疫浸润/CIBERSORT/关键基因与CIBERSORT免疫细胞相关性.csv",check.names = F)
corCnaExpr = na.omit(corCnaExpr)


my_palette <- colorRampPalette(c(blue,"white",orange), alpha=TRUE)(n=128)

ggplot(corCnaExpr, aes(x=gene,y=immune_cells)) +
  geom_point(aes(size=-log10(p.value), color=cor)) +
  scale_color_gradientn('Correlation', 
                        colors=my_palette) + 
  scale_size_continuous(range = c(1,12)) + 
  
  theme_bw() +
  theme(#panel.grid.minor = element_blank(), 
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 12, color = rep(c(blue,blue),c(26,10))),
    axis.title = element_blank(),
    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
    legend.position = "bottom",
    plot.margin = unit(c(1,1,1,1), "lines"))

ggsave("模型基因与ssGSEA差异细胞相关性.pdf", width = 8,height = 9)


####棒棒糖图

library(ggplot2)

Corr <- read.csv("6.immune-analysis/关键基因与免疫细胞相关性.csv",check.names = F)


p <- list()
for( i in gene){
  
  # i = "PPDPF"
  
  cor = Corr[Corr$gene == i,]
  cor = cor[order(cor$cor),]
  cor$Correlation = cor$cor
  cor = cor[order(cor$Correlation),]
  
  # pdf(file = paste0('9.immune-analysis//',i,'与免疫细胞相关性.pdf'),width=5,height=5)
  print(i)
  p[[i]] <- ggplot(cor,aes(x = Correlation,y=reorder(immune_cells,Correlation),fill=p.value)) +
    ylab('immune cell')+
    geom_segment(aes(xend=0,yend=immune_cells))+
    
    geom_point(shape=19,aes(size = Correlation,colour=p.value))+
    
    scale_fill_gradient2(mid = "#ffcc25",high = "darkgreen")+
    scale_color_gradient2(mid = "#ffcc25",high = "darkgreen") +
    
    scale_size_continuous()+
    ggtitle(paste0(i))+
    theme(
      axis.title=element_text(size=13,face="plain",color="black"),
      axis.text = element_text(size=13,face="plain",color="black"),
      legend.title=element_text(size=12,face="plain",color="black"),
      legend.background = element_blank(),
      legend.position = "left"
    )+theme_bw()
  # theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
  
  # print(p)
  # 
  # dev.off()
  
}

pdf(file = '6.immune-analysis/关键基因与免疫细胞棒棒糖图.pdf',width = 8,height = 7) 
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],
          print = TRUE,nrow = 2,ncol = 2)
# print(box,nrow = ,ncol = 4)
# graph_one
dev.off()


for( i in gene){
  
  # i = "PTX3"  
  data = Corr[Corr$gene==i,]
  data <- data %>% 
    mutate(Correlation = ifelse(cor>0, "1-positive", "negative")) #setting group
  
  fig4 <- ggplot(data, aes(x=reorder(immune_cells,cor), y=cor,fill = Correlation)) +

    
    geom_point(shape=21,aes(size = -log10(p.value),colour=Correlation)) + 
    
    geom_text(aes(label =sprintf("%.2f",cor)), color = "black", size = 2.4)+ 
    xlab("Immunocytes") +
    scale_y_continuous(paste0(i,"-Correlation"),breaks  = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)) +
    coord_flip() +
    theme_minimal()+
    theme(
      axis.title=element_text(size=13,face="plain",color="black"),
      axis.text = element_text(size=13,face="plain",color="black"),
      legend.title=element_text(size=12,face="plain",color="black"),
      legend.background = element_blank(),
      legend.position = "left"
    )+theme_bw()
  
  
  print(fig4)
  
  ggsave(file=paste0(i,"与ssGSEA免疫细胞相关棒棒糖图.pdf"),width = 5,height = 5)
  dev.off()
  
}



###################6. Differential immune cell correlation scatter plot#######

library(ggpubr)
library(ggExtra)
library(scales)
# install.packages('ggExtra')

rt = read.csv('6.immune-analysis/mcp_counter.csv',check.names = F)

x1 <- c('T cell','T cell CD8+')

box = list()
for(i in x1){
  # i = 'REPS2'
  df1=as.data.frame(rt[,c("Neutrophil",i)])
  colnames(df1) <- c("y","x")
  p1=ggplot(df1, aes(x, y)) + 
    xlab(paste0(i)) + ylab("Neutrophil")+ 
    geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  print(i)
  box[[i]] <- ggMarginal(p1, type="density", xparams=list(fill = "#ed000099"),
                         yparams=list(fill = "#0099b499"))
  
  # pdf(file=paste0("6-REPS2","-Neutrophil-cor.pdf"), width=4.2, height=4)
  # print(p2)
  # dev.off()
  
}

x2 <- c('T cell CD8+')
aox = list()
for(i in x2){
  # i = 'REPS2'
  df1=as.data.frame(rt[,c("T cell",i)])
  colnames(df1) <- c("y","x")
  p1=ggplot(df1, aes(x, y)) + 
    xlab(paste0(i)) + ylab("T cell")+ 
    geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  print(i)
  aox[[i]] <- ggMarginal(p1, type="density", xparams=list(fill = '#ed000099'),
                         yparams=list(fill = "#0099b499"))
  
  # pdf(file=paste0("6-REPS2","-Neutrophil-cor.pdf"), width=4.2, height=4)
  # print(p2)
  # dev.off()
  
}


pdf(file="6.immune-analysis//scatterplot.pdf",width = 9,height = 3)

ggarrange(aox[[1]],
          box[[1]],box[[2]],
          print = TRUE,nrow = 1,ncol = 3)
dev.off()

###
##################7. miRNA######

# dir.create('10. miRNA-TF-network')
library(multiMiR)

###predict mRNA
example1 <- get_multimir(org = "hsa",           #  hsa/mmu/rno
                         mirna = 'hsa-miR-182-5p', 
                         table = "validated",   # validated/predicted/all
                         summary = TRUE)        


#####predict miRNA
target = hub
example3 <- get_multimir(org = "hsa",                   #  hsa/mmu/rno
                         target = target,               # target gene
                         table = "predicted",           # miRNA-target
                         summary = TRUE, 
                         predicted.cutoff = 30,          # score TOP 35%
                         predicted.cutoff.type = "p",    # score TOP 35%
                         predicted.site = "all")         # conserved/nonconserved/all
table(example3@data$type)  
example3_result <- example3@data
head(example3_result)
table(example3_result$database)
write.csv(example3_result,'7.hubgene-network//miRNA-gene.csv')

###Venn

a<- read.csv('7.hubgene-network/miRNA-gene.csv',check.names = F)
table(a$database)
a1 = a[a$database =='mirdb',]
a2 = a[a$database == 'miranda',]
a3 = a[a$database =='targetscan',]

pdf("10. miRNA-TF/miRNA-venn.pdf", width = 6, height = 6)

venn::venn(list(`mirDB`= a1$mature_mirna_id,
                `miranda`= a2$mature_mirna_id,
                `targetscan`= a3$mature_mirna_id),
           zcolor= c("#3F60AA","#CC340C",'darkgreen'),box=T,ilabels = ,sncs = 2,ilcs = 2,
           opacity = 0.3)
dev.off()

x1 = intersect(a1$mature_mirna_id,a2$mature_mirna_id)
x2 = intersect(x1,a3$mature_mirna_id)
aa = a[a$mature_mirna_id %in% x2,]

lnc <- read.csv('d:/human.all6.lnc.csv',check.names = F)
aa$miRNA = aa$mature_mirna_id
lnc_mi = merge(aa,lnc,by ='miRNA')

table(lnc_mi$target_symbol)
write.csv(lnc_mi,'7.hubgene-network//lncrna_miRNA-gene.csv')


##################7. DE-TF#####
a <- read.csv('7.hubgene-network/TF-gene.csv',check.names = F)
b <- read.csv('2.modulegene/limmaOut.csv',check.names = F,row.names = 1)
b = b[b$adj.P.Val < 0.05,]

ab = a[a$node2 %in% rownames(b),]
write.csv(ab,'7.hubgene-network/TF-gene.csv')
