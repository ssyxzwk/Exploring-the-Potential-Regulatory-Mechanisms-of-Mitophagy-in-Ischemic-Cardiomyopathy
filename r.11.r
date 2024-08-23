setwd("/data/nas1/dailihui/re-project/01.SJZZK-451-3")

rm(list = ls())

library(dplyr)

data <- read.csv('./input.csv',check.names = F,row.names = 1)

df.exp <- t(data) %>% as.data.frame()



# df.exp <- read.csv("./dat.GSE42955.csv")

a <- colnames(df.exp)
b <- ifelse(grepl("ICM",a),"ICM","NF")

df.pheno <- data.frame(sample=a,group=b)

hub.genes = data.frame(symbol=c('PPDPF','DPEP2','LTBP1','SOCS2'))
# rownames(df.exp) <- 
# df.exp$gene <- NULL
hub.exp = df.exp[as.character(hub.genes$symbol),]

hub.exp$gene = rownames(hub.exp)

library(data.table)
df.exp.plot = melt(hub.exp, id.vars = "gene", variable.name = "sample", value.name = "value")
df.exp.plot$group = df.pheno$group[match(df.exp.plot$sample, df.pheno$sample)]



# df.diff = lc.ttest(hub.exp[-ncol(hub.exp)], df.pheno$group)$Result
# df.diff$p = df.diff$`Control-Case.p.value` %>% signif(3) %>% paste0("p = ", .)
# df.diff$label = paste0(rownames(df.diff),"\n",df.diff$p)



hub.exp$gene <- NULL
df.pheno$group <- factor(df.pheno$group,levels = c("ICM","HF"))

library(limma)

# 1. 创建设计矩阵（假设 df.pheno$group 是一个因子变量，包含 "Control" 和 "Case"）
design <- model.matrix(~0 + df.pheno$group)
colnames(design) <- levels(df.pheno$group)

# 2. 进行线性模型拟合
fit <- lmFit(hub.exp, design)

# 3. 创建对比矩阵（这里进行 "Case" 对比 "Control" 的比较）
contrast.matrix <- makeContrasts(ICM-HF, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# 4. 应用贝叶斯方法计算差异表达的统计量
fit2 <- eBayes(fit2)

# 5. 提取结果
df.diff <- topTable(fit2, adjust="BH", number=Inf)

# 6. 添加 p 值标注
df.diff$p = signif(df.diff$P.Value, 3) %>% paste0("p = ", .)
df.diff$label = paste0(rownames(df.diff), "\n", df.diff$p)

# 7. 打印结果
print(df.diff)
















df.exp.plot$label = df.diff$label[match(df.exp.plot$gene, rownames(df.diff))]
df.exp.plot$group = factor(df.exp.plot$group, levels = c("NF", "ICM"))
library(ggplot2)
ggplot() +
  geom_boxplot(data = df.exp.plot, aes(x = group, y = value, color = group)) + 
  facet_wrap(label ~ ., scales = "free", ncol = 4) +
  theme_bw() + 
  xlab(NULL) + 
  ylab("Expression") + scale_color_manual(values = c("darkgreen", "darkorange")) +
  theme(text = element_text(size = 16), legend.title = element_blank())
ggsave("01.Validation11.png", width = 12, height = 5, units = "in", dpi = 300)
ggsave("01.Validation11.pdf", width = 12, height = 5, units = "in", dpi = 300)


library(pROC)
res.roc = list()
hub.exp = hub.exp[df.pheno$sample] %>% t %>% as.data.frame() 
colnames(hub.exp) = hub.genes$symbol

df.pheno$group <- ifelse(df.pheno$group=="ICM","case","control")

# 创建一个空的列表来存储 ROC 结果
res.roc <- list()

# 使用 apply 函数来计算每列的 ROC 曲线
apply(hub.exp, 2, function(x) {
  res <- roc(df.pheno$group, x)
  res.roc[[length(res.roc) + 1]] <<- res
})

names(res.roc) = colnames(hub.exp)

res.auc = sapply(res.roc, auc) %>% data.frame
res.auc = cbind(Gene = rownames(res.auc), AUC = res.auc[[1]]) %>% as.data.frame 
res.auc$AUC = round(res.auc$AUC, 3)
res.auc.plot = res.auc
res.auc.plot = res.auc.plot[order(res.auc.plot$AUC, decreasing = T),]

library(patchwork)
library(gridExtra)
p1 = ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + scale_color_manual(values = 1:28, guide = guide_legend(ncol = 4)) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"), 
        legend.position = "bottom", legend.text = element_text(size = 13))

p2 = ggplot() + 
  geom_point(aes(x = c(0,1), y = c(0,1)), color = "transparent") +
  annotation_custom(grob = tableGrob(res.auc.plot, rows = NULL), xmin = 0, xmax = 0.1, ymin = -0.2, ymax = 1) +
  xlim(c(0,0.1)) + ylim(c(0,1)) + theme_bw() + 
  theme(text = element_blank(), line = element_blank(), rect = element_blank())

p = p1 + p2
p
ggsave("02.Validation.ROC.png", width = 9, height = 9, dpi = 300, units = "in", bg = "white")
ggsave("02.Validation.ROC.pdf", width = 9, height = 9, dpi = 300, units = "in", bg = "white")
