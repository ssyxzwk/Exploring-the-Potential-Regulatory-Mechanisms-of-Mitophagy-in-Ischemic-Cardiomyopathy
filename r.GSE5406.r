###验证集GSE5406---------------------------------
gset<-getGEO("GSE5406",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL96",destdir = '/data/nas1/dailihui/pipeline/GPL/')
# gpl <- GPL96_soft
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
gpl$`Gene Symbol`
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('Gene Symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
colnames(probe2symobl) <- c("ID","symbol")
probe2symobl=probe2symobl[probe2symobl$symbol!='',]

dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%
  dplyr::select(symbol,everything())%>%
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%
  arrange(desc(rowMean))%>%
  distinct(symbol,.keep_all = T)%>%
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames(colnames(.)[1])

###mrna
protein <- read.table("/data/nas1/dailihui/pipeline/PCG/protein_coding.txt",header = T,sep = "\t")
dat <- dat[rownames(dat)%in%protein$gene_name,]

###分组信息
pd<-pData(a)
head(pd)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
table(group$group)
group$group[grep("Ischemic", group$group)] <- "Ischemic"
group$group[grep("NonFailing", group$group)] <- "Control"
group$group[grep("Idiopathic", group$group)] <- "Idiopathic"

library(stringr)
group <- subset(group,group%in%c("Ischemic","Control"))
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.csv(dat,file = 'dat(GSE5406).csv')
write.csv(group,file = 'group(GSE5406).csv')



hub.genes = data.frame(symbol=c('PPDPF','DPEP2','LTBP1','SOCS2'))
# rownames(df.exp) <- 

hub.exp = dat[as.character(hub.genes$symbol),]%>%lc.tableToNum()

hub.exp$gene = rownames(hub.exp)

library(data.table)
df.exp.plot = melt(hub.exp, id.vars = "gene", variable.name = "sample", value.name = "value")
df.exp.plot$group = group$group[match(df.exp.plot$sample, group$sample)]

df.diff = lc.ttest(hub.exp[-ncol(hub.exp)], group$group)$Result
df.diff$p = df.diff$`Control-Ischemic.p.value` %>% signif(3) %>% paste0("p = ", .)
df.diff$label = paste0(rownames(df.diff),"\n",df.diff$p)

df.exp.plot = lc.tableToNum(df.exp.plot)
df.exp.plot$label = df.diff$label[match(df.exp.plot$gene, rownames(df.diff))]
df.exp.plot$group = factor(df.exp.plot$group, levels = c("Control", "Case"))



library(pROC)
res.roc = list()
hub.exp = hub.exp[group$sample] %>% t %>% as.data.frame() %>% lc.tableToNum()
colnames(hub.exp) = hub.genes$symbol

apply(hub.exp, 2, function(x){
  res = roc(group$group, x)
  res.roc[[length(res.roc) + 1]] <<- res
})
names(res.roc) = colnames(hub.exp)

res.auc = sapply(res.roc, auc) %>% data.frame
res.auc = cbind(Gene = rownames(res.auc), AUC = res.auc[[1]]) %>% as.data.frame %>% lc.tableToNum
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
ggsave("02.Validation_GSE5406.ROC.png", width = 9, height = 9, dpi = 300, units = "in", bg = "white")
ggsave("02.Validation_GSE5406.ROC.pdf", width = 9, height = 9, dpi = 300, units = "in", bg = "white")
