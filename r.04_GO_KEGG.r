rm(list = ls())
setwd(('/data/nas1/dailihui/re-project/01.SJZZK-451-3'))
if (! dir.exists("./04_GO_KEGG")){
  dir.create("./04_GO_KEGG")
}
setwd("./04_GO_KEGG")
getwd()


library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(tidyverse)

##通路富集分析------

DEGs1 <- read.csv('../up_gene.csv',row.names=1)


gene_transform <- bitr(rownames(DEGs1),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")

# GO富集------
gc()
ego1 <- enrichGO(gene = gene_transform$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pAdjustMethod = "none",   #数据不足把BH改为NONE，即不做p矫正
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1,
                 readable = TRUE)

go_result1 <- data.frame(ego1)
table(go_result1$ONTOLOGY)

#    BP  CC  MF 
# 252  30  34 
write.csv(go_result1,file = "01.GO.csv",row.names = T)
go_result <- go_result1

#提取各组数据需要展示的数量#
display_number = c(5, 5, 5)  ##这三个数字分别代表选取的BP、CC、MF的数量
BP <- go_result[which(go_result$ONTOLOGY=='BP'),]
CC <- go_result[which(go_result$ONTOLOGY=='CC'),]
MF <- go_result[which(go_result$ONTOLOGY=='MF'),]
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]

#将提取的各组数据进行整合
go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),  #指定ego_result_BP、ego_result_CC、ego_result_MFID为ID                        
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  Count=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count), #指定ego_result_BP、ego_result_CC、ego_result_MF的Count为GeneNumber
  type=factor(c(rep("Biological Process", display_number[1]), #设置biological process、cellular component、molecular function 的展示顺序
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))
##设置GO term名字长短，过长则设置相应长度
for(i in 1:nrow(go_enrich)){
  description_splite=strsplit(go_enrich$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")  #选择前五个单词GO term名字
  go_enrich$Description[i]=description_collapse
  go_enrich$Description=gsub(pattern = "NA","",go_enrich$Description)  #gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。gsub(“目标字符”, “替换字符”, 对象)
}
go_enrich <- na.omit(go_enrich)
#转成因子，防止重新排列
go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)
head(go_enrich)
p <- ggplot(go_enrich,
            aes(x=type_order,y=Count, fill=type)) +  #x、y轴定义；根据type填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("cadetblue1", "#FF6666", "#b2e7cb") ) + #柱状图填充颜色
  xlab("") + #x轴标签
  ylab("Counts") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() +
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 75,vjust = 1, hjust = 1 )) 
p
ggsave(filename = '02._up_GO.pdf',p,w=10,h=6)
ggsave(filename = '02.up_GO.png',p,w=10,h=6,units = 'in',dpi = 300)
dev.off()

## KEGG富集分析------
##KEGG气泡图---
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",   #大鼠的富集分析通路库
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "KEGG.csv",quote = F,row.names = F)  ##真正的输出结果3条
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
barplot(kk,showCategory = 10)

ggsave('02.KEGG_dot.png',kk_dot,width =8,height = 6)
ggsave('02.KEGG_dot.pdf',kk_dot,width =8,height = 6)

kk_bar<- barplot(kk,showCategory = 10)
kk_bar

ggsave('02.KEGG_bar.png',kk_bar,width =8,height = 6)
ggsave('02.KEGG_bar.pdf',kk_bar,width =8,height = 6)

# cnetplot(kk,categorySize="pvalue",foldChange = gene_transform$ENTREZID,colorEdge = T)
kk_net <- cnetplot(kk,foldChange = gene_transform$ENTREZID,circular=T,colorEdge=T)

ggsave('02.KEGG_net.png',kk_net,width =8,height = 6)
ggsave('02.KEGG_net.pdf',kk_net,width =8,height = 6)



##down------
DEGs2 <- read.csv('../down_gene.csv',row.names=1)


gene_transform <- bitr(rownames(DEGs2),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")

# GO富集------
gc()
ego1 <- enrichGO(gene = gene_transform$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pAdjustMethod = "none",   #数据不足把BH改为NONE，即不做p矫正
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1,
                 readable = TRUE)

go_result1 <- data.frame(ego1)
table(go_result1$ONTOLOGY)

#    BP  CC  MF 
# 137  25  46 
write.csv(go_result1,file = "01.down_GO.csv",row.names = T)
go_result <- go_result1

#提取各组数据需要展示的数量#
display_number = c(5, 5, 5)  ##这三个数字分别代表选取的BP、CC、MF的数量
BP <- go_result[which(go_result$ONTOLOGY=='BP'),]
CC <- go_result[which(go_result$ONTOLOGY=='CC'),]
MF <- go_result[which(go_result$ONTOLOGY=='MF'),]
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]

#将提取的各组数据进行整合
go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),  #指定ego_result_BP、ego_result_CC、ego_result_MFID为ID                        
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  Count=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count), #指定ego_result_BP、ego_result_CC、ego_result_MF的Count为GeneNumber
  type=factor(c(rep("Biological Process", display_number[1]), #设置biological process、cellular component、molecular function 的展示顺序
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))
##设置GO term名字长短，过长则设置相应长度
for(i in 1:nrow(go_enrich)){
  description_splite=strsplit(go_enrich$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")  #选择前五个单词GO term名字
  go_enrich$Description[i]=description_collapse
  go_enrich$Description=gsub(pattern = "NA","",go_enrich$Description)  #gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。gsub(“目标字符”, “替换字符”, 对象)
}
go_enrich <- na.omit(go_enrich)
#转成因子，防止重新排列
go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)
head(go_enrich)
p <- ggplot(go_enrich,
            aes(x=type_order,y=Count, fill=type)) +  #x、y轴定义；根据type填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("cadetblue1", "#FF6666", "#b2e7cb") ) + #柱状图填充颜色
  xlab("") + #x轴标签
  ylab("Counts") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() +
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 75,vjust = 1, hjust = 1 )) 
p
ggsave(filename = '02.down_GO.pdf',p,w=10,h=6)
ggsave(filename = '02.down_GO.png',p,w=10,h=6,units = 'in',dpi = 300)
dev.off()

## KEGG富集分析------
##KEGG气泡图---
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",   #大鼠的富集分析通路库
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "KEGG.csv",quote = F,row.names = F)  ##真正的输出结果3条
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
barplot(kk,showCategory = 10)

ggsave('02.KEGG_dot.png',kk_dot,width =8,height = 6)
ggsave('02.KEGG_dot.pdf',kk_dot,width =8,height = 6)

kk_bar<- barplot(kk,showCategory = 10)
kk_bar

ggsave('02.KEGG_bar.png',kk_bar,width =8,height = 6)
ggsave('02.KEGG_bar.pdf',kk_bar,width =8,height = 6)

# cnetplot(kk,categorySize="pvalue",foldChange = gene_transform$ENTREZID,colorEdge = T)
kk_net <- cnetplot(kk,foldChange = gene_transform$ENTREZID,circular=T,colorEdge=T)

ggsave('02.down_KEGG_net.png',kk_net,width =8,height = 6)
ggsave('02.down_KEGG_net.pdf',kk_net,width =8,height = 6)
