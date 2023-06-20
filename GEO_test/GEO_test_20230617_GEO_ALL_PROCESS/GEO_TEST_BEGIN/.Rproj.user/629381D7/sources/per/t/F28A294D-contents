rm(list = ls())  
load(file = 'step4output.Rdata')
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)

# 1.GO 富集分析-----

#(1)输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']

#(2)富集
#以下步骤耗时很长，设置了存在即跳过
if(!file.exists(paste0(gse_number,"_GO.Rdata"))){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,file = paste0(gse_number,"_GO.Rdata"))
}
load(paste0(gse_number,"_GO.Rdata"))

#(3)可视化
#条带图
barplot(ego)
#气泡图
dotplot(ego)
#个性化设置
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

#geneList 用于设置下面图的颜色
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(ego, showCategory = 3,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map,这个函数最近更新过，版本不同代码会不同
Biobase::package.version("enrichplot")

if(T){
  emapplot(pairwise_termsim(ego)) #新版本
}else{
  emapplot(ego)#老版本
}

#(4)展示通路关系 https://zhuanlan.zhihu.com/p/99789859
#goplot(ego)

#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList,showCategory = 8,label_format = 30)

# 2.KEGG pathway analysis-----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']

#（2）对上调/下调/所有差异基因进行富集分析
if(!file.exists(paste0(gse_number,"_KEGG.Rdata"))){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa')
  save(kk.diff,kk.down,kk.up,file = paste0(gse_number,"_KEGG.Rdata"))
}
load(paste0(gse_number,"_KEGG.Rdata"))

#(3)看看富集到了吗？https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$p.adjust<0.05)
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)

#(4)按照pvalue筛选通路
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)

#(5)可视化
source("GO_KEGG_GSEA_New_functions_20230617.R")
g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg
#g_kegg +scale_y_continuous(labels = c(10,5,0,5,10))，改坐标轴用
ggsave(g_kegg,filename = 'kegg_up_down.png')

# 3.gsea作kegg和GO富集分析----
## https://www.jianshu.com/p/c5b7b7dbf29b

#(1)查看示例数据
data(geneList, package="DOSE")
#(2)将我们的数据转换成示例数据的格式
geneList=deg$logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)
#(3)富集分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
#(4)可视化
g2 = kegg_plot(up_kegg,down_kegg)
g2
ggsave(g2,filename = 'gsea_kegg_up_down.png')