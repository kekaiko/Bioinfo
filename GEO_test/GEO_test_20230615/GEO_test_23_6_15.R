suppressPackageStartupMessages(library(GEOquery))  #屏蔽加载的红字

#下载数据
if(F){
  gset <- getGEO("GSE42589",destdir = ".",
                 # AnnotGPL = F,
                 # getGPL = F,   #这两行用于避免下载注释文件，太大了
                 
                 
  )
}


# a <- read.table("GSE42589_series_matrix.txt.gz",comment.char = "!",header = T,sep = "\t",fill = T)
# head(a)
#使用GPL文件进行转化(没学会)  
# if(F){
#   GPL_test <- getGEO("GPL6244",destdir = ".")  #生成了一个对象
#   class(GPL_test)
#   str(GPL_test)
#   colnames(Table(GPL_test)) 
#   #[1] "ID"              "GB_LIST"         "SPOT_ID"         "seqname"         "RANGE_GB"        "RANGE_STRAND"    "RANGE_START"    
#   #[8] "RANGE_STOP"      "total_probes"    "gene_assignment" "mrna_assignment" "category"  
#   head(Table(GPL_test)[,c(1,9,10)])}



#把下载GEO数据集的表达矩阵这个过程包装成了函数
if(F){
  downGSE <- function(studyID = "GSE1009", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir)
  
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
}
downGSE(studyID = "GSE42589")
}


exprSet <- read.csv("GSE42589_exprSet.csv",header = T,row.names = 1)
pdata <- read.csv("GSE42589_metadata.csv",header = T,row.names = 1)


#载入GPL6244对应的R包注释
library("hugene10sttranscriptcluster.db")

ls("package:hugene10sttranscriptcluster.db") 
#转化芯片ID
ids=toTable(hugene10sttranscriptclusterSYMBOL)
save(ids,exprSet,pdata,file = 'input.Rdata')
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))


plot(table(sort(table(ids$symbol))))

table(rownames(exprSet) %in% ids$probe_id) #判断rownames(exprSet)和hugene10sttranscriptcluster.db包中的probe_id有多少重合
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,] #取rownames(exprSet)和hugene10sttranscriptcluster.db包中的probe_id重合的探针
dim(exprSet)

ids=ids[match(rownames(exprSet),ids$probe_id),] #重新排序，按照匹配上的排序
head(ids)
exprSet[1:5,1:5]
#转化成sampleid
if(T){
  tmp = by(exprSet,ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  dim(exprSet)
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  dim(exprSet)
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2] #将exprSet数据框的行名（基因探针ID）更新为ids数据框中匹配的probe_id列对应的第二列，即更新为基因的符号（symbol）标识。这一步的目的是将基因表达数据的行名更改为基因符号，以提供更直观的标识。
  exprSet[1:5,1:5]  
}


group_list <- pdata[,1]
# 分隔字符串并提取所需部分
group_list <- sapply(strsplit(group_list, " sample "), function(x) x[1])
#我们使用strsplit()函数来分隔字符串。strsplit()函数用于将一个字符串分割为多个子字符串，并返回一个列表。我们将group_list作为要分割的字符串，并将" sample "作为分隔符。这样，每个字符串都会被分割为两部分，即样本类型和样本编号。
group_list <- gsub("NSCL/P", "process", group_list)# 将"NSCL/P"替换为"process"

## ggplot2作图
if(F){

library(reshape2)
library("dplyr")
copy_exprSet <- exprSet

copy_exprSet <- tibble::rownames_to_column(copy_exprSet, "SampleID")
exprSet_L=melt(copy_exprSet)


colnames(exprSet_L)=c('sampleID','sample','value')
exprSet_L$group=rep(group_list,each=nrow(exprSet))
head(exprSet_L)
### ggplot2 

  library(ggplot2)
  p=ggplot(exprSet_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density() 
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p)
}

## hclust 
if(F){
  colnames(exprSet)=paste(group_list,1:13,sep='')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10)) 
  plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
  
}
## PCA 
if(F){
  # BiocManager::install('ggfortify')
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list 
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
  
  library("FactoMineR")#画主成分分析图需要加载这两个包
  library("factoextra") 
  df=as.data.frame(t(exprSet))
  dat.pca <- PCA(df, graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
}

## t.test
if(F){
  dat = exprSet
  group_list=as.factor(group_list)
  group1 = which(group_list == levels(group_list)[1])
  group2 = which(group_list == levels(group_list)[2])
  dat1 = dat[, group1]
  dat2 = dat[, group2]
  dat = cbind(dat1, dat2)
  pvals = apply(exprSet, 1, function(x){
    t.test(as.numeric(x)~group_list)$p.value
  })
  p.adj = p.adjust(pvals, method = "BH")
  avg_1 = rowMeans(dat1)
  avg_2 = rowMeans(dat2)
  log2FC = avg_2-avg_1
  DEG_t.test = cbind(avg_1, avg_2, log2FC, pvals, p.adj)
  DEG_t.test=DEG_t.test[order(DEG_t.test[,4]),]
  DEG_t.test=as.data.frame(DEG_t.test)
  head(DEG_t.test)
}



## 差异分析   limma包
if(T){
  dim(exprSet)
  library(limma)
  
  #简单做个 QC检验
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  boxplot(exprSet, col = cols,main="expression value",las=2) #通过这些boxplot可以看到各个芯片直接数据还算整齐，可以进行差异比较
  
  #制作分组矩阵
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  #制作差异比较矩阵
  #contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
  contrast.matrix<-makeContrasts("control-process",levels = design)  #确保control 在前？
  contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
  
  #进行limma 差异分析
  limma_deg <- function(exprSet,design,contrast.matrix)
  {
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf,adjust="BH")#adjust='BH'表示对p值进行Benjamini-Hochberg多重检验校正。
    #coef = 1: 提取第一个系数对应的差异表达信息。这通常用于比较不同处理条件之间的差异。
    #coef = 2: 提取第二个系数对应的差异表达信息。这通常用于比较不同处理条件之间的差异，或者用于比较两个基因组之间的差异。
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  nrDEG_results<- limma_deg(exprSet,design,contrast.matrix)
  #nrDEG<- nrDEG_results[order(nrDEG_results$logFC, decreasing = TRUE), ]  ## 按logFC列对数据框进行降序排序
  nrDEG<- nrDEG_results #不排序
  head(nrDEG)
}

# ## 做火山图 用EnhancedVolcano包
if(F){
  library(EnhancedVolcano)
  #The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6，
  EnhancedVolcano(nrDEG,
                  lab = rownames(nrDEG),
                  x = 'logFC',
                  y = 'P.Value')
}



## for volcano 
if(F){

  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable',
              ifelse( df$logFC >1,'up',
                      ifelse( df$logFC < -1,'down','stable') )
  )
  table(df$g)#此时获得的上调基因为20，下调基因为145
  df$symbol=rownames(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "symbol", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = rownames(head(head(nrDEG))),
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  
  
}

##### volcano plot3(比较好看)
if(F){
  ## volcano plot
  DEG=nrDEG
  logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )#在上述代码中，with(DEG, mean(abs(logFC)) + 2*sd(abs(logFC))) 的作用是在 DEG 数据框的环境中计算 mean(abs(logFC)) + 2*sd(abs(logFC))。也就是说，logFC 列的平均值和两倍标准差是在 DEG 数据框的上下文中计算的，而无需每次使用 DEG$ 来引用列名。这样可以简化代码并提高可读性。
  DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT') 
                         #第一个 ifelse 的条件是 DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff，即 DEG 的 P 值小于 0.05 并且 DEG 的对数折叠变化的绝对值大于 logFC_cutoff。如果这个条件为真，表示 DEG 是显著差异表达基因。
                         
                         #在第一个 ifelse 的替代分支中，又嵌套了一个 ifelse，其条件是 DEG$logFC > logFC_cutoff，即 DEG 的对数折叠变化大于 logFC_cutoff。如果这个条件为真，表示 DEG 是上调基因（'UP'），否则表示 DEG 是下调基因（'DOWN'）。
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
  )
  this_tile
  head(DEG)
  g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
  print(g)
}



## 热图：for heatmap 
if(F){ 
  
  x=nrDEG$logFC
  names(x)=rownames(nrDEG)
  cg=c(names(head(sort(x),25)),#老大的代码中是100、100，由于我获得上下调基因少，就改了 
       names(tail(sort(x),25)))
  library(pheatmap)
  pheatmap(exprSet[cg,],show_colnames =F,show_rownames = F)
  n=t(scale(t(exprSet[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,
           show_rownames = T,
           cluster_cols = F,
           annotation_col=ac)
  
  
}

## heatmap 简单版
if(F){
  library(pheatmap)
  choose_gene=head(rownames(nrDEG),25)
  choose_matrix=exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix)
}

#GO 和 KEGG 分析

options(tidyverse.quiet = TRUE)
library(tidyverse)
options(data.table.quiet = TRUE)
library(data.table)

## 定义差异基因标准
deg <- nrDEG %>% 
  filter(abs(logFC) > 1.1) %>% 
  filter(P.Value < 0.05)

# deg <- deg %>% 
#   column_to_rownames("gene_id")  把第一列gene_id转为rownames



## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=0
deg$g=ifelse(deg$P.Value > 0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC< -logFC_t,'DOWN','stable') ))


table(deg$g)
head(deg)
deg$symbol=rownames(deg)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#从symbol转化为entrezid
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)
#合并DEG 和 df  通过SYMBOL 和 symbol
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)   
#gene_up = DEG[DEG$g == 'UP', 'ENTREZID']：这行代码使用逻辑条件筛选出数据框 DEG 中满足条件 DEG$g == 'UP' 的行，并提取这些行中的 ENTREZID 列，将结果赋值给 gene_up
gene_up= DEG[DEG$g == 'UP','ENTREZID']  

gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

#GO富集分析
if(F){
  go.up <- enrichGO(gene = gene_up,
                    OrgDb = "org.Hs.eg.db", 
                    ont="all",
                    universe = gene_all,
                    pvalueCutoff = 0.05, ##根据具体情况调整
                    qvalueCutoff =0.05, ####根据具体情况调整
                    readable=F) 
  
  head(go.up)[,1:6]
  dotplot(go.up )
  #ggsave('go.up.dotplot.pdf')
  
  go.down <- enrichGO(gene = gene_down,
                      OrgDb = "org.Hs.eg.db", 
                      ont="all",
                      universe = gene_all,
                      pvalueCutoff = 0.05, ##根据具体情况调整
                      qvalueCutoff =0.05,##根据具体情况调整
                      readable=F) 
  head(go.down)[,1:6]
  dotplot(go.down )
  #ggsave('go.down.dotplot.pdf')
  
  go.diff <- enrichGO(gene = gene_diff,
                      OrgDb = "org.Hs.eg.db", 
                      ont="all",
                      universe = gene_all,
                      pvalueCutoff = 1, ##根据具体情况调整
                      qvalueCutoff = 1)  ##根据具体情况调整
  head(go.diff)[,1:6]
  dotplot(go.diff )
  #ggsave('go.diff.dotplot.pdf')
  
  go_diff_dt <- as.data.frame(go.diff)
  go_down_dt <- as.data.frame(go.down)
  go_up_dt <- as.data.frame(go.up)
  
  ## 上面用的是富集分析中的pvalue,条目太多了，这里p.adjust
  ## 注意函数New.functions.R也要将pvalue改为p.adjust
  down_go<-go_down_dt[go_down_dt$p.adjust<0.0001,]##根据具体情况调整
  down_go$group=-1
  up_go<-go_up_dt[go_up_dt$p.adjust<0.0001,]##根据具体情况调整
  up_go$group=1
  #加载封装好的函数
  source('GO_KEGG_GSEA_New_functions.R')
  
  ## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
  go_figure=go.kegg_plot(up.data=up_go,down.data=down_go)
  print(go_figure)
  
  ggsave(go_figure,filename = 'GO_up_down_test_230616.png')
  #这个GO分析图,左边的英文名称代表差异基因富集到的 GO Terms，X轴表示-log10(dat$p.adjust)乘上分组，
  # 如：橙色的GO terms（其X轴大于0），表示上调基因富集的条目,X轴越大p.adjust越显著；
  # 蓝色的GO terms（其X轴小于0），表示下调基因富集的条目，X轴越小p.adjust越显著。
}

#KEGG富集分析
if(F){
  kegg.up <- enrichKEGG(gene = gene_up,
                        keyType = "kegg",
                        organism = "hsa",
                    universe = gene_all,  #background genes. 
                    pvalueCutoff = 0.9, ##根据具体情况调整
                    qvalueCutoff =0.9, ####根据具体情况调整
                    use_internal_data=T
  ) 
  
  head(kegg.up)[,1:6]
  dotplot(kegg.up )
  #ggsave('kegg.up.dotplot.pdf')
  
  kegg.down <- enrichKEGG(gene = gene_down,
                          keyType = "kegg",
                          organism = "hsa",
                          universe = gene_all,  #background genes. 
                          pvalueCutoff = 0.9, ##根据具体情况调整
                          qvalueCutoff =0.9, ####根据具体情况调整
                          use_internal_data=T
  ) 
  head(kegg.down)[,1:6]
  dotplot(kegg.down )
  #ggsave('kegg.down.dotplot.pdf')
  
  kegg.diff <- enrichKEGG(gene = gene_diff,
                        keyType = "kegg",
                        organism = "hsa",
                        universe = gene_all,  #background genes. 
                        pvalueCutoff = 0.9, ##根据具体情况调整
                        qvalueCutoff =0.9, ####根据具体情况调整
                        use_internal_data=T
  )
  head(kegg.diff)[,1:6]
  dotplot(kegg.diff )
  #ggsave('kegg.diff.dotplot.pdf')
  
  kegg_diff_dt <- as.data.frame(kegg.diff)
  kegg_down_dt <- as.data.frame(kegg.down)
  kegg_up_dt <- as.data.frame(kegg.up)
  
  ## 上面用的是富集分析中的pvalue,条目太多了，这里p.adjust
  ## 注意函数New.functions.R也要将pvalue改为p.adjust
  down_kegg<-kegg_down_dt[kegg_down_dt$p.adjust<1,]##根据具体情况调整
  down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$p.adjust<1,]##根据具体情况调整
  up_kegg$group=1
  #加载封装好的函数
  source('GO_KEGG_GSEA_New_functions.R')
  
  ## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
  kegg_figure=go.kegg_plot(up.data=up_kegg,down.data=down_kegg)
  print(kegg_figure)
  
  ggsave(kegg_figure,filename = 'KEGG_up_down_test_230616.png')
  #这个kegg分析图,左边的英文名称代表差异基因富集到的 GO Terms，X轴表示-log10(dat$p.adjust)乘上分组，
  # 如：橙色的GO terms（其X轴大于0），表示上调基因富集的条目,X轴越大p.adjust越显著；
  # 蓝色的GO terms（其X轴小于0），表示下调基因富集的条目，X轴越小p.adjust越显著。
}

#对于GSEA
if(F){
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',

                    nPerm        = 1000,##根据具体情况调整
                    minGSSize    = 70,##根据具体情况调整
                    pvalueCutoff = 0.9, ##根据具体情况调整
                    verbose      = FALSE,
                    by="DOSE"
                    )
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$p.adjust<0.05 & kk_gse$enrichmentScore < 0,]
  down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$p.adjust<0.05 & kk_gse$enrichmentScore > 0,]
  up_kegg$group=1
  
  ## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
  g_kegg=go.kegg_plot(up.data=up_kegg,down.data=down_kegg)
  #g_kegg=go.kegg_plot(up_kegg,down_kegg) #还可以更直接些
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea_test_230616.pdf')
  
}




