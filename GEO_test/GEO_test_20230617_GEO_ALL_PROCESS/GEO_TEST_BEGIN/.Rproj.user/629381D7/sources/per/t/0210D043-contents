rm(list = ls()) 
load(file = "step1output.Rdata")
load(file = "step4output.Rdata")
#1.火山图-----
library(dplyr)
library(ggplot2)
dat  = deg
#ggplot是可以赋值的
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

#加label，设置想要加上的label用if(T)
if(F){
  #自选基因
  for_label <- dat%>% 
    filter(symbol %in% c("HADHA","LRRFIP1")) 
}
if(F){
  #p值最小的10个
  for_label <- dat %>% head(10)
}
if(T) {
  #p值最小的前3下调和前3上调
  x1 = dat %>% 
    filter(change == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}
#实现在图上加label
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = volcano_plot,filename = paste0(gse_number,"_volcano.png"))

#2.差异基因热图-----

load(file = 'step2output.Rdata')
if(F){
  #全部差异基因
  cg = deg$probe_id[deg$change !="stable"]
  length(cg)
}else{
  #取前30上调和前30下调
  x=deg$logFC[deg$change !="stable"] 
  names(x)=deg$probe_id[deg$change !="stable"] 
  cg=names(c(head(sort(x),30),tail(sort(x),30)))
  length(cg)
}
n=exp[cg,]
dim(n)
#差异基因热图
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
heatmap_plot <- pheatmap(n,show_colnames =F,
                         show_rownames = F,
                         scale = "row",
                         cluster_cols = F, #处理组和对照组不聚类，将聚类去掉
                         annotation_col=annotation_col,
                         breaks = seq(-3,3,length.out = 100)
) 
heatmap_plot
ggsave(heatmap_plot,filename = paste0(gse_number,"_heatmap.png"))
load("pca_plot.Rdata")
library(patchwork)
library(ggplotify)
patchwork_plot <- ((pca_plot + volcano_plot) /as.ggplot(heatmap_plot))
patchwork_plot
ggsave(patchwork_plot,filename = paste0(gse_number,"_patchwork.png"))
