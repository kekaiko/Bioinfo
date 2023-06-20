rm(list = ls())  
load(file = "step1output.Rdata")
load(file = "step2output.Rdata")
#输入数据：exp和Group
#Principal Component Analysis的基础知识
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials

# 1.PCA 图----
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot
ggsave(plot = pca_plot,filename = paste0(gse_number,"_PCA.png"))
save(pca_plot,file = "pca_plot.Rdata")

# 2.top 1000 sd 热图---- 
cg=names(tail(sort(apply(exp,1,sd)),1000))
n=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=Group)#注解
rownames(annotation_col)=colnames(n) 
# 用标准化的数据画热图，两种方法的比较：https://mp.weixin.qq.com/s/jW59ujbmsKcZ2_CM5qRuAg，效果是一样的
## 1.使用热图参数
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
        # color =colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),#默认配色
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-3,3,length.out = 100)#length.out = 100和color参数有关
) #breaks 参数解读在上面链接
dev.off()
## 2.自行标准化再画热图
n2 = t(scale(t(n)))
pheatmap(n2,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         breaks = seq(-3,3,length.out = 100)
)
dev.off()