rm(list = ls())
options(stringsAsFactors = F)
library(stringr)

data<-read.table("D:\\Study\\Shanghaitech1\\Rotation4\\RNASeq\\RNASeq_Matrix_Summary\\brainE13.5\\brainE13.5_230509_gene-count-matrix.txt",
                 header = T,
                 row.names = "gene") 
colnames_data<- colnames(data)   #原始数据的列名
new_colnames_data <- vector() #新的列名
#for循环用于提取新的行名   substr函数可以从一个字符串中提取子串，可以指定提取的起始位置和长度，用于从字符串中获取指定的子串内容。
#gregexpr("pairing.", colnames_data[i])[[1]][1] + 8是表示将从colnames_data[i]中匹配到的"pairing."的位置截取字符串，加上8个字符扩展截取范围，以便把"pairing."字符也截取进去。
#gregexpr("Aligned.", colnames_data[i])[[1]][1] - 1表示从字符串c[i]中查找“Aligned.”并获取第一次出现的位置，然后再减一。这意味着，从此位置到字符串的末尾，我们将获取所需的字符串，即从“pairing.”开始到“Aligned.”之前的字符串。
for (i in 1:length(colnames_data)) {
  new_colnames_data[i]<- sub("_Aligned.*", "", colnames_data[i]) #这里使用了sub函数，它可以通过正则表达式将字符串中的一部分替换成另一部分。正则表达式"_Aligned.*"会匹配下划线后面的所有字符，包括下划线本身。
}
colnames(data) <- new_colnames_data   #更改data 的列名为新列名

rownames_data <- rownames(data)
new_rownames_data <- vector() #新的行名
for (i in 1:length(rownames_data)) {
  new_rownames_data[i]<-sub("\\.\\d+$", "", rownames_data[i])   #去除ensembl id 的版本号
}
rownames(data) <- new_rownames_data
data<- data[rowSums(data) >= 5*length(colnames_data),]   #筛选

df_nomalize <- t(scale(t(data)))  #矩阵数据的标准化

#num_rows <- nrow(df_nomalize)
# num_cols <- ncol(df_nomalize)
# colnames(df_nomalize)<-c(1:num_cols)
head(data)
head(df_nomalize)

#对数据矩阵进行转置；
expt <- t(data[1:20,])
chemt <- t(df_nomalize[1:20,])
num_cols_chemt <- ncol(chemt)
colnames(chemt)<-c(1:num_cols_chemt)
#预览转置后的数据；
expt[1:9,1:6]
chemt[1:9,1:10]

#2.计算相关性
#安装psych包；
#install.packages("psych")
#载入R包；
library(psych)
#计算基因表达量之间的pearson相关性；
#ct1 <- corr.test(expt,chemt,method = "pearson")
ct1 <- corr.test(expt,chemt,method = "spearman")
#提取相关性系数矩阵；
r1 <- ct1$r
#提取pvalue值矩阵；
p1 <- round(ct1$p,5)
#预览转置后的相关性系数矩阵和pvalue矩阵；
r2 <- round(t(r1),3)
p2 <- t(p1)
r2[1:9,1:8]
p2[1:9,1:6]


#使用显著性星号标记进行替换；
p2[p2>=0 & p2 < 0.001] <- "***"
p2[p2>=0.001 & p2 < 0.01] <- "**"
p2[p2>=0.01 & p2 < 0.05] <- "*"
p2[p2>=0.05 & p2 <= 1] <- ""
#预览替换后的矩阵；
p2[1:9,1:9]



#使用ComplexHeatmap包绘制热图；
library(BiocManager)
#安装相关R包；
#install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

#创建颜色映射函数，数据范围与颜色相对应；
range(r2)
#[1] -0.808 0.883
#这里创建4种不同的配色方案；
col_fun1 = colorRamp2(c(range(r2)[1], 0, range(r2)[2]), c("#0f86a9", "white", "#FC8452"))
col_fun2 = colorRamp2(c(range(r2)[1], 0, range(r2)[2]), c("#A5CC26", "white", "#FF7BAC"))
col_fun3 = colorRamp2(c(range(r2)[1], 0, range(r2)[2]), c("#3FA9F5", "white", "#FF931E"))
col_fun4 = colorRamp2(c(-1, 0, 1), c("#ffa500", "white", "#B3A9EB"))
#查看颜色函数效果；
col_fun1(seq(-2, 2))
#[1] "#0F86A9FF" "#0F86A9FF" "#FFFFFFFF" "#FC8452FF" "#FC8452FF"

#热图格子大小设置；
cellwidth = 0.7
cellheight = 0.7
cn = dim(r2)[2]
rn = dim(r2)[1]
w=cellwidth*cn
h=cellheight*rn

#绘制相关性热图；
Heatmap(r2,name ="r", col = col_fun2,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "gray30",lwd = 1.4),
        row_dend_gp = gpar(col = "gray30",lwd = 1.4),
        #设置聚类gap数量小；
        row_split = 2, column_split = 2,
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)))




#绘制热图，显示相关性系数，保留两位小数；
png(filename = "r-data-heatmap.png", width = 800, height = 800, res = 100)

#png("r-data-heatmap.png")
Heatmap(r2,name ="r", col = col_fun1,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        border_gp = gpar(col = "#0f86a9",lty = 2,lwd = 1.2),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        row_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        #设置聚类gap数量和大小；
        row_split = 2, column_split = 2,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        #显示数值设置；
        #i,j对应数据矩阵的行和列索引；
        #x,y对应热图cell中心点坐标；
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", r2[i, j]), x, y,
                    gp = gpar(fontsize = 6))
        })
dev.off()


#绘制热图显示显著性星号标记；
Heatmap(r2,name ="r", col = col_fun2,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        border_gp = gpar(col = "#0f86a9",lty = 2,lwd = 1.2),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        row_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        #设置聚类gap数量和大小；
        row_split = 2, column_split = 2,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        #显示星号标记设置；
        #vjust垂直微调星号的位置；
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(p2[i, j], x, y, vjust = 0.7,
                    gp = gpar(fontsize = 13,col="white"))
        })




#只显示大于某一范围的数值；
Heatmap(r2,name ="r", col = col_fun3,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        border_gp = gpar(col = "#0f86a9",lty = 2,lwd = 1.2),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        row_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        #设置聚类gap数量和大小；
        row_split = 2, column_split = 2,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        #显示绝对值大于0.5的相关性系数；
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(abs(r2[i, j]) > 0.5){
            grid.text(sprintf("%.2f", r2[i, j]), x, y,
                      gp = gpar(fontsize = 6))}
        })
  

#使用加减号表示正负相关；
Heatmap(r2,name ="r", col = col_fun1,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        border_gp = gpar(col = "#0f86a9",lty = 2,lwd = 1.2),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        row_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        #设置聚类gap数量和大小；
        row_split = 2, column_split = 2,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        #使用加减号表示正负相关；
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(r2[i, j] > 0.5){
            grid.text("+", x, y,
                      gp = gpar(fontsize = 12,
                                fontface = "plain",
                                col="gray10"))
          }else if(r2[i, j] < -0.5){
            grid.text("-", x, y, vjust = 0.4,
                      gp = gpar(fontsize = 13,
                                fontface = "plain",
                                col="gray10"))
          }
        })

