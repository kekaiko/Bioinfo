setwd("D:\\Study\\Shanghaitech1\\Rotation4\\Bioinfo\\GEO_test\\GEO_test_20230616")
suppressPackageStartupMessages(library(GEOquery))  #屏蔽加载的红字


#把下载GEO数据集的表达矩阵这个过程包装成了函数  #下载GSE162903
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
  downGSE(studyID = "GSE162903")
}
# a <- read.table("GSE162903_series_matrix.txt.gz",comment.char = "!",header = T,sep = "\t",fill = T)
# # head(a)
exprSet <- read.csv("GSE162903_RNA_counts_filtered_matrix.csv.gz",header = T,row.names = 1)
pdata <- read.csv("GSE162903_metadata.csv",header = T,row.names = 1)

dim(exprSet)
