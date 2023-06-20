# Group(实验分组)和ids(探针注释)
rm(list = ls())  
load(file = "step1output.Rdata")
library(stringr)
# 1.Group-----
# 第一类，有现成的可以用来分组的列
if(F) {Group = pd$`disease state:ch1` 
}
#第二类，自己生成（不建议用）
if(F){
  Group=c(rep("RA",times=13),
          rep("control",times=9))
}
rep(c("RA","control"),times = c(13,9))

#第三类，匹配关键词，自行分类
Group=ifelse(str_detect(pd$title,"control"),"control","patient")

#设置参考水平，指定levels，对照组在前，处理组在后（很重要）
Group = factor(Group,
               levels = c("control","patient"))
Group

# 注意levels与因子内容必须对应一致
# Group = pd$`disease state:ch1`
# Group = factor(Group,
#                 levels = c("healthy control","rheumatoid arthritis"))

#2.ids-----
#方法1 BioconductorR包(最常用)
gpl_number 
#http://www.bio-info-trainee.com/1399.html
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)

# 方法2 读取GPL平台的soft文件，按列取子集（还没用过）
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
if(F){
  #注：soft文件列名不统一，活学活用，有的GPL平台没有提供注释，如GPL16956
  a = getGEO(gpl_number,destdir = ".")
  b = a@dataTable@table
  colnames(b)
  ids2 = b[,c("ID","Gene Symbol")]
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]
}

# 方法3 官网下载,文件读取
if(F){
##http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
  }
# 方法4 自主注释 
if(F){
  #https://mp.weixin.qq.com/s/mrtjpN8yDKUdCSvSUuUwcA
}

save(exp,Group,ids,gse_number,file = "step2output.Rdata")
