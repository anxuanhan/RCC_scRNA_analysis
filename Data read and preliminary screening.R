#数据处理
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)

sam.name <- "mergedata"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#读入数据，行是基因，列是细胞
control1 <- Read10X("./control1/") #36601行 6794880列
control2 <- Read10X("./control2/")
control3 <- Read10X("./control3/")

colnames(control1) <- paste(colnames(control1),"ctrl1",sep = "_")
colnames(control2) <- paste(colnames(control2),"ctrl2",sep = "_")
colnames(control3) <- paste(colnames(control3),"ctrl3",sep = "_")

ctrl.data <- cbind(control1,control2,control3)

# 创建单细胞分析对象
ctrl.aggregate <- CreateSeuratObject(
  ctrl.data,
  project = "mergedata",
  min.cells = 10, #滤掉在少于 10 个细胞中表达的基因
  min.features = 200, #过滤掉表达少于 200 基因的细胞。
  names.field = 2,
  names.delim = "_"
  
)

save(ctrl.aggregate,file = paste0("./",sam.name,"/",
                                  sam.name,"_raw_ctrlSeuratObject.Rdata"))

patient1 <- Read10X("./patient1/")
patient2 <- Read10X("./patient2/")
patient3 <- Read10X("./patient3/")

colnames(patient1) <- paste(colnames(patient1),"pt1",sep = "_")
colnames(patient2) <- paste(colnames(patient2),"pt2",sep = "_")
colnames(patient3) <- paste(colnames(patient3),"pt3",sep = "_")

tumor.data <- cbind(patient1,patient2,patient3)


# 创建单细胞分析对象
tumor.aggregate <- CreateSeuratObject(
  tumor.data,
  project = "mergedata",
  min.cells = 10, #滤掉在少于 10 个细胞中表达的基因
  min.features = 200, #过滤掉表达少于 200 基因的细胞。
  names.field = 2,
  names.delim = "_"
  
)

save(tumor.aggregate,file = paste0("./",sam.name,"/",
                                  sam.name,"_raw_tumorSeuratObject.Rdata"))
