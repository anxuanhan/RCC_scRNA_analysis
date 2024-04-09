library(dplyr)
library(Seurat)
library(patchwork)

data <- load("./mergedata/mergedata_raw_tumorSeuratObject.Rdata")

###步骤2：质控
###线粒体基因比例
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") 
#head(data@meta.data, 5)

###质控可视化
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
##图1：nFeature_RNA表达数，nCount_RNA细胞数
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
##图2：nFeature_RNA表达数，nCount_RNA细胞数，percent.mt线粒体比例三者间的关系

###根据图1和图2 去除表达小于200或大于7500或>20%线粒体基因
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA <4500 & percent.mt < 10)

###写出matrix_filtrate
#write.table(data@assays$RNA@counts,file="GBM filtrate qupici.txt",quote=F,sep="\t")

###步骤3：标准化
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
#data[["RNA"]]@data

###步骤4：寻找细胞间的高异质性基因
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)  #top2000差异基因
top10 <- head(VariableFeatures(data), 10)   #top10差异基因
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave("myplot.pdf", plot1 + plot2, width = 10, height = 10)
##图3：展示top2000差异基因中Top10高异质性基因

###步骤5：Scaling the data
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

###步骤6：PCA 分析
data <- RunPCA(data, features = VariableFeatures(object = data))
##表1：PCA中各PC所包含的基因

#print(data[["pca"]], dims = 1:5, nfeatures = 5) #只显示各PC中前5个基因

VizDimLoadings(data, dims = 1:12, reduction = "pca")
##图4：展示12个PC中特异性基因

DimPlot(data, reduction = "pca")
##图5：PCA结果展示

DimHeatmap(data, dims = 1:1, cells = 500, balanced = TRUE)
##图6：PCA中PC1所包含的基因

DimHeatmap(data, dims = 1:12, cells = 500, balanced = TRUE)
##图7：PCA中各PC所包含的基因

###步骤7：确定数据维度
#方法1 耗时4min
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
##图8：绘制PC曲线及p值，确定数据维度

#方法2
ElbowPlot(data)
##图9：根据拐点确定数据维度

###步骤8：根据上一步中的合适维度对细胞进行分群
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.3)#resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
#head(Idents(data), 4) #查看前4类

###步骤9：非线性分类(UMAP/tSNE)

#使用umap进行分群
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap", label = TRUE) 

##图10：umap分群图
#saveRDS(data, file = "E:/Desktop/table2 umap.rds")
##表2：umap分群结果保存
#readRDS(data, file = "E:/Desktop/GSCT/GBM/table2 umap.rds")
##读取上一步保存的umap分群结果

###写出umap cluster结果
#write.table(data@meta.data,file="GBM 10X umap.txt",quote=F,sep="\t")

# #使用tsne进行分群
#data <- RunTSNE(data, dims = 1:10)
#DimPlot(data, reduction = "tsne", label = TRUE)
##图11：tsne分群图
#saveRDS(data, file = "E:/Desktop/table3 tsne.rds")
##表3：tsne分群结果保存

###步骤10：差异基因
#cluster1.markers <- FindMarkers(data, ident.1 = 1, min.pct = 0.25) #查看cluster1的差异markers
#head(cluster1.markers, n = 5)
#cluster5.markers <- FindMarkers(data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)#查看cluster5与cluster0与3的差异markers
#head(cluster5.markers, n = 5)
#每个cluster与其他各cluster和所有clusters的差异，
#min.pct：在两组细胞中的任何一组中检测到的最小百分
#thresh.test：在两组细胞间以一定数量的差异表达（平均）
#max.cells.per.ident：通过降低每个类的采样值，提高计算速度
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#cluster1.markers <- FindMarkers(data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.table(data.markers,file="table 4 cluster markers.txt",quote=F,sep="\t")
##表4：分群特征性基因

VlnPlot(data, features = c('TRAC', 'IL7R', 'CD3D', 'CD3E', 'CD8A','CCR7','SELL','CD44','FAS','FASLG','CCL5','GZMK','CD8B' ))#CD8 T
VlnPlot(data, features = c('CD4','CD40LG','IL7R','CD3D', 'CD3E')) #2.CD4 T
VlnPlot(data, features = c('SELL','CXCL8', 'FCGR3B', 'MNDA','CXCR2'))  #3.Neutrophil
VlnPlot(data, features = c('MS4A1', 'CD40', 'IGHM', 'CD19', 'HLA-DRA','CD79A','MZB1','IGKC','JCHAIN')) #4.B cell
VlnPlot(data, features = c('TPBS2', 'KIT', 'ENPP3', 'SLC18A2','CPA3','TPSAB1'))#5.Mast cell
VlnPlot(data, features = c('IL3RA','THBD','CD1C','CD68','CD83','FCER1A','CST3','HLA-DQB2','BIRC3'))#6.Dendritic cell
VlnPlot(data, features = c('CD3D', 'TRAC', 'TRBC2', 'CD52', 'CD40LG','CCR7', 'CD3E', 'IL7R','CCR6', 'IL32'))#7.Th cell
VlnPlot(data, features = c('CD3D', 'TRAC', 'TRBC2', 'FOXP3','IL2RA', 'IL32', 'CTLA4', 'TNFRSF18')) #8.T regulatory
VlnPlot(data, features = c('CD56','CD16', 'NKG2D','NKP30','NKP44','NKP46','KLRD1'))#9.NK cell
VlnPlot(data, features = c('GNLY', 'CD3D', 'PRF1', 'CD8A', 'TRGC2','TRAC', 'KLRD1', 'IL7R', 'CD3E', 'FCGR3A'))#10.NKT cell
VlnPlot(data, features = c('KIT', 'TPSB2', 'FCER1A', 'SLC18A2','TPSAB1'))#11.Mast cell
VlnPlot(data, features = c('IGKC', 'JCHAIN', 'IGHA2', 'IGHG3', 'MZB1')) #12.Plasma cell
VlnPlot(data, features = c('CA1','HBB', 'HBE1', 'AHSP', 'ALAS2'))#13.Erythroblast
VlnPlot(data, features = c('CSF1R','CD68','CD163','CD14'))#14.Macrophage
VlnPlot(data, features = c('HLA-DR','LIN','CD44','CD11B'))#15.MDSC
VlnPlot(data, features = c('FAP','PDGFR','Vimentin','PDPN','CD70'))#16.CAF
VlnPlot(data, features = c('NDUFA4L2', 'CA9', 'SLC17A3'))#17.Renal cell carcinoma
VlnPlot(data, features = c('ITGAX', 'HLA-DR', 'Ly6c','CD1A','ITGB3','CD9','CCR6','TEK'))#18.Monocyte
##图12：展示某些基因在各组中的表达情况

#VlnPlot(data, features = c("CD3E", "CD4", "CD8A", "IL7R","PARP1","CD44"), slot = "counts", log = TRUE)
##图13：展示某些基因在各组中的表达情况（原始counts）

FeaturePlot(data, features = c("CCL2","CD74","CA9"), min.cutoff = "q10", max.cutoff = "q99")
FeaturePlot(data, features = c("CA9"), min.cutoff = "q10", max.cutoff = "q99")
##图14：展示某些基因在各组中的表达及分别情况，设置最低与最高显示值（也可用具体表达值，min.cutoff = 1, max.cutoff = 3）
#细胞周期 "MCM6","TOP2A","MKI67","PCNA"
#免疫标志物 "CD3E", "CD4", "CD8A", "IL7R","PARP1","CD28"

top10 <- data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(data, features = top10$gene) + NoLegend()
##图15：展示在各组中的特异性基因

RidgePlot(data, features =  c("CCL2","CD74"), ncol = 2)
##图16：特异性基因在各组中的表达

DotPlot(data, features =  c("CCL2","CD74")) + RotatedAxis()
##图17：特异性基因在各组中的表达


###步骤11：对分群进行标注
new.cluster.ids <- c("Proliferating CD8+ T cell", "CD4+T cell", "NK cell", "Macrophage", "T cell", "RCC", 
                     "Healthy cell", "Dendritic cell", "B cell", "CD8+T cell", "Mast cell", "Treg cell", 
                     "Plasma cell", "Endothelium","Neutrophil","16")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
##图18：Cluster增加标注
saveRDS(data, file = "./细胞分群new/Cluster Annotion")
##表5：存储上述结果

###步骤12：系统发育分析（Phylogenetic Analysis of Identity Classes）
data<-BuildClusterTree(data)
Tool(object = data, slot = 'BuildClusterTree')
PlotClusterTree(data)
##图19：Cluster Tree View

###步骤13：细胞周期分析
data <- CellCycleScoring(object = data,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)#+scale_color_npg()
##图20：各组细胞周期


###步骤14：提取某一cluster  并进行PCA分析
data  #查看总的细胞数
table(data@active.ident) #查看每个分组的细胞数
# Generate the table
active.ident.table <- table(data@active.ident)

# Write it to a CSV file
write.csv(active.ident.table, file = "./cellcluster_count.csv")

subdata<-subset(x = data,idents=c("RCC"))  #提取的Cluster名
subdata
###返回步骤6：subdata的PCA分析

##确定表达CCR中表达CCL2的细胞数量
ccl2_expression = subdata@assays$RNA$counts["CCL2", ]
num_cells_expressing_ccl2 = sum(ccl2_expression > 0)
num_cells_expressing_ccl2

##确定表达CCR中表达CC74的细胞数量
cd74_expression = subdata@assays$RNA$counts["CD74", ]
num_cells_expressing_cd74 = sum(cd74_expression > 0)
num_cells_expressing_cd74

##确定表达CCR中表达CA9的细胞数量
ca9_expression = subdata@assays$RNA$counts["CA9", ]
num_cells_expressing_ca9 = sum(ca9_expression > 0)
num_cells_expressing_ca9

##筛选同时表达CCL2和CD74的细胞
coexpression_index = (subdata@assays$RNA$counts["CCL2", ] > 0 & subdata@assays$RNA$counts["CD74", ] > 0)
coexpression_cells = subdata[, coexpression_index]
coexpression_cells

##筛选同时表达CCL2和CD74和CA9的细胞
coexpression_index = (subdata@assays$RNA$counts["CCL2", ] > 0 & subdata@assays$RNA$counts["CD74", ] > 0&subdata@assays$RNA$counts["CA9", ] > 0)
coexpression_cells = subdata[, coexpression_index]
coexpression_cells

# 获取CA9的表达数据
ca9_expression = subdata@assays$RNA$counts["CA9", ]

# 创建CCL2+CD74+CA9-细胞的子集
neg_index = (ccl2_expression > 0 & cd74_expression > 0 & ca9_expression == 0)
neg_cells = subdata[, neg_index]

# 创建CCL2+CD74+CA9+细胞的子集
pos_index = (ccl2_expression > 0 & cd74_expression > 0 & ca9_expression > 0)
pos_cells = subdata[, pos_index]

# 找出在neg_cells中表达，而在pos_cells中不表达的基因
expressed_in_neg = rowSums(neg_cells@assays$RNA$counts > 0) > 0
expressed_in_pos = rowSums(pos_cells@assays$RNA$counts > 0) > 0
genes_expressed_only_in_neg = names(expressed_in_neg[expressed_in_neg & !expressed_in_pos])

# 汇总neg_cells中每个基因的表达量
expression_sums = rowSums(neg_cells@assays$RNA$counts[genes_expressed_only_in_neg, ])

# 降序排序
sorted_genes = sort(expression_sums, decreasing = TRUE)

# 保存到CSV文件
write.csv(sorted_genes, file = "sorted_genes.csv")
