setwd('./VIPER/single-cell-pipeline/')

source('functions/process-utils.R')
source('functions/cluster-functions.R')
source('functions/viper-utils.R')
library(ggplot2)
library(ggpubr)
library(viper)
library(pheatmap)
library(RColorBrewer)
library(MUDAN)
library(umap)

full_data <- read.csv("/home/nghlm/Gao/Code/data/venetoclaxdata_matrix.csv", row.names = 1, check.names = FALSE)

expr <- full_data[, 1:(ncol(full_data) - 5)]
raw.mat <-  t(expr)  

pheno <- data.frame(CellType = full_data[, ncol(full_data) - 1])
rownames(pheno) <- rownames(full_data) 


mt.genes <- read.table('mt-genes.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
hum.mt <- mt.genes$hum.ensg
QCPlots(raw.mat, hum.mt)
mt.mat <- MTFilter(raw.mat, hum.mt)
filt.mat <- QCTransform(mt.mat)
cpm.mat <- CPMTransform(filt.mat)

# saveRDS(cpm.mat, file = "custom/venetoclax-cpm.rds")
cpm.mat <- readRDS("custom/venetoclax-cpm.rds")
rank.mat <- RankTransform(cpm.mat)

ARACNeTable(cpm.mat, 'custom/venetoclax-cpm.tsv')



#### 已有的两组标签（比如“对照 vs 处理”）直接在 PISCES 中比较蛋白活性差异，跳过聚类和meta-cell步骤，直接聚焦在：利用已有分组跑 ARACNe 网络,用 VIPER 计算蛋白活性 比较两组的蛋白活性差异（找 Master Regulators）

RegProcess('custom/network.txt', cpm.mat, out.dir = 'custom/', out.name = 'mydata-net-')
net <- readRDS('custom/mydata-net-pruned.rds')

pAct <- viper(rank.mat, net, method ='scale')

group <- pheno[ colnames(pAct), "CellType" ]

pAct <- Ensemble2GeneName(pAct) # 推荐转基因名便于解读

dim(pAct)
print(table(group)) 
print(all(names(group) == colnames(pAct))) # 检查 group 和 pAct 对齐

MRs <- BTTestMRs(pAct, group)
summary(MRs)



library(pheatmap)
library(RColorBrewer)
library(matrixStats)  # 用于 rowSds

# 1) 取出 20 个 MR，先裁剪极端值，避免颜色被几对极端数挤压
genes <- MR_UnWrap(MRs, top = 10)
pAct.sub <- pAct[genes, ]
pAct.sub[pAct.sub > 5]  <- 5   # 可根据需要调整阈值
pAct.sub[pAct.sub < -5] <- -5

# 2) 丢掉方差为 0 的基因，避免 scale="row" 产生全 NA
var_rows <- rowSds(as.matrix(pAct.sub)) > 0
pAct.sub  <- pAct.sub[var_rows, ]

# 3) 建注释、配颜色
ann <- data.frame(Group = factor(group, levels = c("resistant", "sensitive")))
rownames(ann) <- names(group)
ann_colors <- list(
  Group = c(resistant = "#E41A1C",  # 红
            sensitive = "#377EB8")  # 蓝
)

# 4) 画图（不再用 scale="row"；如果仍想行标准化，确保前面方差>0 即可）
pheatmap(
  pAct.sub,
  annotation_col  = ann,
  annotation_colors = ann_colors,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  main            = "Group Difference Master Regulators",
  color           = colorRampPalette(brewer.pal(11, "RdBu"))(100)
)





#### 聚类和meta-cell步骤之后在 PISCES 中比较蛋白活性差异

####  此时运行ARACNe算法 ####


RegProcess('custom/network.txt', cpm.mat, out.dir = 'custom/', out.name = 'venetoclax_r1-net-')

r1.net <- readRDS('custom/venetoclax_r1-net-pruned.rds')
r1.pAct <- viper(rank.mat, r1.net, method = 'none')

r1.viperDist <- as.dist(viperSimilarity(r1.pAct))
r1.clusts <- PamKRange(r1.viperDist, kmin = 2, kmax = 10)
r1.clustSil <- SilScoreEval(r1.clusts, r1.viperDist)
plot.dat <- data.frame('k' = 2:10, 'Silhouette.Scores' = r1.clustSil)
ggplot(plot.dat, aes(x = k, y = Silhouette.Scores)) + geom_point() + geom_line() +
  ggtitle('1.1 Clustering Silhouette Scores') + theme_bw()



#### 修改函数使其适应小样本 ####
MakeCMfA <- function(dat.mat, dist.mat, numNeighbors = 5, clustering, subSize = 200,
                     out.dir, out.name = '', sizeThresh = 0) {

  if (is.list(clustering) && !is.null(clustering$clustering)) {
    cluster_vec <- clustering$clustering
  } else {
    cluster_vec <- clustering
  }

  clust.mats <- ClusterMatrices(dat.mat, list(clustering = cluster_vec), sizeThresh = sizeThresh)

  k <- length(clust.mats)
  if (k == 0) {
    stop("No clusters met the size threshold. Try setting sizeThresh = 0.")
  }

  meta.mats <- list()
  for (i in 1:k) {
    mat <- clust.mats[[i]]
    meta.mat <- MetaCells(mat, dist.mat, numNeighbors, subSize)
    meta.mat <- CPMTransform(meta.mat)
    file.name <- paste(out.dir, out.name, '_clust-', i, '-metaCells.tsv', sep = '')
    ARACNeTable(meta.mat, file.name, subset = FALSE)
    meta.mats[[i]] <- meta.mat
  }
  return(meta.mats)
}

r1.clustMats <- MakeCMfA(filt.mat, r1.viperDist,
                         clustering = r1.clusts$k3,
                         out.dir = 'custom/',
                         out.name = 'venetoclax-r1-clusts',
                         subSize = 200) 


custom.clustering <- pheno$CellType
names(custom.clustering) <- rownames(pheno)  # 设置命名
all(names(custom.clustering) %in% colnames(filt.mat))  # 应该为 TRUE
# 1 resistant 2 sensitive
r1.clustMats <- MakeCMfA(filt.mat,
                         r1.viperDist,
                         clustering = custom.clustering,
                         out.dir = 'custom/',
                         out.name = 'venetoclax-r1-clusts',
                         subSize = 200)

####  对各个簇运行ARACNe算法 ####
####  删除arcane算法生成network文件的列名  ####
c1.net <- RegProcess('custom/r2-nets/venetoclax-r2-c1_finalNet.txt', r1.clustMats[[1]], 
                     out.dir = 'custom/r2-nets/', out.name = 'venetoclax-r2-c1_')
c2.net <- RegProcess('custom/r2-nets/venetoclax-r2-c2_finalNet.txt', r1.clustMats[[2]], 
                     out.dir = 'custom/r2-nets/', out.name = 'venetoclax-r2-c2_')
# c3.net <- RegProcess('custom/r2-nets/venetoclax-r2-c3_finalNet.tsv', r1.clustMats[[3]], 
                    #  out.dir = 'custom/r2-nets/', out.name = 'venetoclax-r2-c3_')

                     
# load in networks
c1.net <- readRDS('custom/r2-nets/venetoclax-r2-c1_pruned.rds')
c2.net <- readRDS('custom/r2-nets/venetoclax-r2-c2_pruned.rds')
# c3.net <- readRDS('tutorial/pbmc-r2-c3_pruned.rds')


# infer protein activity
r2.pAct <- viper(rank.mat, list('c1' = c1.net, 'c2' = c2.net), method = 'none')
# r2.pAct <- viper(rank.mat, list('c1' = c1.net, 'c2' = c2.net, 'c3' = c3.net), method = 'none')

r2.cbcMRs <- CBCMRs(r2.pAct) # identify the most representative proteins
r2.pAct.cbc <- r2.pAct[ r2.cbcMRs ,] # filter the protein activity matrix
r2.louvain <- LouvainClust(r2.pAct.cbc) # perform clustering analysis 
#默认参数无法分群,所以限制代表性蛋白数量（比如只保留 top 100）
length(r2.cbcMRs)
# 计算每个蛋白活性的标准差来
std.dev <- apply(r2.pAct[r2.cbcMRs, ], 1, sd)
top.50 <- names(sort(std.dev, decreasing = TRUE))[1:230]

r2.pAct.cbc <- r2.pAct[top.50, ]

r2.louvain <- LouvainClust(r2.pAct.cbc)
table(r2.louvain) #此时可以分群

### 3.1 Identifying Differentially activated Proteins

r2.pAct <- Ensemble2GeneName(r2.pAct)
r2.MRs <- BTTestMRs(r2.pAct, r2.louvain)

r2.pAct <- Ensemble2GeneName(r2.pAct)

ClusterHeatmap(r2.pAct[ MR_UnWrap(r2.MRs, top = 10) , ], clust = r2.louvain, plotTitle = 'Louvain Clustering: Differentially Activated Proteins')


top.mrs <- MR_UnWrap(r2.MRs, top = 10)
length(top.mrs)
print(top.mrs)







markers <- c('CCR7', 'CD8A', 'MS4A1', 'PPBP', 'MS4A7', 'IL7R', 'CD14', 'FCER1A', 'FCGR3A')
MarkerGrid(r2.cbcUMAP, r2.louvain, r2.pAct, markers, 'PBMC Marker Activity')

ClusterHeatmap(r2.pAct[ MR_UnWrap(r2.MRs, top = 10) , ], clust = r2.louvain, plotTitle = 'Louvain Clustering: Differentially Activated Proteins')

markers <- c('CCR7', 'CD8A', 'MS4A1', 'PPBP', 'MS4A7', 'IL7R', 'CD14', 'FCER1A', 'FCGR3A')
MarkerGrid(r2.cbcUMAP, r2.louvain, r2.pAct, markers, 'PBMC Marker Activity')