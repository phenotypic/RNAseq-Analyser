library(Seurat)
library(patchwork)
library(cowplot)
library(dplyr)

rk5.data <- Read10X(data.dir = "~/Downloads/RK5")
rk6.data <- Read10X(data.dir = "~/Downloads/RK6")

rk5 <- CreateSeuratObject(counts = rk5.data, min.cells = 3, min.features = 200, project = "rk5")
rk5$rk6 <- "rk5"
rk5[["percent.mt"]] <- PercentageFeatureSet(object = rk5, pattern = "^MT-")
rk5 <- subset(x=rk5, subset= nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <5)
rk5 <- NormalizeData(object = rk5, verbose = FALSE)
rk5 <- FindVariableFeatures(object = rk5, selection.method = "vst", nfeatures = 2000)

rk6 <- CreateSeuratObject(counts =rk6.data, min.cells = 3, min.features = 200, project = "rk6")
rk6$rk6 <- "rk6"
rk6[["percent.mt"]] <- PercentageFeatureSet(object = rk6, pattern = "^MT-")
rk6 <- subset(x=rk6, subset= nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <5)
rk6 <- NormalizeData(object = rk6, verbose = FALSE)
rk6 <- FindVariableFeatures(object = rk6, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(rk5, rk6), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

immune.combined[["clusterID"]] <- Idents(object = immune.combined)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "rk6")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "rk6")

ann_comb_markers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, logfc.threshold = 0.25)  
ann_comb_markers <- ann_comb_markers %>% dplyr::arrange(cluster, p_val_adj)
top5_comb <- ann_comb_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5_comb)

write.csv(ann_comb_markers, file = "combined_all_markers.csv", quote = FALSE, row.names = FALSE)
write.csv(top5_comb, file = "top_for_each_cluster.csv", quote = FALSE, row.names = FALSE)


FeaturePlot(immune.combined, features = c("aldh1a3", "sult2st3", "abcb5", "ebf3a.1", "lmx1ba", "pax7b", "cdh6", "CDK18", "si:dkey-7c18.24", "alk", "zgc:111983", "slc8a4b"), min.cutoff = "q9")
