rm(list = ls()) # remove environment variables
library(tidyverse)
library(Seurat)
library(glmGamPoi)
library(DoubletFinder)

# build folders
figure_folder <- "/home/weihu/rb_revision_project/figures"
if(!dir.exists(figure_folder)){
  dir.create(figure_folder)
}

result_folder <- "/home/weihu/rb_revision_project/results"
if(!dir.exists(result_folder)){
  dir.create(result_folder)
}

# get sample list
dir_matrix <- "/home/weihu/rb_revision_project/data/matrix"
samples = list.files(dir_matrix)

######## Create Seurat objects ########

samList = lapply(samples,function(sp){
  folder=file.path(dir_matrix, sp)
  CreateSeuratObject(counts = Read10X(folder),
                     project = sp,
                     min.cells = 10,
                     min.features = 200)
})

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# manually adjust human cell cycle genes
s.genes <- cc.genes$s.genes
s.genes[s.genes == "MLF1IP"] <- "CENPU"
g2m.genes <- cc.genes$g2m.genes
g2m.genes[g2m.genes == "FAM64A"] <- "PIMREG"
g2m.genes[g2m.genes == "HN1"] <- "JPT1"

######## Filter cells ########

rt_summary <- data.frame()
names(samList)=samples
for(i in samples){
  samList[[i]] <- PercentageFeatureSet(object = samList[[i]],  pattern = "^MT-", col.name = "percent.mt") #rat mitochondrial genes start with Mt- 
  cell_before <- nrow(samList[[i]]@meta.data)
  
  nCount=samList[[i]]@meta.data$nCount_RNA
  nFeature=samList[[i]]@meta.data$nFeature_RNA
  mt=samList[[i]]@meta.data$percent.mt
  samList[[i]] <- AddMetaData(samList[[i]], log10(nFeature), col.name = "log.nFeature")
  samList[[i]] <- AddMetaData(samList[[i]], log10(nCount), col.name = "log.nCount_RNA")
  samList[[i]] <- subset(samList[[i]],
                         subset =
                           log.nFeature > median(log10(nFeature))-3*mad(log10(nFeature)) &
                           log.nCount_RNA > median(log10(nCount))-3*mad(log10(nCount)) &
                           percent.mt < median(mt) + 3*mad(mt))
  cell_after <- nrow(samList[[i]]@meta.data)
  
  # summary
  rt_summary_i <- data.frame(Sample = i, Before_Filter = cell_before, After_Filter = cell_after)
  rt_summary <- rbind(rt_summary, rt_summary_i)
  
  cat(i, "Before filter :", cell_before, "cells; ", "After filter :", cell_after,"cells\n")
  
}

######## Remove Doublets ########

# run DoubletFinder for individual sample
for(i in names(samList)){
  # Pre-process Seurat object
  pc.num=1:15
  samList[[i]] <- NormalizeData(samList[[i]])
  samList[[i]] <- FindVariableFeatures(samList[[i]], selection.method = "vst", nfeatures = 2000)
  samList[[i]] <- ScaleData(samList[[i]], verbose = FALSE)
  samList[[i]] <- RunPCA(samList[[i]], verbose=FALSE)
  samList[[i]] <- RunUMAP(samList[[i]], dims = pc.num)
  samList[[i]] <- FindNeighbors(samList[[i]], dims = pc.num) %>% FindClusters(resolution = 0.3)
  
  # pK Identification
  x.res.list <- paramSweep_v3(samList[[i]], PCs = pc.num, sct = FALSE)
  x.stats <- summarizeSweep(x.res.list, GT = FALSE)
  bcmvn <- find.pK(x.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  # Homotypic Doublet Proportion Estimate
  DoubletRate <- nrow(samList[[i]]@meta.data)*8*1e-6
  homotypic.prop <- modelHomotypic(samList[[i]]$seurat_clusters)
  nExp_poi <- round(DoubletRate*nrow(samList[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  samList[[i]] <- doubletFinder_v3(samList[[i]], PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  DF.name = colnames(samList[[i]]@meta.data)[grepl("DF.classification", colnames(samList[[i]]@meta.data))]
  rt_summary$Doublet_Rate[rt_summary$Sample == samList[[i]]@project.name] <- round(sum(samList[[i]]@meta.data[, colnames(samList[[i]]@meta.data) %in% DF.name] == "Doublet")*100/ncol(samList[[i]]), 2)
  
  samList[[i]] <- samList[[i]][, samList[[i]]@meta.data[, DF.name] == "Singlet"]
  rt_summary$After_Remove_Doublet[rt_summary$Sample == samList[[i]]@project.name] <- ncol(samList[[i]])
}

# calculate cell reduction rate
rt_summary$Cell_Reduction_Rate = round((rt_summary$Before_Filter - rt_summary$After_Remove_Doublet)*100/rt_summary$Before_Filter, 2)

samList <- lapply(X = samList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F) # Assign Cell-Cycle Scores
})
features <- SelectIntegrationFeatures(object.list = samList)
samList <- lapply(X = samList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


######## Perform integration ########

samList.anchors <- FindIntegrationAnchors(object.list = samList, anchor.features = features, reduction = "rpca")
sce.mer <- IntegrateData(anchorset = samList.anchors)

######## Quality control ########

setwd('/home/weihu/rb_revision_project/figures')
# QC1：the percentage of mitochondrial gene
pdf("QC1-VlnPlot.pdf",width = 8,height = 4.5)
VlnPlot(sce.mer, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0,
        ncol = 3)
dev.off()

# QC2：the correlation between gene numbers or mitochondrial gene and RNA numbers
plot1 <- FeatureScatter(sce.mer, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.mer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("QC2-FeatureScatter.pdf",width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

######## Perform clustering ########

DefaultAssay(sce.mer) <- "integrated"
sce.mer <- ScaleData(sce.mer, vars.to.regress = c("S.Score", "G2M.Score"), verbose=FALSE)
sce.mer <- RunPCA(sce.mer, npcs=30,verbose=FALSE)
sce.mer <- RunUMAP(sce.mer,reduction='pca',dims=1:30)
sce.mer <- FindNeighbors(sce.mer, reduction='pca',dims = 1:30)
sce.mer <- FindClusters(sce.mer, resolution = 0.8)

# export data
setwd('/home/weihu/rb_revision_project/results')
write.table(rt_summary, file = "filtering_summary.txt", sep = "\t", row.names = F, quote = F)
saveRDS(sce.mer,file='rb_integration_rmCycle_final.rds')

######## plot basic umap ########

setwd('/home/weihu/rb_revision_project/figures')

# orig.ident
pdf('rb_UMAP_orig.ident.pdf',width=7,height=5.5)
DimPlot(sce.mer,reduction='umap',group.by="orig.ident")
dev.off()

# cluster
pdf('rb_UMAP_all_clusters.pdf',width=7,height=5.5)
DimPlot(sce.mer,reduction='umap',label=T)
dev.off()

# cell cycle phase
pdf("rb_UMAP_Phase.pdf",width=7,height=5.5)
DimPlot(sce.mer, group.by="Phase", label=T, label.size=5, reduction='umap')
dev.off()

######## Find cluster markers ########
DefaultAssay(sce.mer) <- "RNA"
all_markers <- FindAllMarkers(sce.mer, only.pos = TRUE)

setwd('/home/weihu/rb_revision_project/results')
write.table(all_markers, "total_cluster_marker_genes.txt", sep="\t", quote = F, row.names = F)