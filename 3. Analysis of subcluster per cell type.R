rm(list = ls()) # remove environment variables
library(tidyverse)
library(Seurat)
library(ggsci)

# build folders
Dir <- "~/rb_revision_2_project/"
figure_folder <- paste0(Dir, "figures/figure_3/3_subcluster_cell_type")
if(!dir.exists(figure_folder)){
  dir.create(figure_folder)
}
result_folder <- paste0(Dir, "results/3_subcluster_cell_type")
if(!dir.exists(result_folder)){
  dir.create(result_folder)
}

# import data
sce.mer <- readRDS("~/rb_revision_2_project/results/base/rb_integration_rmCycle_final.rds")
sce.mer <- subset(sce.mer, sample != "Normal") # remove normal samples
table(sce.mer$sample)
table(sce.mer$manual_annotation)

# extract cell types
Cell_types <- c("MKI67_photoreceptorness_decreased", "retinoma_like", "Cone_precursor_like")

for(Cell in Cell_types){
  sce_sub <- subset(sce.mer, manual_annotation == Cell)
  DefaultAssay(sce_sub) <- "integrated"
  
  # cluster cell types
  sce_sub <- RunPCA(sce_sub, npcs=30,verbose=FALSE)
  sce_sub <- RunUMAP(sce_sub,reduction="pca",dims=1:30)
  sce_sub <- FindNeighbors(sce_sub, reduction="pca",dims = 1:30)
  sce_sub <- FindClusters(sce_sub, resolution = 0.4)
  
  setwd(result_folder)
  saveRDS(sce_sub, file = paste0("rb_", Cell, ".rds"))
  
  # cluster markers
  DefaultAssay(sce_sub) <- "RNA"
  marker <- FindAllMarkers(sce_sub, only.pos = TRUE)
  write.table(marker, paste0(Cell, "_cluster_marker_genes.txt"), sep = "\t", quote = F, row.names = F)
  
  # dot plot
  marker_sig <- marker %>% 
    filter(abs(avg_log2FC) > 0.25, p_val_adj < 0.05)
  marker_sig <- marker_sig[!grepl("MT-", marker_sig$gene),]
  top_markers <- marker_sig %>% 
    group_by(cluster) %>% 
    top_n(n = 5, wt = avg_log2FC) # top cluster markers
  
  setwd(figure_folder)
  p <- DotPlot(object = sce_sub, features = unique(top_markers$gene)) +
    scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))
  pdf(paste0(Cell, "_cluster_markers.pdf"), height = 8, width = 14)
  print(p)
  dev.off()
  
  # UMAP cluster 
  p <- DimPlot(sce_sub, reduction = "umap", label=T)
  pdf(paste0("rb_UMAP_", Cell, "_clusters.pdf"), width = 7, height = 5.5)
  print(p)
  dev.off()
  
  # UMAP sample 
  p <- DimPlot(sce_sub, reduction = "umap", group.by = "sample", label=T)
  pdf(paste0("rb_UMAP_", Cell, "_sample.pdf"), width = 7, height = 5.5)
  print(p)
  dev.off()
  
  # boxplot: cell type ratio 
  # prepare data
  data <- sce_sub@meta.data %>% 
    group_by(seurat_clusters, sample) %>% 
    summarise(Number = n())
  data$sample <- factor(data$sample, levels = c("Intraocular", "Extraocular"))
  
  ggplot(data = data, aes(x = seurat_clusters, y = Number, fill = sample)) +
    geom_bar(stat="identity", position = "fill") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
    scale_fill_npg() +
    labs(x = "Cluster", y = "Proportions (%)", fill = "Sample") +
    theme_test() +
    theme(axis.text = element_text(colour = "black"))
  ggsave(paste0(Cell, "_distribution_cell_type_sample_barplot.pdf"), width = 7, height = 6)
}