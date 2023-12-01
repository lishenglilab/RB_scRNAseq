rm(list = ls()) # remove environment variables
library(tidyverse)
library(monocle3)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# extract cell types
Cell_types <- c("Cone_precursor_like", "MKI67_photoreceptorness_decreased", "retinoma_like")

for(Cell in Cell_types){
  
  # build folders
  Dir <- "~/rb_revision_2_project/"
  figure_folder <- paste0(Dir, "figures/figure_4/4_monocle3_", Cell)
  if(!dir.exists(figure_folder)){
    dir.create(figure_folder)
  }
  result_folder <- paste0(Dir, "results/4_monocle3_", Cell)
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  
  # import data
  sce.mer <- readRDS(paste0(Dir, "results/3_subcluster_cell_type/", "rb_", Cell, ".rds"))
  
  # Store data in a cell_data_set object
  matrix <- GetAssayData(sce.mer, assay = "RNA", slot = "data")
  gene_annotation <- data.frame(gene_short_name = rownames(matrix))
  rownames(gene_annotation) <- rownames(matrix)
  
  cds <- new_cell_data_set(matrix,
                           cell_metadata = sce.mer@meta.data,
                           gene_metadata = gene_annotation)
  
  # Normalize and pre-process the data
  cds <- preprocess_cds(cds, num_dim = 100, norm_method = "none")
  # cds <- align_cds(cds, alignment_group = "orig.ident")
  
  # Reduce dimensionality and visualize the results
  cds <- reduce_dimension(cds)
  
  # replace dimensionality data
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(sce.mer, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  # Cluster your cells
  cds <- cluster_cells(cds)
  plot_cells(cds, color_cells_by = "partition")
  
  # Learn the trajectory graph
  cds <- learn_graph(cds)
  plot_cells(cds,
             color_cells_by = "seurat_clusters",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  # Order the cells in pseudotime
  # a helper function to identify the root principal points:
  
  samples <- c("Intraocular")
  roots <- c()
  for( i in samples){
    get_earliest_principal_node <- function(cds, sample=i){
      cell_ids <- which(colData(cds)[, "sample"] == sample)
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
    roots <- c(roots, get_earliest_principal_node(cds))
  }
  
  cds <- order_cells(cds, root_pr_nodes=roots)
  cds@colData@listData$sample <- factor(cds@colData@listData$sample, levels = c("Intraocular", "Extraocular"))
  
  # plotting
  setwd(figure_folder)
  group_number <- length(unique(cds@colData@listData$sample))
  pdf("trajectory_pseudotime.pdf", width = 6, height = 4.5)
  p1 <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=3)
  print(p1)
  dev.off()
  
  pdf("trajectory_each_sample.pdf", width = 4* group_number, height = 4)
  p2 <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=3)  +
    facet_wrap(~sample, ncol = 3)
  print(p2)
  dev.off()
  
  invassive_gene <- data.table::fread("/home/weihu/rb_revision_2_project/data/侵袭转移相关基因list.txt", header = F, data.table = F, check.names = F)
  
  cds_subset <- cds[rowData(cds)$gene_short_name %in% invassive_gene$V1,]
  partition_summary <- data.frame(table(partitions(cds_subset)))
  cds_subset <- cds_subset[, partitions(cds_subset) %in% partition_summary$Var1[partition_summary$Freq == max(partition_summary$Freq)]]
  
  gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime", expression_family = "negbinomial")
  fit_coefs <- coefficient_table(gene_fits)
  sig_time_terms <- fit_coefs %>% filter(term == "pseudotime", q_value < 0.001) %>%
    select(gene_short_name, term, q_value, estimate)
  
  out <- fit_coefs %>% filter(term == "pseudotime") %>%
    select(gene_short_name, term, q_value, estimate)
  
  setwd(result_folder)
  write.table(out, "differential_invasive_gene_pseudotime.txt" ,row.names = F, sep="\t", quote = F)
  saveRDS(cds, paste0(Cell,"_monocle3.rds"))
  
  # plot selected differential results
  genes <- sig_time_terms$gene_short_name
  
  pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
  
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes
  
  #K means with 6 groups
  setwd(figure_folder)
  pdf("pseudotime_heatma_diff.pdf", width=9, height=8)
  htkm <- Heatmap(
    pt.matrix,
    name                         = "z-score",
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = TRUE,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    km = 4,
    row_title_rot                = 0,
    cluster_rows                 = TRUE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE)
  print(htkm)
  dev.off()
  
}
