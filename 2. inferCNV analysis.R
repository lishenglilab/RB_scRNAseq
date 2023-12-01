rm(list = ls()) # remove environment variables
library(tidyverse)
library(Seurat)
library(infercnv)
library(Hmisc)

# import data
sce.mer <- readRDS("~/rb_revision_2_project/results/base/rb_integration_rmCycle_final.rds")
sce.mer <- sce.mer <- subset(sce.mer, manual_annotation %nin% c("Müller_glia", "Microglia"))
table(sce.mer$sample, sce.mer$manual_annotation)

cells <- c("Cone_precursor_like", "retinoma_like", "MKI67_photoreceptorness_decreased")
for(i in cells){
  print(i)
  
  # extract cells
  sce.sub <- subset(sce.mer, manual_annotation %in% i | sample %in% "Normal")
  
  # build folders
  Dir <- "/home/weihu/rb_revision_2_project/results/1_inferCNV_tumor/"
  rt_folder <- paste0(Dir, i, "/")
  if(!dir.exists(rt_folder)){
    dir.create(rt_folder)
  }
  setwd(rt_folder)
  
  # prepare input data
  # annotation file
  cellAnnota <- data.frame(V1 = rownames(sce.sub@meta.data), V2 = sce.sub@meta.data$sample)
  write.table(cellAnnota, "cellAnnota.txt", col.names = F, row.names = F, sep = "\t", quote = F)
  
  # matrix file
  DefaultAssay(sce.sub) <- "RNA"
  exprMatrix <- as.matrix(GetAssayData(sce.sub, slot="counts"))
  write.table(exprMatrix, "exprMatrix.txt", sep = "\t", quote = F)

  # start inferCNV
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix= "exprMatrix.txt",
                                      annotations_file= "cellAnnota.txt",
                                      delim="\t",
                                      gene_order_file="/home/public/reference/inferCNV/gencode_V35_position.txt",
                                      ref_group_names = "Normal") 
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir= "./", 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=FALSE,
                               output_format="pdf")
}