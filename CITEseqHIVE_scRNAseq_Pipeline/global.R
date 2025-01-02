library(tools)
library(Seurat)
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)
library(scCustomize)

load_seurat_obj <- function(path){
  errors <- c()
  if (!tolower(tools::file_ext(path)) == 'rds'){
    errors <- c(errors, "Invalid rds file")
    return(errors)
  }
  tryCatch(
    {
      obj <- readRDS(path)
    },
    error = function(e){
      errors <- c(errors, "Invalid rds file")
      return(errors)
    }
  )
  if (!inherits(obj, "Seurat")){
    errors <- c(errors, "File is not a Seurat object")
  }
  return(obj)
}

# Create UMAPs for metadata
create_metadata_UMAP <- function(obj, col){
  if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mito")){
    col_df <- data.frame(obj@reductions$umap@cell.embeddings, data = obj@meta.data[,col])
    umap <- ggplot(data = col_df) +
      geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(data)), size = 0.01) +
      scale_colour_gradientn(colours = rainbow(7))
  } else if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
  } else {
    umap <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(umap)
}

# Compare UMAPs
compare.umaps <- function(obj, group.by = NULL, split.by = NULL) {
  valid_splits <- c("HIVE", "PID", "VISIT", "STIM", "VISIT_STIM", "group")
  
  # Check if split.by is valid
  if (!is.null(split.by) && split.by %in% valid_splits) {
    umaps <- DimPlot_scCustom(obj, 
                     reduction = "umap", 
                     label = TRUE, 
                     group.by = group.by,
                     split.by = split.by,
                     pt.size = 0.5,
                     num_columns = 2)
  } else {
    umaps <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Comparison parameter invalid"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(umaps)
}

# Create feature plots for all genes
create_feature_plot <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
  } else {
    FP <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
}

# Set list of antigen specific genes and proteins
AS.genes_RNA <- c("CD3G", "CD4", "CD8A", "CD69", "CD40LG", "TNFRSF9", "IL2RA", "TNFRSF4", "HLA-DRA")
AS.genes_ADT <- c("CD3-UCHT1", "CD4-RPA.T4", "CD8", "CD69", "CD154", "CD137", "CD25", "CD134", "HLA.DR")

# Create a function to add margins to plots
add_margins_and_subtitle <- function(plot, subtitle) {
  plot + 
    labs(subtitle = subtitle) + 
    theme(plot.margin = unit(c(1, 1, 1, 3), "cm"), 
          plot.subtitle = element_text(hjust = 0.5))
}

# Create violin plots for antigen-specific genes
create_AS.vlnplots_RNA <- function(obj_T, obj_M72_T, obj_UNSTIM_T, AS_gene, layer.gene) {
  DefaultAssay(obj_T) <- "RNA"
  DefaultAssay(obj_M72_T) <- "RNA"
  DefaultAssay(obj_UNSTIM_T) <- "RNA"
  
  if (AS_gene %in% AS.genes_RNA){
    vln_RNA <- VlnPlot(obj_T, 
                       features = AS_gene,
                       split.by = "VISIT",
                       group.by = "sctype_classification",
                       layer = layer.gene,
                       log = TRUE,
                       pt.size = 0,
                       combine = FALSE)
    vln_RNA <- lapply(vln_RNA, add_margins_and_subtitle, subtitle = "All RNA")
    
    vln_RNA_M72 <- VlnPlot(obj_M72_T, 
                           features = AS_gene,
                           split.by = "VISIT",
                           group.by = "sctype_classification",
                           layer = layer.gene,
                           log = TRUE,
                           pt.size = 0,
                           combine = FALSE)
    vln_RNA_M72 <- lapply(vln_RNA_M72, add_margins_and_subtitle, subtitle = "M72 RNA")
    
    vln_RNA_UNSTIM <- VlnPlot(obj_UNSTIM_T, 
                              features = AS_gene,
                              split.by = "VISIT",
                              group.by = "sctype_classification",
                              layer = layer.gene,
                              log = TRUE,
                              pt.size = 0,
                              combine = FALSE)
    vln_RNA_UNSTIM <- lapply(vln_RNA_UNSTIM, add_margins_and_subtitle, subtitle = "UNSTIM RNA")
    
    combined_vln_RNA <- wrap_plots(
      plotlist = c(vln_RNA, vln_RNA_M72, vln_RNA_UNSTIM), ncol = 1)
    
  } else {
    combined_vln_RNA <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(combined_vln_RNA)
}

# Create violin plots for antigen-specific proteins
create_AS.vlnplots_ADT <- function(obj_T, obj_M72_T, obj_UNSTIM_T, AS_protein, layer.protein) {
  DefaultAssay(obj_T) <- "ADT"
  DefaultAssay(obj_M72_T) <- "ADT"
  DefaultAssay(obj_UNSTIM_T) <- "ADT"
  
  if (AS_protein %in% AS.genes_ADT){
    vln_ADT <- VlnPlot(obj_T, 
                       features = AS_protein,
                       split.by = "VISIT",
                       group.by = "sctype_classification",
                       layer = layer.protein,
                       log = TRUE,
                       pt.size = 0,
                       combine = FALSE)
    vln_ADT <- lapply(vln_ADT, add_margins_and_subtitle, subtitle = "All ADT")
    
    vln_ADT_M72 <- VlnPlot(obj_M72_T, 
                           features = AS_protein,
                           split.by = "VISIT",
                           group.by = "sctype_classification",
                           layer = layer.protein,
                           log = TRUE,
                           pt.size = 0,
                           combine = FALSE)
    vln_ADT_M72 <- lapply(vln_ADT_M72, add_margins_and_subtitle, subtitle = "M72 ADT")
    
    vln_ADT_UNSTIM <- VlnPlot(obj_UNSTIM_T, 
                              features = AS_protein,
                              split.by = "VISIT",
                              group.by = "sctype_classification",
                              layer = layer.protein,
                              log = TRUE,
                              pt.size = 0,
                              combine = FALSE)
    vln_ADT_UNSTIM <- lapply(vln_ADT_UNSTIM, add_margins_and_subtitle, subtitle = "UNSTIM ADT")
    
    combined_vln_ADT <- wrap_plots(
      plotlist = c(vln_ADT, vln_ADT_M72, vln_ADT_UNSTIM), ncol = 1)
    
  } else {
    combined_vln_ADT <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(combined_vln_ADT)
}

# Create feature plots for antigen-specific genes
create_AS.featureplots_RNA <- function(obj_T, AS_gene) {
  DefaultAssay(obj_T) <- "RNA"
  
  if (AS_gene %in% AS.genes_RNA){
    feature_RNA <- FeaturePlot(obj_T, 
                               features = AS_gene,
                               label = TRUE,
                               split.by = "VISIT_STIM",
                               pt.size = 0,
                               combine = FALSE)
    combined_feature_RNA <- wrap_plots(feature_RNA, ncol = 2) 
    
  } else {
    feature_RNA <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(combined_feature_RNA)
}

# Create feature plots for antigen-specific proteins
create_AS.featureplots_ADT <- function(obj_T, AS_protein) {
  DefaultAssay(obj_T) <- "ADT"
  
  if (AS_protein %in% AS.genes_ADT){
    feature_ADT <- FeaturePlot(obj_T, 
                               features = AS_protein,
                               label = TRUE,
                               split.by = "VISIT_STIM",
                               pt.size = 0,
                               combine = FALSE)
    combined_feature_ADT <- wrap_plots(feature_ADT, ncol = 2) 
    
  } else {
    feature_ADT <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(combined_feature_ADT)
}
