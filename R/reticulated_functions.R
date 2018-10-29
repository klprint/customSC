#' Find clusters using HDBSCAN
#'
#' This function will use the scale.data slot of the Seurat object to find clusters using the HDBSCAN algorithm.
#' HDBSCAN is implemented in python and connected to R using the reticulate package.
#'
#' @param adata Seurat object with var.genes slot filled.
#' @return A character vector of class assignments. The class "-1" is used to label cells which could not be assigned to a cluster.
#' @details Make sure you have installed reticulate and enabled the right environment, before running this function.
#' The following python (3.6) packages need to be installed:
#'
#'   - umap
#'
#'   - hdbscan
hdbscan_clustering <- function(adata){

  if(is.null(adata)){
    stop("Find HVG genes first and save them in @var.genes slot of the Seurat object.")
  }
  mod.file <- system.file("py", package = "customSC")

  csc_py <- reticulate::import_from_path("hdbscan_clustering", mod.file)

  labels <- csc_py$hdbscan_cluster(t(adata@scale.data[adata@var.genes, ]))

  return(as.character(labels))
}
