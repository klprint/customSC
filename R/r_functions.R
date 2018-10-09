# Reading 10x cellranger data using the output directory
get.10x.data <- function(path, cell.name.addition = NULL){

  barcodes <- read.table(file.path(path, "barcodes.tsv"), comment.char = "%")
  genes <- read.table(file.path(path, "genes.tsv"), comment.char = "%", sep = "\t")
  mtx <- read.table(file.path(path, "matrix.mtx"), comment.char = "%", sep = " ")[-1,]

  A <- Matrix::sparseMatrix(i = mtx[,1], j = mtx[,2], x = mtx[,3])

  row.names(A) <- genes[1:range(mtx[,1])[2],1]
  colnames(A) <- paste0(cell.name.addition, barcodes[1:range(mtx[,2])[2],1])

  return(A)
}

# Print the left top corner of a given dataframe or matrix
#' Corner
#'
#' Getting the top left corner of a dataframe, tibble or matrix.
#'
#' @param df Dataframe, matrix or tibble to be shown
corner <- function(df) df[1:min(c(5,nrow(df))),
                          1:min(c(5,ncol(df)))]


#' Generates a seurat object of 10x raw data
#'
#' Using the 10x output path with the genes, barcodes and matrix files.
#'
#' @param path Path to the 10x output
#' @param cell.name.addition Prefix for each cell barcode
#' @param min.genes Minimal number of genes which need to be detected per barcode to keep the barcode
#' @return Seurat object
make.seurat.object <- function(path,
                               cell.name.addition = NULL, min.genes = 100){
  A <- get.10x.data(path, cell.name.addition)

  gene.annot <- read.table(file.path(path, "genes.tsv"), comment.char = "%", sep = "\t")
  colnames(gene.annot) <- c("GeneID", "GeneName")

  A <- Seurat::CreateSeuratObject(raw.data = A, min.genes = min.genes)

  A@misc <- list(geneAnnotation = gene.annot)

  return(A)
}


# Produces a meta data table for each cell
# Can be used for the novel nuclei droplet identification
generate.frac.intron.meta <- function(sobj.ft, sobj.exon){
  ft <- sobj.ft@raw.data
  ex <- sobj.exon@raw.data

  keep.cells <- intersect(colnames(ft), colnames(ex))

  cat("Number of cells in both sets:",length(keep.cells), "\n")
  ex <- ex[,keep.cells]
  ft <- ft[,keep.cells]

  ex.sum <- as.data.frame(Matrix::colSums(ex))
  colnames(ex.sum) = "exonic"

  ft.sum <- as.data.frame(Matrix::colSums(ft))
  colnames(ft.sum) <- "fullTranscript"

  cells.meta <- merge(ft.sum, ex.sum, by = 0)
  colnames(cells.meta)[1] <- "CellID"

  cells.meta <- filter(cells.meta, fullTranscript > 0)
  cells.meta <- cells.meta %>% arrange(desc(fullTranscript)) %>%
    mutate(CellRank = 1:nrow(cells.meta)) %>%
    mutate(fracIntronic = 1 - (exonic / fullTranscript))

  sobj.exon@misc$cells.meta <- cells.meta

  meta <- sobj.exon@meta.data
  meta <- merge(meta, cells.meta, by.x = 0, by.y = "CellID", all.x = T)
  row.names(meta) = meta$Row.names
  meta <- meta[,-1]

  sobj.exon@meta.data <- meta
  return(sobj.exon)
}

# Plot the novel nuclei droplet identification plot
frac.intron.plot <- function(exon=NULL, fulltranscript=NULL, frac.cut = NULL, rank.cut = NULL){
  if(!is.null(fulltranscript)){
    exon <- generate.frac.intron.meta(fulltranscript, exon)
  }

  p <- ggplot(exon@misc$cells.meta, aes(x = CellRank, y = fracIntronic)) +
    geom_point(size = 0.1) +
    theme_gray() +
    ylab("Fraction intronic UMI") +
    xlab("Full transcript UMI cell rank")

  if(!is.null(frac.cut)){
    p <- p + geom_hline(yintercept = frac.cut)
  }

  if(!is.null(rank.cut)){
    p <- p + geom_vline(xintercept = rank.cut)
  }

  return(p)

}

#' Filtering droplets
#'
#' Filters droplets based on their properties, calculated using generate.frac.intron.meta
#'
#' @param sobj.exon Seurat object, pretreated with generate.frac.intron.meta
filter.nuc.droplets <- function(sobj.exon, min.frac.intronic = 0.5, max.rank = NULL){
  meta.cell <- sobj.exon@misc$cells.meta

  meta.cell <- meta.cell %>% filter(fracIntronic >= min.frac.intronic)

  if(!is.null(max.rank)){
    meta.cell <- meta.cell %>% filter(CellRank <= max.rank)
  }

  meta.cell <- filter(meta.cell, CellID %in% rownames(sobj.exon@meta.data))

  return(SubsetData(sobj.exon, cells.use = meta.cell$CellID))
}

#' Classical kneepoint plot
#'
#' Generates a kneepoint plot of the seurat object.
make.knee.plot <- function(sobj, max_rank = 2e4){
  meta <- sobj@meta.data

  meta <- meta %>% arrange(desc(nUMI)) %>%
    mutate(cumSum = cumsum(nUMI)) %>%
    mutate(UMIrank = 1:nrow(meta)) %>%
    filter(UMIrank < max_rank)

  ggplot(meta, aes(x = UMIrank, y = cumSum)) +
    geom_point() +
    theme_gray()
}
