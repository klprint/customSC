
# Define the class itself
corner <- function(df) df[1:min(5,nrow(df)), 1:min(5,ncol(df))]

createExperiment <- setClass("scExperiment",
                             slots = list(
                               counts = "dgCMatrix",
                               disp.counts = "dgCMatrix",
                               disp.exprs = "dgCMatrix",
                               description = "character",
                               sf = "numeric",
                               obs = "list",
                               disp.genes = "character",
                               var = "data.frame",
                               anscombe = "dgCMatrix",
                               dr = "list"
                             ))

## Define the initialization method
setMethod("initialize", "scExperiment",
          function(.Object,
                   counts = numeric(0),
                   description = character(""),
                   sf = numeric(0),
                   obs = list(),
                   var = data.frame(),
                   anscombe = numeric(0),
                   disp.counts = numeric(0),
                   disp.exprs = numeric(0),
                   dr = list()){
            .Object@counts <- Matrix::Matrix(counts, sparse = T)
            .Object@description <- description
            .Object@sf <- colSums(counts) / mean(colSums(counts))
            .Object@obs <- list("n_genes" = apply(counts > 0, 2, sum))
            .Object@var <- data.frame(gene = row.names(counts),
                                      n_cells = apply(counts > 0, 1, sum),
                                      stringsAsFactors = F)

            .Object
          })
## How is the method printed
setMethod("show", "scExperiment",
          function(object){
            cat(object@description,"\n")
            cat("\ngenes x cells\n")
            print(dim(object@counts))
            cat("Corner of the count matrix:\n\n")
            print(corner(object@counts))
          })

# Plotting a gene histogram
setGeneric(name = "plotGenesDist",
           def = function(experiment, n_bins = 30){
             standardGeneric("plotGenesDist")
           })

setMethod("plotGenesDist", "scExperiment",
          function(experiment, n_bins){
            hist(experiment@obs$n_genes, main="Number of genes per observation", breaks = n_bins)
          })

# Anscombe transformation of the data
setGeneric(name = "doAnscombeTransform",
           def = function(experiment, sf = NULL){
             standardGeneric("doAnscombeTransform")
           })
setMethod("doAnscombeTransform", "scExperiment",
          function(experiment, sf = NULL){
            ans <- apply(experiment@counts, 2, function(x){
              sqrt(x + 3/8)-sqrt(3/8)})

            if(is.null(sf)){
              ans <- ans / sqrt(experiment@sf)
            }else{
              ans <- ans / sqrt(sf)
            }

            experiment@anscombe <- Matrix::Matrix(ans, sparse = T)

            return(experiment)
          })

setMethod("doAnscombeTransform", "seurat",
          function(experiment, sf = NULL){
            if(is.null(sf)){
              cat("Calculate SF\n")
              sf <- experiment@meta.data$nUMI / mean(experiment@meta.data$nUMI)
            }

            cat("Transform data\n")
            # ans <- apply(experiment@raw.data, 2, function(x){
            #   sqrt(x + 3/8) - sqrt(3/8)
            # })
            ans <- NULL
            n_cells <-  ncol(experiment@data)

            # for(i in 1:ncol(experiment@raw.data)){
            #   z <- experiment@raw.data[,i]
            #   z <- sqrt(z + 3/8) - sqrt(3/8)
            #   ans <- cbind(ans, Matrix::Matrix(z, sparse = T))
            #   setTxtProgressBar(pb, i)
            # }
            cuts <- split(1:n_cells, cut(1:n_cells, quantile(1:n_cells, probs = seq(0,1,by=0.1)), include.lowest = T))
            pb <- pbmcapply::progressBar(min = 1, max = length(cuts))
            for(j in 1:length(cuts)){
              tmp <- experiment@data[,cuts[[j]]]
              ans <- cbind(ans, Matrix::Matrix(apply(tmp, 2, function(x) sqrt(x + 3/8) - sqrt(3/8)), sparse = T))
              setTxtProgressBar(pb, j)
            }
            close(pb)
            cat("Removing SF\n")
            print(dim(ans))
            ans <- ans / sqrt(sf)


            #experiment@scale.data <- Matrix::Matrix(ans, sparse = T)
            experiment@scale.data <- ans
            experiment@data <- ans
            experiment@meta.data$sf <- sf

            return(experiment)
          })

# Dispersed genes identification
setGeneric(name = "findDispGenes",
           def = function(experiment,
                          mean.cutoff = 0.25,
                          var.per.mean.cutoff = 0.95){
             standardGeneric("findDispGenes")
           })
setMethod("findDispGenes", "scExperiment",
          function(experiment,
                   mean.cutoff = 0.25,
                   var.per.mean.cutoff = 0.95){
            means <- apply(experiment@counts / experiment@sf, 1, mean)
            vars <- apply(experiment@counts / experiment@sf, 1, var)
            vars.per.means = vars / means

            disp.genes <- names(means[means > quantile(means, mean.cutoff)])
            disp.genes <- intersect(disp.genes, names(vars.per.means[vars.per.means > quantile(vars.per.means, var.per.mean.cutoff)]))

            experiment@disp.genes <- disp.genes
            experiment@var$is.dispersed <- names(means) %in% disp.genes
            experiment@disp.counts <- Matrix::Matrix(experiment@counts[disp.genes,], sparse = T)

            if(ncol(experiment@anscombe) == 0){
              experiment <- doAnscombeTransform(experiment)
            }
            experiment@disp.exprs <- experiment@anscombe[disp.genes, ]


            return(experiment)
          })

setMethod("findDispGenes", "seurat",
          function(experiment,
                   mean.cutoff = 0.25,
                   var.per.mean.cutoff = 0.95){
            require(ggplot2)

            if(!("sf" %in% colnames(experiment@meta.data))){
              experiment@meta.data$sf <- experiment@meta.data$nUMI / mean(experiment@meta.data$nUMI)
            }

            means <- apply(experiment@data / experiment@meta.data$sf, 1, mean)
            vars <- apply(experiment@data / experiment@meta.data$sf, 1, var)
            vars.per.means = vars / means

            disp.genes <- names(means[means > quantile(means, mean.cutoff, na.rm = T)])
            disp.genes <- intersect(disp.genes,
                                    names(vars.per.means[vars.per.means > quantile(vars.per.means,
                                                                                   var.per.mean.cutoff,
                                                                                   na.rm = T)]))

            experiment@var.genes <- disp.genes
            experiment@hvg.info <- data.frame(gene.mean = means,
                                              gene.var = vars,
                                              gene.dispersion = vars.per.means,
                                              gene.dispersion.scaled = vars.per.means,
                                              is.dispersed = names(means) %in% disp.genes)
            row.names(experiment@hvg.info) <- names(means)
            #experiment@var.genes <- names(means) %in% disp.genes
            #experiment@disp.counts <- Matrix::Matrix(experiment@counts[disp.genes,], sparse = T)

            # if(ncol(experiment@anscombe) == 0){
            #   experiment <- doAnscombeTransform(experiment)
            # }
            # experiment@disp.exprs <- experiment@anscombe[disp.genes, ]
            print(
            ggplot(experiment@hvg.info, aes(x = gene.mean, y = gene.dispersion, color = is.dispersed) ) +
              geom_point(size = 0.5, alpha = 0.5) +
              scale_y_log10() +
              scale_x_log10() +
              ylab("variance / mean") +
              xlab("mean") +
              geom_hline(yintercept = mean(1/experiment@meta.data$sf))
            )
            return(experiment)
          })

multiFindDispGenes <- function(multiseurat, mean.cutoff = .25, var.per.mean.cutoff = .95){
  disp.genes <- list()
  for(sample in levels(multiseurat@meta.data$orig.ident)){
    cat("Processig", sample, "\n")
    tmp.s <- SubsetData(multiseurat, subset.name = "orig.ident", accept.value = sample)
    tmp.s <- findDispGenes(tmp.s, mean.cutoff = mean.cutoff, var.per.mean.cutoff = var.per.mean.cutoff)
    if(is.logical(multiseurat@var.genes)){
      multiseurat@var.genes <- tmp.s@var.genes
    }else{
      multiseurat@var.genes <- intersect(multiseurat@var.genes, tmp.s@var.genes)
    }
  }
  return(multiseurat)
}

# Do PCA
setGeneric(name = "doPCA",
           def = function(experiment, do.scale = F){
             standardGeneric("doPCA")
           })

setMethod("doPCA", "scExperiment",
          function(experiment, do.scale = F){
            if(ncol(experiment@disp.exprs) == 0){
              experiment <- findDispGenes(experiment)
            }
            mat <- t(as.matrix(experiment@disp.exprs))
            if(do.scale){
              mat <- scale(mat)
            }
            pca <- prcomp(mat)
            experiment@dr$pca = pca

            return(experiment)
          })
setGeneric(name = "plotPCA",
           def = function(experiment){
             standardGeneric("plotPCA")
           })
setMethod("plotPCA", "scExperiment",
          function(experiment){
            pca = as.data.frame(experiment@dr$pca$x)
            ggplot(pca, aes(x = PC1, y = PC2)) +
              geom_point(size = 0.1)
          })

# Do UMAP
setGeneric(name = "doUMAP",
           function(experiment, do.plot = F, on.pca = F, n_pcs = 10, do.scale = F){
             standardGeneric("doUMAP")
           })

setMethod("doUMAP", "scExperiment",
          function(experiment, do.plot = F, on.pca = F, n_pcs = 10, do.scale = F){
            mat <- as.matrix(experiment@disp.exprs)

            if(nrow(mat) == 0){
              stop("No dispersed genes were identified!\nRun findDispGenes first!")
            }

            if(on.pca){
              mat <- experiment@dr$pca$x[,1:n_pcs]
            }else{
              mat <- t(mat)
            }

            if(do.scale & !on.pca){
              mat <- scale(mat)
            }

            umap.coords <- umap(mat)

            experiment@dr$umap <- umap.coords

            return(experiment)
          })

setGeneric(name = "plotUMAP",
           def = function(experiment, color = NULL){
             standardGeneric("plotUMAP")
           })
setMethod("plotUMAP", "scExperiment",
          function(experiment){
            df = as.data.frame(experiment@dr$umap$layout)
            ggplot(df, aes(x = V1, y = V2)) +
              geom_point(size = 0.1)
          })



setGeneric(name = "plotGene",
           def = function(experiment, gene, dr = "UMAP"){
             standardGeneric("plotGene")
           })
setMethod("plotGene", "scExperiment",
          function(experiment, gene, dr = "UMAP"){
            if(dr == "UMAP"){
              plotUMAP(experiment) +
                geom_point(aes(color = experiment@anscombe[gene,]), size = 0.1) +
                scale_color_continuous(name = gene, high = "red", low = "gray")
            }else if(dr == "pca"){
              plotPCA(experiment = experiment) +
                geom_point(aes(color = experiment@anscombe[gene,]), size = 0.1) +
                scale_color_continuous(name = gene, high = "red", low = "gray")
            }
          })


plotKneepoint <- function(cnts, sample.name = NULL){
  require(tidyverse)
  nUMI.c <- colSums(cnts)
  dc <- data.frame(cell = colnames(cnts),
                   nUMI = nUMI.c)

  dc <- dc %>% arrange(desc(nUMI))
  dc <- dc %>% mutate(id = 1:nrow(dc))
  dc <- dc %>% filter(id <= 2.5e4)
  dc <- dc %>% mutate(cum.nUMI = cumsum(nUMI))

  p1 <- ggplot(dc, aes(x = id, y = cum.nUMI)) +
    geom_point()

  if(!is.null(sample.name)){
    p1 <- p1 + ggtitle(sample.name)
  }

  print(
    p1
  )

  return(list(plot = p1, data = dc))
}

plotMultiKnee <- function(mcnts, sep = "_"){
  splitted.names.df <- do.call(
    rbind,
    strsplit(colnames(mcnts), split = sep)
  )

  samples <- levels(as.factor(splitted.names.df[,1]))

  out <- list()
  for(sample in samples){
    tmp <- mcnts[,startsWith(colnames(mcnts), sample)]
    out[[sample]] <- plotKneepoint(tmp, sample.name = sample)
  }

  return(out)
}

diffMultiSeurat <- function(ms.obj, lvl1, lvl2){
  s.lvl1 <- SubsetData(ms.obj, subset.name = "orig.ident", accept.value = lvl1)
  s.lvl2 <- SubsetData(ms.obj, subset.name = "orig.ident", accept.value = lvl2)

  gene.lvl1 <- rowSums(s.lvl1@data)
  gene.lvl2 <- rowSums(s.lvl2@data)

  plot(gene.lvl1, gene.lvl2, log = "xy", xlab = lvl1, ylab = lvl2)

  return(log2(gene.lvl1 / gene.lvl2))
}
