MyDoubletFinder <- function (seu, expected.doublets = 0, proportion.artificial = 0.25, 
                             proportion.NN = 0.01) 
{
  if (expected.doublets == 0) {
    stop("Need to set number of expected doublets...")
  }
  print("Creating artificial doublets...")
  data <- seu@assays$RNA@counts
  real.cells <- colnames(seu)
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1 - proportion.artificial) - 
                        n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("DOUBLET", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  print("Normalizing...")
  norm.data <- Seurat::LogNormalize(data = data_wdoublets, scale.factor = seu@commands$NormalizeData.RNA@params$scale.factor,
                                    verbose = F)
  print("Finding variable genes...")
  clip.max <- sqrt(x = ncol(x = norm.data))
  hvf.info <- data.frame(mean = rowMeans(x = as.matrix(norm.data)))
  hvf.info$variance <- .Call('_Seurat_SparseRowVar2', PACKAGE = 'Seurat', norm.data,hvf.info$mean,F)
  hvf.info$variance.expected <- 0
  not.const <- hvf.info$variance > 0
  fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
               data = hvf.info[not.const, ], span = 0.3)
  hvf.info$variance.expected[not.const] <- 10^fit$fitted
  hvf.info$variance.standardized <- .Call('_Seurat_SparseRowVarStd', PACKAGE = 'Seurat',
                                          norm.data,
                                          hvf.info$mean,
                                          sqrt(hvf.info$variance.expected),
                                          clip.max,F)
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 
                               0), ]
  hvf.info <- hvf.info[order(hvf.info$variance.standardized, 
                             decreasing = TRUE), , drop = FALSE]
  var_genes <- head(x = rownames(x = hvf.info), n = seu@commands$FindVariableFeatures.RNA@params$nfeatures)
  print("Running PCA...")
  rm(clip.max,hvf.info,fit)
  norm.data=.Call('_Seurat_FastSparseRowScale', PACKAGE = 'Seurat', norm.data[var_genes,],T,T,10,F)
  PCA.doublet=irlba::prcomp_irlba(t(norm.data),n = 10)
  cell.names <- colnames(data_wdoublets)
  nCells <- length(cell.names)
  print("Calculating PC distance matrix...")
  pca.coord <- PCA.doublet$x
  rm(PCA.doublet,norm.data,data_wdoublets)
  gc()
  PCdist=proxy::dist(pca.coord[1:n_real.cells,],pca.coord)
  pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
  rownames(pANN) <- real.cells
  colnames(pANN) <- "pANN"
  k <- round(nCells * proportion.NN)
  for (i in 1:n_real.cells) {
    neighbors <- order(PCdist[i,])[2:(k+1)]
    neighbor.names <- cell.names[neighbors]
    pANN[i, 1] <- length(grep("^DOUBLET", neighbor.names))/k
  }
  rm(PCdist)
  gc()
  seu <- Seurat::AddMetaData(seu, metadata = pANN, col.name = "pANN")
  predictions <- as.data.frame(rep("Singlet", n_real.cells),
                               ncol = 1, stringsAsFactors = FALSE)
  rownames(predictions) <- real.cells
  doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pANN,
                                                       decreasing = TRUE)]
  doublet.predictions <- doublet.predictions[1:expected.doublets]
  predictions[doublet.predictions, ] <- "Doublet"
  colnames(predictions) <- "pANNPredictions"
  seu <- Seurat::AddMetaData(seu, metadata = predictions, col.name = "pANNPredictions")
  return(seu)
}

MySeuratDR2Gg <- function(seurat.data,
                       x.dim = 1,
                       y.dim = 2,
                       reduction.use="tsne",
                       reduction.key="tSNE",
                       estimate.variation.explain.percentage=F,
                       color_manual = c(brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"))){
  library(ggplot2)
  coord <- seurat.data@reductions[[reduction.use]]@cell.embeddings
  coord <- as.data.frame(coord)
  
  sample.inf <- seurat.data@meta.data
  df <- cbind(x.pos = coord[,x.dim],
              y.pos = coord[,y.dim],
              sample.inf[row.names(coord),])
  
  p <- ggplot(data = df,
              mapping = aes(x = x.pos,
                            y = y.pos,
                            label = rownames(sample.inf)))
  
  if(estimate.variation.explain.percentage){
    coord.x.VEP=(seurat.data@reductions[[reduction.use]]@stdev[x.dim])^2/
      seurat.data@reductions[[reduction.use]]@misc$total.variance*100
    coord.y.VEP=(seurat.data@reductions[[reduction.use]]@stdev[y.dim])^2/
      seurat.data@reductions[[reduction.use]]@misc$total.variance*100
    coord.x.VEP=paste("(",round(coord.x.VEP,1),"%)",sep = "")
    coord.y.VEP=paste("(",round(coord.y.VEP,1),"%)",sep = "")
  }else{
    coord.x.VEP=NULL
    coord.y.VEP=NULL
  }
  p <- p + xlab(paste(reduction.key,x.dim,coord.x.VEP,sep = ""))
  p <- p + ylab(paste(reduction.key,y.dim,coord.y.VEP,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=1)
  return(p)
}
