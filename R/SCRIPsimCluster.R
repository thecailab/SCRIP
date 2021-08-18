#' #SCRIP group simulation function for clustering analysis
#'
#' @export
simu.VEGs <- function(counts.matrix, params=params, base_allcellmeans, mode="GP-trendedBCV", nCells, nfeatures=1000){
  message("Starting simulating SCRIP")
  rownames(counts.matrix) <- paste0("Gene",1:nrow(counts.matrix))
  colnames(counts.matrix) <- paste0("Cell",1:ncol(counts.matrix))

  pbmc <- Seurat::CreateSeuratObject(counts = counts.matrix, project = "pbmc3k", min.cells = 3, min.features = 200)
  seurat_stim <- Seurat::NormalizeData(pbmc,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)
  stim <- Seurat::FindVariableFeatures(object = seurat_stim,
                               selection.method = "vst",
                               nfeatures = nfeatures)
  top <- head(VariableFeatures(stim),nfeatures)
  EGcounts <- counts.matrix[top,]

  com_base_cellmeans <- base_allcellmeans
  names(com_base_cellmeans) <- rownames(counts.matrix)

  com_base_cellmeans[top] <- apply(EGcounts,1,mean)

  #### SCRIP ####
  if (mode=="BP-trendedBCV"){
    message("Starting simulating SCRIP BP-trendedBCV model")
    SCRIP.Trend.burst <- SCRIPsimu(data=counts.matrix, params=params, batchCells=nCells,
                                   base_allcellmeans_SC=com_base_cellmeans,
                                   mode="BP-trendedBCV")
    res <- counts(SCRIP.Trend.burst)
  }

  if (mode=="BP-commonBCV"){
    message("Starting simulating SCRIP BP-commonBCV model")
    SCRIP.common.burst <-  SCRIPsimu(data=counts.matrix, params=params, batchCells=nCells,
                                   base_allcellmeans_SC=com_base_cellmeans,
                                   mode="BP-commonBCV")
    res <- counts(SCRIP.common.burst)
  }


  if (mode=="GP-trendedBCV"){
    message("Starting simulating SCRIP GP-trendedBCV model")
    SCRIP.Trend.noburst <- SCRIPsimu(data=counts.matrix, params=params, batchCells=nCells,
                                   base_allcellmeans_SC=com_base_cellmeans,
                                   mode="GP-trendedBCV")
    res <- counts(SCRIP.Trend.noburst)
  }


  if (mode=="GP-commonBCV"){
    message("Starting simulating SCRIP GP-commonBCV model")
    SCRIP.common.noburst <-  SCRIPsimu(data=counts.matrix, params=params, batchCells=nCells,
                                     base_allcellmeans_SC=com_base_cellmeans,
                                     mode="GP-commonBCV")
    res <- counts(SCRIP.common.noburst)
  }

  return(res)
}



simu_cluster <- function(expre_data, pheno_data, CT, mode, nfeatures, seed=2021){
  set.seed(seed)
  params <- splatter::splatEstimate(expre_data)
  parnGene <- params@nGenes
  parshape <- params@mean.shape
  parrate <- params@mean.rate
  base_allcellmeans=rgamma(parnGene, shape=parshape, rate=parrate)

  final <- CT.infor <- NULL
  for (CT in CTlist){
    counts <- expre_data[,which(pheno_data$cellType==CT)]
    message(paste0("starting simulating cell type: ", CT))
    res <- simu.VEGs(counts.matrix = counts,params=params, base_allcellmeans=base_allcellmeans, nCells=ncol(counts), mode=mode, nfeatures = nfeatures)
    final <- cbind(final, res)
    CT.infor <- c(CT.infor, rep(CT, ncol(counts)))
  }
  final.list <- list(final=final, CT.infor=CT.infor)
  return(final.list)
}




