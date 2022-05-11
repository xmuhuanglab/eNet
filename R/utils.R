# Nomalized data
centerCounts <- function(obj,
                         doInChunks=TRUE,
                         chunkSize=1000){
  if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  
  if(ncol(obj) > 10000)
    doInChunks <- TRUE
  
  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }
  
  counts.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }
    
    cat("Computing centered counts for cells: ",beginning," to ", ending,"..\n")
    
    if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    counts.l[[i]] <- cCounts
    
    gc()
  }
  
  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")
  
  if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}

#Gene peak overlaps
genePeakOv <- function(ATAC.se, # SummarizedExperiment object of scATAC data
                       RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                       genome, # Must be one of "hg19", "mm10", or "hg38"
                       geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                       windowPadSize=100000, # base pairs padded on either side of gene TSS
                       proPadSize = 2000 # base pairs padded on either side of gene TSS for enhancer
                        
) {
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  ATACmat <- assay(ATAC.se) # Rownames preserved
  
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATAC.se) # Peak ranges
  names(peakRanges) <- rownames(ATAC.se)
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  
  
  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  # switch(genome, hg19 = {
  #   TSSg <- hg19TSSRanges
  # }, hg38 = {
  #   TSSg <- hg38TSSRanges
  # }, mm10 = {
  #   TSSg <- mm10TSSRanges
  # })
  if(genome %in% c("hg19", "hg38", "mm10")){
    load(paste0('../data/', genome, '_refSeq.Rdata'))
    TSSg <- get(paste0(genome, 'TSSRanges'))
  }
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Pad promoter by this much *either side*
  proflank <- GenomicRanges::flank(TSSg,
                                   width = proPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping promoter-gene pairs ..\n")
  geneProOv <- findOverlaps(query = proflank,subject = peakSummits)
  numpgPairs <- length(geneProOv)
  
  cat("Found ",numpgPairs,"total promoter-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene promoter window: ",length(unique(subjectHits(geneProOv))),"\n")
  cat("Number of gene promoter windows that overlap any peak summit: ",length(unique(queryHits(geneProOv))),"\n\n")
  peakSummits <- peakSummits[-unique(subjectHits(geneProOv))]
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlapPairs(query = TSSflank,subject = peakSummits)
  genePeakOv <- data.frame(genePeakOv@first,genePeakOv@second)
  numPairs <- nrow(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(genePeakOv$Peak)),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(genePeakOv$gene_name)),"\n\n")
  
  return(genePeakOv)

}

# Gene peak correlation
PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        peakRanges,
                        seed = 2022,
                        mtd="spearman") {
  
  
  geneIndices <- as.character(OV$gene_name)
  peakIndices <- as.character(OV$Peak)
  
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  
  
  idy = uniquepeaks
  ovlpn = length(idy)
  idyexcl = which(!1:length(peakRanges) %in% idy)
  gene.name = uniquegenes
  
  gene.val = RNA[uniquegenes,];
  
  data.raw = ATAC
  peak.raw = peakRanges
  # use set
  data.use = ATAC[idy,];
  peak.use = peak.raw[uniquepeaks];
  
  # shuf set
  data.shuf <- data.use
  peak.shuf <- peak.use
  colheader <- colnames(data.shuf)
  #Shuffle row-wise:
  set.seed(seed)
  data.shuf <- data.shuf[sample(nrow(data.shuf)),]
  #Shuffle column-wise:
  set.seed(seed)
  data.shuf <- data.shuf[,sample(ncol(data.shuf)), with=F]
  colnames(data.shuf) <- colheader
  
  # rdm set
  set.seed(seed)
  
  #    rdm <- sample(idyexcl, 1000)
  rdm <- sample(idyexcl, ovlpn)
  data.rdm = data.raw[rdm,];
  peak.rdm = peak.raw[rdm];
  
  
  # rdm shuf set
  data.rdmshuf <- data.rdm
  peak.rdmshuf <- peak.rdm
  colheader <- colnames(data.rdmshuf)
  #Shuffle row-wise:
  set.seed(seed)
  data.rdmshuf <- data.rdmshuf[sample(nrow(data.rdmshuf)),]
  #Shuffle column-wise:
  set.seed(seed)
  data.rdmshuf <- data.rdmshuf[,sample(ncol(data.rdmshuf)), with=F]
  colnames(data.rdmshuf) <- colheader
  
  corrFunc <- function(var1, var2, method) {
    result = cor.test(var1, var2, method = method)
    data.frame(result[c("estimate","p.value","statistic")],
               stringsAsFactors=FALSE)
  }
  
  #-----------------------------------
  #cat("Real set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.use));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.use[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.use$estimate <- corr[, "estimate"]
  peak.use$statistic <- corr[, "statistic"]
  peak.use$method <- mtd
  peak.use$Pval <- corr[, "p.value"]
  peak.use$FDR <- p.adjust(peak.use$Pval, method = "BH")
  peak.use$class <- "corr"
  peak.use$Gene <- gene.name
  
  #cat("Done ... \n", file = stderr())
  
  #------------------------------
  #cat("Shuffle set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.shuf));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.shuf[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.shuf$estimate <- corr[, "estimate"]
  peak.shuf$statistic <- corr[, "statistic"]
  peak.shuf$method <- mtd
  peak.shuf$Pval <- corr[, "p.value"]
  peak.shuf$FDR <- p.adjust(peak.shuf$Pval, method = "BH")
  peak.shuf$class <- "shuf"
  peak.shuf$Gene <- gene.name
  
  #cat("Done ... \n", file = stderr())
  
  #--------------------------------
  #cat("Random set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.rdm));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.rdm[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  
  peak.rdm$estimate <- corr[, "estimate"]
  peak.rdm$statistic <- corr[, "statistic"]
  peak.rdm$method <- mtd
  peak.rdm$Pval <- corr[, "p.value"]
  peak.rdm$FDR <- p.adjust(peak.rdm$Pval, method = "BH")
  peak.rdm$class <- "random"
  peak.rdm$Gene <- gene.name
  #cat("Done ... \n", file = stderr())
  
  
  #----------------------
  #cat("Rdmshuf set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.rdmshuf));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.rdmshuf[t,]), mtd)
  })
  corr <- do.call(rbind, corr)
  peak.rdmshuf$estimate <- corr[, "estimate"]
  peak.rdmshuf$statistic <- corr[, "statistic"]
  peak.rdmshuf$method <- mtd
  peak.rdmshuf$Pval <- corr[, "p.value"]
  peak.rdmshuf$FDR <- p.adjust(peak.rdmshuf$Pval, method = "BH")
  peak.rdmshuf$class <- "rdmShuf"
  peak.rdmshuf$Gene <- gene.name
  #cat("Done ... \n", file = stderr())
  
  #############################
  # OUTPUT
  peak.use.df <- as.data.frame(peak.use)
  peak.shuf.df <- as.data.frame(peak.shuf)
  peak.rdm.df <- as.data.frame(peak.rdm)
  peak.rdmshuf.df <- as.data.frame(peak.rdmshuf)
  peak.df <- rbind(peak.use.df, peak.shuf.df)
  peak.df <- rbind(peak.df, peak.rdm.df)
  peak.df <- rbind(peak.df, peak.rdmshuf.df)
  
  return(peak.df);
  
}


# build enhancer network for each gene
Network <- function(conns=conns, # A data frame of co-accessibility scores, also the output of Step.3
                    GPTab=GPTabFilt,  # A list of enhancer cluster, the output of Step.2
                    cutoff=0.1 # The cutoff of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
){
  se <- as.character(GPTab$Peak)
  se <- se[order(se)]
  conn <- conns %>% filter(Peak1 %in% se & Peak2 %in% se)
  data <- matrix(nrow = length(se),ncol = length(se))
  rownames(data) <- se
  colnames(data) <- se
  conn$coaccess[which(is.na(conn$coaccess))] <- 0
  conn$coaccess[which(conn$coaccess < cutoff)] <- 0 
  conn$coaccess[which(conn$coaccess >= cutoff)] <- 1
  for (m in 1:nrow(data)) {
    se1 <- se[m]
    co <- conn[which(conn$Peak1 == se1),]
    rownames(co) <- co$Peak2
    co <- co[se[-m],]
    data[m,-m] <- co$coaccess
    data[m,m] <- 0
  }
  data[is.na(data)] <- 0
  rownames(data) <- colnames(data) <- 1:length(rownames(data))
  if(max(data)>0){
    net <- graph_from_adjacency_matrix(as.matrix(data), mode="undirected")
    return(net)
    #names(net) <- gene
  }
  
}

# calculate several metrics for each enhancer network
NetworkComplexity <- function(net=net,
                              gene = gene,
                              GPTab=GPTabFilt)
  {
  if(class(net) == "NULL"){
    data <- data.frame(gene=gene,
                       EdgeNum=0, 
                       NetworkSize=nrow(GPTab),
                       MaxDegree=0, 
                       NetworkConnectivity=0, 
                       row.names = gene)
  }else{
    EdgeNum <- length(E(net))*2
    NetworkSize <- length(V(net))
    MaxDegree <- max(degree(net, v = V(net), mode = 'all'))
    NetworkConnectivity <- EdgeNum/NetworkSize
    data <- data.frame(gene=gene,
                       EdgeNum=EdgeNum, 
                       NetworkSize=NetworkSize,
                       MaxDegree=MaxDegree, 
                       NetworkConnectivity=NetworkConnectivity, 
                       row.names = gene)

  }
  
  return(data)
}
