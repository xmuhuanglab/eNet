# Tutorial for eNet
Five input files are needed for eNet:
 1. scATAC-seq matrix;
 2. scRNA-seq matrix;
 3. UMAP coordinates; # A data frame with columns representing the coordinates of each cell in reduced dimension space
 4. Cell metadata; # A data frame containing attributes of individual cells
 5. Reference genome. # Obtained from https://github.com/buenrostrolab/stimATAC_analyses_code
### Step 1. Preparing input data (Input)
```r
setwd("workdir")
rm(list = ls())
library(Signac)
library(IRanges)
library(Seurat)
library(reshape2)
library(TFBSTools)
library(dplyr)
library(SummarizedExperiment)
library(Matrix)

# ---------- Load scRNA matrix and scATAC matrix
load("../data/cre.mat.Rdata") # peak-cell matrix
load("../data/rna.mat.Rdata") # scRNA matrix
# Note that the scATAC-seq and scRNA-seq matrix must have the same columns.
load("../data/dcluster_coords.Rdata") # UMAP coordinates
load("../data/metadata.Rdata") # Cell metadata
source('./MainFunc.R')
source('./utils.R')
```
### Step 2. Identifying the putative enhancer cluster (Node)
```r
# --- 2.1 calculate gene-peak correlation
GPTab <- GPCor(cre.mat=cre.mat,  # peak-cell matrix
               exp.mat=rna.mat,  # scRNA-seq matrix or gene activity matrix
               normalizeRNAMat=T, # if the rna matrix need normalization
               genome = "hg19", # reference genome, must be one of "hg19", "mm10", or "hg38"
               windowPadSize = 100000, # base pairs padded on either side of gene TSS
               proPadSize = 2000, # base pairs padded on either side of gene TSS for enhancer
               nCores=8 # How many registerCores to use
)
save(GPTab, file = 'genePeakTab.Rdata')
# --- 2.2 obtain significantly correlated peak-gene pairs
# For determining the threshold value of estimate, namely peak-gene correlation, we advice to choose the value of quantile 95% of overall estimate using quantile(GPTab$estimate, seq(0,1,0.05))[['95%']]
GPTabFilt <- FindNode(GPTab = GPTab, # data frame of gene-peak correlation
                      genome = "hg19", # reference genome, must be one of "hg19", "mm10", or "hg38"
                      estimate = 0, # the threshold value of peak-gene correlation to determine whether the peak-gene pair is significantly correlated
                      proPadSize = 2000, # base pairs padded on either side of gene TSS for enhancer
                      FDR = 0.05 # the threshold value of FDR to determine whether the peak-gene pair is significantly correlated
)
save(GPTabFilt, file = 'genePeakTabFilt.Rdata')
# --- 2.3 convert the table of gene-peak pairs to a list of enhancer cluster
GPPair <- list()
for(i in unique(GPTabFilt$Gene)){
  tmp <- subset(GPTabFilt, Gene == i)
  GPPair[[i]] <- as.character(tmp$Peak)
}
save(GPPair, file = 'GPPair.Rdata')
```
### Step 3. Identifying the predicted enhancer interactions (Edge)
```r
library(monocle3)
library(Signac)
library(Seurat)
library(cicero)
conns <- FindEdge(peaks.mat=cre.mat,  # peak-cell matrix
                  GPPair=GPPair, # A list of enhancer cluster, the output of Step.2
                  cellinfo=metadata,  # A data frame containing attributes of individual cells.
                  k=50, # Number of cells to aggregate per bin when generate an aggregated input CDS for cicero
                  coords=dcluster_coords, # A data frame with columns representing the coordinates of each cell in reduced dimension space (generally 2-3 dimensions).
                  genome='hg19' # reference genome, must be one of "hg19", "mm10", or "hg38"
)
save(conns, file = 'conns_GPPairPeaks.Rdata')
```
### Step 4. Building enhancer networks (Network)
```r
library(igraph)
library(dplyr)
library(doParallel)
# To determine the threshold value of co-accssibility score, we advice to choose the value around quantile 90%-95% using quantile(subset(conns, coaccess>0)$coaccess,seq(0,1,0.05),na.rm = T)
NetworkList <- BuildNetwork(conns=conns, # A data frame of co-accessibility scores, also the output of Step.3
                            GPTab=GPTabFilt,  # a list of enhancer cluster, the output of Step.2
                            cutoff=0.1, # the cutoff of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
                            nCores=8 # How many cores to use
)
save(NetworkList, file = "NetworkList.Rdata")
# Visualize the enhancer network
plot.igraph(NetworkList[["ESPN"]])
```
### Step 5.Calculating network complexity(Network complexity)
```r
# To determine the threshold value of co-accssibility score, we advice to choose the value around quantile 90%-95% using quantile(subset(conns, coaccess>0)$coaccess,seq(0,1,0.05),na.rm = T)
Networkinfo <- NetComplexity(conns=conns,  # enhancer-enhancer co-accessibilty calculated using Cicero
                             GPTab=GPTabFilt,   # a list of enhancer cluster, also the output of Step.2
                             cutoff=0.1,  # the threshold value of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
                             nCores=8 # How many cores to use
)
save(Networkinfo, file = "Networkinfo.Rdata")
```
### Step 6. Classification of enhancer networks (Mode)
```r
# To determine the threshold value of network size, namely SizeCutoff in the function NetworkMode, we advice to choose the median of network size for all genes.
# To determine the threshold value of network connectivity, namely ConnectivityCutoff in the function NetworkMode, we advice to use value 1 to ensure that on average each enhancer has at least one interaction with other enhancer.
Mode <- NetworkMode(Networkinfo=Networkinfo,  # Network information, the output file in Step.5
                    SizeCutoff=5, # this parameter represents the threshold value of network size, which can be used to distinguish Simple and Network/Multiple mode
                    ConnectivityCutoff=1 # this parameter represents the threshold value of network connectivity, which can be used to distinguish Network and Multiple mode
)
save(Mode, file = 'NetworkMode.Rdata')

# Visualize the three modes of enhancer networks 
library(ggplot2)
ggplot(Mode, aes(x=NetworkSize, y=NetworkConnectivity, color=Mode), alpha=0.8) + geom_point() +
  theme_classic() + xlab('log2(NetworkSize)') + ylab('Network connectivity')
```
