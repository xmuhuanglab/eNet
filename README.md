# eNet: an algorithm to build enhancer networks based on scATAC-seq and scRNA-seq data
## Introduction
eNet is an algorithm designed to integrate single-cell chromatin accessibility and gene expression profiles and build enhancer networks, delineating how multiple enhancers interact with each other in gene regulation. 

## Workflow
### Step1. Preparing input matrix (Input)
Two matrices are needed for the input of eNet. 1) scATAC-seq matrix (peak-cell); 2) scRNA-seq matrix (gene-cell).
```r
load('./cre.mat.Rdata')
load('./rna.mat.Rdata')
```
### Step2. Identifying the putative enhancer cluster (Node)
```r
GPTab <- GPCor(cre.mat=cre.mat,  
               exp.mat=rna.mat,  
               normalizeRNAMat=T, 
               genome = "hg19", 
               windowPadSize = 100000, 
               proPadSize = 2000,
               nCores=8
)
GPTabFilt <- FindNode(GPTab = GPTab, 
                      genome = "hg19", 
                      estimate = 0, 
                      proPadSize = 2000, 
                      FDR = 0.05 # 
)
GPPair <- list()
for(i in unique(GPTabFilt$Gene)){
  tmp <- subset(GPTabFilt, Gene == i)
  GPPair[[i]] <- as.character(tmp$Peak)
}
```
### Step3. Identifying the predicted enhancer interactions (Edge)
```r
conns <- FindEdge(peaks.mat=cre.mat, 
                  GPPair=GPPair,
                  cellinfo=metadata, 
                  k=50, 
                  coords=dcluster_coords, 
                  genome='hg19' 
)
```
### Step4. Building enhancer networks (Network)
```r
NetworkList <- BuildNetwork(conns=conns, 
                            GPTab=GPTabFilt,  
                            cutoff=0.1, 
                            nCores=8 
)
```
### Step5. Calculating network complexity (Network complexity)
```r
Networkinfo <- NetComplexity(conns=conns,  
                             GPTab=GPTabFilt,   
                             cutoff=0.1,
                             nCores=8 
)
```
### Step6. Classification of enhancer networks (Mode)
```r
Mode <- NetworkMode(Networkinfo=Networkinfo,  
                    SizeCutoff=5, 
                    ConnectivityCutoff=1
)
```

## How to cite eNet
1. Danni Hong#, Hongli Lin#, Lifang Liu, Muya Shu, Jianwu Dai, Falong Lu, Jialiang Huang*. Complexity of enhancer networks predicts cell identity and disease genes. (Submitted)
2. Muya Shu#, Danni Hong#, Hongli Lin, Jixiang Zhang, Zhengnan Luo, Yi Du, Zheng Sun, Man Yin, Yanyun Yin, Shilai Bao, Zhiyong Liu, Falong Lu*, Jialiang Huang*, Jianwu Dai*. Enhancer networks driving mouse spinal cord development revealed by single-cell multi-omics analysis. (Submitted)
