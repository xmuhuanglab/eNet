# eNet: an algorithm to build enhancer networks based on scATAC-seq and scRNA-seq data
## Introduction
eNet is an algorithm designed to integrate single-cell chromatin accessibility and gene expression profiles and build enhancer networks, delineating how multiple enhancers interact with each other in gene regulation. 

## Workflow
### Step1. Preparing input matrix (Input)
Two matrices are needed for the input of eNet. 1) scATAC-seq matrix (peak-cell); 2) scRNA-seq matrix (gene-cell).
### Step2. Identifying the putative enhancer cluster (Node)
GPPair <- FindNode(cre.mat=cre.mat,  # peak-cell matrix
                  exp.mat=rna.mat,  # scRNA-seq matrix
                  normalizeRNAMat=T, # if the rna matrix need normalization
                  genome = "hg19", # reference genome, must be one of "hg19", "mm10", or "hg38"
                  windowPadSize = 100000, # base pairs padded on either side of gene TSS
                  proPadSize = 2000, # base pairs padded on either side of gene TSS for enhancer
                  estimate = 0, # the threshold value of peak-gene correlation to determine whether the peak-gene pair is significantly correlated
                  FDR = 0.05 # the threshold value of FDR to determine whether the peak-gene pair is significantly correlated
)
### Step3. Identifying the predicted enhancer interactions (Edge)
conns <- FindEdge(peaks.mat=cre.mat,  # peak-cell matrix
                  GPPair=GPPair, # A list of enhancer cluster, the output of Step.2
                  cellinfo=metadata,  # A data frame containing attributes of individual cells.
                  k=50, # Number of cells to aggregate per bin when generate an aggregated input CDS for cicero
                  coords=dcluster_coords # A data frame with columns representing the coordinates of each cell in reduced dimension space (generally 2-3 dimensions).
)
### Step4. Building enhancer networks (Network)
NetworkList <- NetworkBuilding(conns=conns, # A data frame of co-accessibility scores, also the output of Step.3
                               GPTab=GPTabFilt,  # A data frame of filtered significantly correlated peak-gene pairs generated in Step.2
                               cutoff=0.1 # the cutoff of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
)
You can visualize the network by running the command: plot.igraph(NetworkList[[gene]])
### Step5. Calculating network complexity (Network complexity)
Networkinfo <- NetComplexity(conns=conns,  # enhancer-enhancer co-accessibilty calculated using Cicero
                             GPTab=GPTabFilt,  # A data frame of filtered significantly correlated peak-gene pairs generated in Step.2
                             cutoff=0.1  # the threshold value of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
)
### Step6. Classification of enhancer networks (Mode)
Mode <- NetworkMode(Networkinfo=Networkinfo,  # Network information, the output file in Step.5
                    SizeCutoff=5, # this parameter represents the threshold value of network size, which can be used to distinguish Simple and Network/Multiple mode
                    ConnectivityCutoff=1 # this parameter represents the threshold value of network connectivity, which can be used to distinguish Network and Multiple mode
)

## How to cite eNet
1. Danni Hong#, Hongli Lin#, Lifang Liu, Muya Shu, Jianwu Dai, Falong Lu, Jialiang Huang*. Complexity of enhancer networks predicts cell identity and disease genes. (Submitted)
2. Muya Shu#, Danni Hong#, Hongli Lin, Jixiang Zhang, Zhengnan Luo, Yi Du, Zheng Sun, Man Yin, Yanyun Yin, Shilai Bao, Zhiyong Liu, Falong Lu*, Jialiang Huang*, Jianwu Dai*. Enhancer networks driving mouse spinal cord development revealed by single-cell multi-omics analysis. (Submitted)
