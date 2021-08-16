# SCRIP
SCRIP: an accurate simulator for single-cell RNA sequencing data

## Author
Fei Qin, Xizhi Luo, Feifei Xiao, Guoshuai Cai

## Description
SCRIP provides a flexible Gamma-Poisson mixture and a Beta-Poisson mixture framework to simulate scRNA-seq data. SCRIP package was built based on the framework of splatter. Both Gamma-Poisson and Beta-Poisson distribution model the over dispersion of scRNA-seq data. Specifically, Beta-Poisson model was used to model bursting effect. The dispersion was accurately simulated by fitting the mean-BCV dependency using generalized additive model (GAM). Other key characteristics of scRNA-seq data including library size, zero inflation and outliers were also modeled by SCRIP. With its flexible modeling, SCIRP enables various application for different experimental designs and goals including DE analysis, clustering analysis, trajectory-based analysis and bursting analysis

## Installation
```r
install.packages("devtools")
library(devtools)
BiocManager::install("splatter")
install_github("thecailab/SCRIP")
```

## Running SCRIP
### Examples

```r
library(SCRIP)
library(Biobase)
library(splatter)

EMTAB.eset <- EMTABesethealthy
expre_data = exprs(EMTAB.eset)
pheno_data = pData(EMTAB.eset)

dim(expre_data)
[1] 25453  1097

expre_data[1:6,1:6]
        AZ_A10 AZ_A11 AZ_A12 AZ_A2 AZ_A5 AZ_A6
SGIP1        0      0      0    32     0     0
AZIN2        0      0      0     0     0     0
CLIC4        3      0      0     1     0     0
AGBL4        0      0      0     0     0     0
NECAP2       0      0      0     0     0     0
SLC45A1      0      0      0     0     0     0

head(pheno_data)
 sampleID SubjectName cellTypeID cellType
AZ_A10        1   Non T2D 1          5    delta
AZ_A11        1   Non T2D 1          2    alpha
AZ_A12        1   Non T2D 1          5    delta
AZ_A2         1   Non T2D 1          9    gamma
AZ_A5         1   Non T2D 1          6   ductal
AZ_A6         1   Non T2D 1          2    alpha
 

# To simplify computation here, only acinar celltype in sample 5 was utilized here. 
sub_expre_data=expre_data[,which((pheno_data$sampleID %in% c(5)) & 
                                   (pheno_data[,4] == "acinar"))]
params <- splatEstimate(sub_expre_data)
sim_trend <-  SCRIPsimu(data=sub_expre_data, params=params, mode="GP-trendedBCV")

sim_trend
class: SingleCellExperiment 
dim: 25453 80 
metadata(16): Params CT.index ... batch.facScale bcv.shrink
assays(5): BatchCellMeans BaseCellMeans CellMeans TrueCounts counts
rownames(25453): Gene1 Gene2 ... Gene25452 Gene25453
rowData names(4): Gene BaseGeneMean OutlierFactor GeneMean
colnames(80): Cell1 Cell2 ... Cell79 Cell80
colData names(3): Cell Batch ExpLibSize
reducedDimNames(0):
altExpNames(0):

```
