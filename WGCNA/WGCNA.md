WGCNA - Network Construction
================
Margarida Ferreira
Last updated July 2021

-   [1. Initial setup](#initial-setup)
-   [2. Data input and pre-processing](#data-input-and-pre-processing)
    -   [2.1 Gene expression data input](#gene-expression-data-input)
    -   [2.2 Quality check of genes and samples](#quality-check-of-genes-and-samples)
    -   [2.3 Trait data input](#trait-data-input)
    -   [2.4 Save expression and trait data](#save-expression-and-trait-data)
-   [3. Gene network constructuin and identification of modules](#gene-network-constructuin-and-identification-of-modules)
    -   [3.1. Choice of soft-thresholding power: analysis of network topology](#choice-of-soft-thresholding-power-analysis-of-network-topology)
        -   [Brain](#brain)
            -   [Fig S6-Brain](#fig-s6-brain)
        -   [Heart](#heart)
            -   [Fig S6-Heart](#fig-s6-heart)
        -   [Muscle](#muscle)
            -   [Fig S6-Muscle](#fig-s6-muscle)
        -   [Liver](#liver)
            -   [Fig S6-Liver](#fig-s6-liver)
        -   [Pancreas](#pancreas)
            -   [Fig S6-Pancreas](#fig-s6-pancreas)
    -   [3.2. Block-wise network construction and module detection](#block-wise-network-construction-and-module-detection)
        -   [Brain](#brain-1)
        -   [Heart](#heart-1)
        -   [Muscle](#muscle-1)
        -   [Liver](#liver-1)
        -   [Pancreas](#pancreas-1)

The WGCNA analysis was carried out using the package WGCNA (Langfelder and Horvath 2008) and mainly followed the [tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/) provided by the authors.

# 1. Initial setup

Loading the required packages:

``` r
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
library('dlfUtils')
```

# 2. Data input and pre-processing

## 2.1 Gene expression data input

We use DESeq2's variance stabilizing transformation (performed in the `Normalization.md`).

``` r
#brain
brainData = read.table("vst_brain_counts.txt", header=TRUE)
dim(brainData)
```

    ## [1] 34164    46

``` r
brainDataExpr = as.data.frame(t(brainData)) #note the transposition of the data frame
dim(brainDataExpr)
```

    ## [1]    46 34164

``` r
#heart
heartData = read.table("vst_heart_counts.txt", header=TRUE)
dim(heartData)
```

    ## [1] 28073    45

``` r
heartDataExpr = as.data.frame(t(heartData)) #note the transposition of the data frame
dim(heartDataExpr)
```

    ## [1]    45 28073

``` r
#muscle
muscleData = read.table("vst_muscle_counts.txt", header=TRUE)
dim(muscleData)
```

    ## [1] 18978    43

``` r
muscleDataExpr = as.data.frame(t(muscleData)) #note the transposition of the data frame
dim(muscleDataExpr)
```

    ## [1]    43 18978

``` r
#liver
liverData = read.table("vst_liver_counts.txt", header=TRUE)
dim(liverData)
```

    ## [1] 20157    45

``` r
liverDataExpr = as.data.frame(t(liverData)) #note the transposition of the data frame
dim(liverDataExpr)
```

    ## [1]    45 20157

``` r
#pancreas
pancreasData = read.table("vst_pancreas_counts.txt", header=TRUE)
dim(pancreasData)
```

    ## [1] 18414    32

``` r
pancreasDataExpr = as.data.frame(t(pancreasData)) #note the transposition of the data frame
dim(pancreasDataExpr)
```

    ## [1]    32 18414

## 2.2 Quality check of genes and samples

As proposed by the authors, We used the function `goodSamplesGenes` to check the data for *"missing entries, entries with weigths below a threshold, and zero-variance genes*. Note that we used this function with default values which can be consulted [here](https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/goodSamplesGenes).

``` r
brain_gsg = goodSamplesGenes(brainDataExpr)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
print("Brain")
```

    ## [1] "Brain"

``` r
brain_gsg$allOK
```

    ## [1] TRUE

``` r
heart_gsg = goodSamplesGenes(heartDataExpr)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
print("Heart")
```

    ## [1] "Heart"

``` r
heart_gsg$allOK
```

    ## [1] TRUE

``` r
muscle_gsg = goodSamplesGenes(muscleDataExpr)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1
    ##   ..step 2

``` r
print("Muscle")
```

    ## [1] "Muscle"

``` r
muscle_gsg$allOK
```

    ## [1] FALSE

``` r
liver_gsg = goodSamplesGenes(liverDataExpr)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
print("Liver")
```

    ## [1] "Liver"

``` r
liver_gsg$allOK
```

    ## [1] TRUE

``` r
pancreas_gsg = goodSamplesGenes(pancreasDataExpr)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1
    ##   ..step 2

``` r
print("Pancreas")
```

    ## [1] "Pancreas"

``` r
pancreas_gsg$allOK
```

    ## [1] FALSE

If `TRUE`, then all genes have passed the filters and, consequently, retained for further analyses. In this case, the muscle and the pancreas have genes that were flagged and need to be removed before moving on in the analysis:

``` r
##MUSCLE
if (!muscle_gsg$allOK)
{
print("Muscle")
#Print the gene and sample names that were removed
if (sum(!muscle_gsg$goodGenes)>0)
print(paste("Removing genes:", paste(names(muscleDataExpr)[!muscle_gsg$goodGenes], collapse = ", ")));
if (sum(!muscle_gsg$goodSamples)>0)
print(paste("Removing samples:", paste(rownames(muscleDataExpr)[!muscle_gsg$goodSamples], collapse = ", ")));
#Remove from the expression matrix the genes and samples flagged as "not good"
muscleDataExpr = muscleDataExpr[muscle_gsg$goodSamples, muscle_gsg$goodGenes]
}
```

    ## [1] "Muscle"
    ## [1] "Removing genes: Mup9"

``` r
dim(muscleDataExpr)
```

    ## [1]    43 18977

``` r
##PANCREAS
if (!pancreas_gsg$allOK)
{
print("Pancreas")
#Print the gene and sample names that were removed
if (sum(!pancreas_gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(pancreasDataExpr)[!pancreas_gsg$goodGenes], collapse = ", ")));
if (sum(!pancreas_gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(pancreasDataExpr)[!pancreas_gsg$goodSamples], collapse = ", ")));
#Remove from the expression matrix the genes and samples flagged as "not good"
pancreasDataExpr = pancreasDataExpr[pancreas_gsg$goodSamples, pancreas_gsg$goodGenes]
}
```

    ## [1] "Pancreas"
    ## Removing genes: Gm31659, Gm6018, Mir208b

``` r
dim(pancreasDataExpr)
```

    ## [1]    32 18411

Notably, very few genes were flagged as not passing the filtering thresholds, which validates the filtering step performed previously in the DEA workflow.

Next in the WGCNA tutorial is the outlier removal step. However, because we've already followed a network approach to identify and remove outlier samples in the DEA worflow, we did not performed this step to avoid the removal of too many samples and eventually remove interesting biological variation.

## 2.3 Trait data input

We now read in the trait data and match the samples for which they were measured to the expression samples.

``` r
##BRAIN
#read in the sample information table
brain_traitData = read.table("coldata_brain.txt", header=TRUE)
dim(brain_traitData)
```

    ## [1] 46  4

``` r
#binarize sex 
brain_traitData$Sex[brain_traitData$Sex == "Female"] <- "0"
brain_traitData$Sex[brain_traitData$Sex == "Male"] <- "1"
brain_traitData$Sex <- as.numeric(brain_traitData$Sex)

##HEART
#read in the sample information table
heart_traitData = read.table("coldata_heart.txt", header=TRUE)
dim(heart_traitData)
```

    ## [1] 45  4

``` r
#binarize sex 
heart_traitData$Sex[heart_traitData$Sex == "Female"] <- "0"
heart_traitData$Sex[heart_traitData$Sex == "Male"] <- "1"
heart_traitData$Sex <- as.numeric(heart_traitData$Sex)

##MUSCLE
#read in the sample information table
muscle_traitData = read.table("coldata_muscle.txt", header=TRUE)
dim(muscle_traitData)
```

    ## [1] 43  4

``` r
#binarize sex 
muscle_traitData$Sex[muscle_traitData$Sex == "Female"] <- "0"
muscle_traitData$Sex[muscle_traitData$Sex == "Male"] <- "1"
muscle_traitData$Sex <- as.numeric(muscle_traitData$Sex)

##LIVER
#read in the sample information table
liver_traitData = read.table("coldata_liver.txt", header=TRUE)
dim(liver_traitData)
```

    ## [1] 45  4

``` r
#binarize sex 
liver_traitData$Sex[liver_traitData$Sex == "Female"] <- "0"
liver_traitData$Sex[liver_traitData$Sex == "Male"] <- "1"
liver_traitData$Sex <- as.numeric(liver_traitData$Sex)

##PANCREAS
#read in the sample information table
pancreas_traitData = read.table("coldata_pancreas.txt", header=TRUE)
dim(pancreas_traitData)
```

    ## [1] 32  4

``` r
#binarize sex 
pancreas_traitData$Sex[pancreas_traitData$Sex == "Female"] <- "0"
pancreas_traitData$Sex[pancreas_traitData$Sex == "Male"] <- "1"
pancreas_traitData$Sex <- as.numeric(pancreas_traitData$Sex)
```

To visualize how the phenotype traits relate to the sample dendogram:

``` r
##BRAIN
#Re-cluster samples
brain_sampleTree = hclust(dist(brainDataExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
brain_traitColors = numbers2colors(brain_traitData[3:4], signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
tiff("sampleDendro_traitHeat_brain.tiff", units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(brain_sampleTree, brain_traitColors,
groupLabels = names(brain_traitData[3:4]),
main = "Brain\nSample dendrogram and trait heatmap")
dev.off()
```

    ## png 
    ##   2

``` r
##HEART
#Re-cluster samples
heart_sampleTree = hclust(dist(heartDataExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
heart_traitColors = numbers2colors(heart_traitData[3:4], signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
tiff("sampleDendro_traitHeat_heart.tiff", units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(heart_sampleTree, heart_traitColors,
groupLabels = names(heart_traitData[3:4]),
main = "Heart\nSample dendrogram and trait heatmap")
dev.off()
```

    ## png 
    ##   2

``` r
##MUSCLE
#Re-cluster samples
muscle_sampleTree = hclust(dist(muscleDataExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
muscle_traitColors = numbers2colors(muscle_traitData[3:4], signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
tiff("sampleDendro_traitHeat_muscle.tiff", units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(muscle_sampleTree, muscle_traitColors,
groupLabels = names(muscle_traitData[3:4]),
main = "Muscle\nSample dendrogram and trait heatmap")
dev.off()
```

    ## png 
    ##   2

``` r
##LIVER
#Re-cluster samples
liver_sampleTree = hclust(dist(liverDataExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
liver_traitColors = numbers2colors(liver_traitData[3:4], signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
tiff("sampleDendro_traitHeat_liver.tiff", units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(liver_sampleTree, liver_traitColors,
groupLabels = names(liver_traitData[3:4]),
main = "Liver\nSample dendrogram and trait heatmap")
dev.off()
```

    ## png 
    ##   2

``` r
##PANCREAS
#Re-cluster samples
pancreas_sampleTree = hclust(dist(pancreasDataExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
pancreas_traitColors = numbers2colors(pancreas_traitData[3:4], signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
tiff("sampleDendro_traitHeat_pancreas.tiff", units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(pancreas_sampleTree, pancreas_traitColors,
groupLabels = names(pancreas_traitData[3:4]),
main = "Pancreas\nSample dendrogram and trait heatmap")
dev.off()
```

    ## png 
    ##   2

## 2.4 Save expression and trait data

This step is itended to save the relevant data files for use in the downstream analyses.

``` r
save(brainDataExpr, brain_traitData, file = "Brain-dataInput.RData")
save(heartDataExpr, heart_traitData, file = "Heart-dataInput.RData")
save(muscleDataExpr, muscle_traitData, file = "Muscle-dataInput.RData")
save(liverDataExpr, liver_traitData, file = "Liver-dataInput.RData")
save(pancreasDataExpr, pancreas_traitData, file = "Pancreas-dataInput.RData")
```

# 3. Gene network constructuin and identification of modules

Before proceeding to the network construction, one needs to choose *"the soft-thresholding power (β) to which co-expression similarity is raised to calculate adjacency."* \[@\] INSERT CITATION HERE.

## 3.1. Choice of soft-thresholding power: analysis of network topology

Zhang and Horvath \[-@\] have proposed the choice of this power β to be such that the network approximates a scale-free topology. The recommendations are that the scale-free topology fit is *"above 0.8 for reasonable powers (less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) "* and the mean connectivity is around [30-50](https://support.bioconductor.org/p/105337/). In case these criteria are not met due to an *"interesting biological variable"*, the choice of the soft-thresholding power can be performed based on the [number of samples and type of network](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html).

### Brain

#### Fig S6-Brain

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
brain_sft = pickSoftThreshold(brainDataExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
```

    ## pickSoftThreshold: will use block size 1309.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 1309 of 34164
    ##    ..working on genes 1310 through 2618 of 34164
    ##    ..working on genes 2619 through 3927 of 34164
    ##    ..working on genes 3928 through 5236 of 34164
    ##    ..working on genes 5237 through 6545 of 34164
    ##    ..working on genes 6546 through 7854 of 34164
    ##    ..working on genes 7855 through 9163 of 34164
    ##    ..working on genes 9164 through 10472 of 34164
    ##    ..working on genes 10473 through 11781 of 34164
    ##    ..working on genes 11782 through 13090 of 34164
    ##    ..working on genes 13091 through 14399 of 34164
    ##    ..working on genes 14400 through 15708 of 34164
    ##    ..working on genes 15709 through 17017 of 34164
    ##    ..working on genes 17018 through 18326 of 34164
    ##    ..working on genes 18327 through 19635 of 34164
    ##    ..working on genes 19636 through 20944 of 34164
    ##    ..working on genes 20945 through 22253 of 34164
    ##    ..working on genes 22254 through 23562 of 34164
    ##    ..working on genes 23563 through 24871 of 34164
    ##    ..working on genes 24872 through 26180 of 34164
    ##    ..working on genes 26181 through 27489 of 34164
    ##    ..working on genes 27490 through 28798 of 34164
    ##    ..working on genes 28799 through 30107 of 34164
    ##    ..working on genes 30108 through 31416 of 34164
    ##    ..working on genes 31417 through 32725 of 34164
    ##    ..working on genes 32726 through 34034 of 34164
    ##    ..working on genes 34035 through 34164 of 34164
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1    0.593 -1.96          0.743 4000.000  3.41e+03 9550.0
    ## 2      2    0.708 -2.16          0.931 1300.000  1.02e+03 5110.0
    ## 3      3    0.772 -2.24          0.977  513.000  3.60e+02 3010.0
    ## 4      4    0.821 -2.29          0.990  230.000  1.42e+02 1890.0
    ## 5      5    0.846 -2.35          0.987  112.000  5.99e+01 1240.0
    ## 6      6    0.870 -2.35          0.992   58.900  2.66e+01  845.0
    ## 7      7    0.886 -2.33          0.992   32.700  1.24e+01  592.0
    ## 8      8    0.893 -2.30          0.985   19.000  5.95e+00  425.0
    ## 9      9    0.888 -2.29          0.969   11.500  2.97e+00  311.0
    ## 10    10    0.901 -2.22          0.979    7.260  1.52e+00  233.0
    ## 11    12    0.895 -2.11          0.964    3.150  4.26e-01  136.0
    ## 12    14    0.918 -1.93          0.984    1.530  1.29e-01   83.6
    ## 13    16    0.944 -1.75          0.995    0.808  4.17e-02   53.7
    ## 14    18    0.944 -1.67          0.988    0.462  1.43e-02   40.6
    ## 15    20    0.945 -1.63          0.989    0.282  5.12e-03   33.5

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(brain_sft$fitIndices[,1], 
     -sign(brain_sft$fitIndices[,3])*brain_sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab=as.expression(bquote("Scale Free Topology Model Fit, signed " ~ R^2 ~ "")),
     type="n",
     main = paste("Scale independence"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(brain_sft$fitIndices[,1], 
     -sign(brain_sft$fitIndices[,3])*brain_sft$fitIndices[,2], 
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(brain_sft$fitIndices[,1],
     brain_sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n", 
     main = paste("Mean connectivity"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(brain_sft$fitIndices[,1], 
     brain_sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
# this line corresponds to a mean connectivty of 50
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
     line2user(line=5, side=3), 'Brain', xpd=NA, cex=2.5, font=2)


tiff('Fig-S7-brain.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5
plot(brain_sft$fitIndices[,1], -sign(brain_sft$fitIndices[,3])*brain_sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab=as.expression(bquote("Scale Free Topology Model Fit, signed " ~ R^2 ~ "")),type="n",main = paste("Scale independence"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(brain_sft$fitIndices[,1], -sign(brain_sft$fitIndices[,3])*brain_sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(brain_sft$fitIndices[,1], brain_sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(brain_sft$fitIndices[,1], brain_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),line2user(line=5, side=3), 'Brain', xpd=NA, cex=2.5, font=2)
dev.off()
```

    ## pdf 
    ##   2

power: 7 (mean connectivity between 30 and 50; scale free topology fit above 0.8)

### Heart

#### Fig S6-Heart

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
heart_sft = pickSoftThreshold(heartDataExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
```

    ## pickSoftThreshold: will use block size 1593.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 1593 of 28073
    ##    ..working on genes 1594 through 3186 of 28073
    ##    ..working on genes 3187 through 4779 of 28073
    ##    ..working on genes 4780 through 6372 of 28073
    ##    ..working on genes 6373 through 7965 of 28073
    ##    ..working on genes 7966 through 9558 of 28073
    ##    ..working on genes 9559 through 11151 of 28073
    ##    ..working on genes 11152 through 12744 of 28073
    ##    ..working on genes 12745 through 14337 of 28073
    ##    ..working on genes 14338 through 15930 of 28073
    ##    ..working on genes 15931 through 17523 of 28073
    ##    ..working on genes 17524 through 19116 of 28073
    ##    ..working on genes 19117 through 20709 of 28073
    ##    ..working on genes 20710 through 22302 of 28073
    ##    ..working on genes 22303 through 23895 of 28073
    ##    ..working on genes 23896 through 25488 of 28073
    ##    ..working on genes 25489 through 27081 of 28073
    ##    ..working on genes 27082 through 28073 of 28073
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1    0.438 -1.17          0.355 3870.000  3.59e+03 8960.0
    ## 2      2    0.645 -1.56          0.809 1400.000  9.92e+02 5190.0
    ## 3      3    0.726 -1.71          0.930  608.000  3.61e+02 3260.0
    ## 4      4    0.770 -1.83          0.968  294.000  1.47e+02 2160.0
    ## 5      5    0.810 -1.90          0.985  153.000  6.32e+01 1480.0
    ## 6      6    0.840 -1.94          0.993   85.000  2.86e+01 1040.0
    ## 7      7    0.863 -1.97          0.994   49.300  1.36e+01  755.0
    ## 8      8    0.878 -1.98          0.994   29.800  6.70e+00  558.0
    ## 9      9    0.891 -1.97          0.996   18.600  3.42e+00  420.0
    ## 10    10    0.899 -1.95          0.997   12.000  1.80e+00  321.0
    ## 11    12    0.915 -1.91          0.998    5.330  5.33e-01  195.0
    ## 12    14    0.922 -1.86          0.994    2.590  1.73e-01  124.0
    ## 13    16    0.928 -1.81          0.996    1.350  5.91e-02   82.2
    ## 14    18    0.928 -1.78          0.990    0.743  2.15e-02   56.6
    ## 15    20    0.930 -1.73          0.991    0.431  8.12e-03   40.0

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(heart_sft$fitIndices[,1], 
     -sign(heart_sft$fitIndices[,3])*heart_sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab=as.expression(bquote("Scale Free Topology Model Fit, signed " ~ R^2 ~ "")),
     type="n",
     main = paste("Scale independence"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(heart_sft$fitIndices[,1], 
     -sign(heart_sft$fitIndices[,3])*heart_sft$fitIndices[,2], 
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(heart_sft$fitIndices[,1],
     heart_sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n", 
     main = paste("Mean connectivity"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(heart_sft$fitIndices[,1], 
     heart_sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
# this line corresponds to a mean connectivty of 50
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
     line2user(line=5, side=3), 'Heart', xpd=NA, cex=2.5, font=2)


tiff('Fig-S7-heart.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5
plot(heart_sft$fitIndices[,1], -sign(heart_sft$fitIndices[,3])*heart_sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab=as.expression(bquote("Scale Free Topology Model Fit, signed " ~ R^2 ~ "")),type="n",main = paste("Scale independence"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(heart_sft$fitIndices[,1], -sign(heart_sft$fitIndices[,3])*heart_sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(heart_sft$fitIndices[,1], heart_sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(heart_sft$fitIndices[,1], heart_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),line2user(line=5, side=3), 'Heart', xpd=NA, cex=2.5, font=2)
dev.off()
```

    ## pdf 
    ##   2

power: 7 (mean connectivity between 30 and 50; scale free topology fit above 0.8)

### Muscle

#### Fig S6-Muscle

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
muscle_sft = pickSoftThreshold(muscleDataExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
```

    ## pickSoftThreshold: will use block size 2357.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2357 of 18977
    ##    ..working on genes 2358 through 4714 of 18977
    ##    ..working on genes 4715 through 7071 of 18977
    ##    ..working on genes 7072 through 9428 of 18977
    ##    ..working on genes 9429 through 11785 of 18977
    ##    ..working on genes 11786 through 14142 of 18977
    ##    ..working on genes 14143 through 16499 of 18977
    ##    ..working on genes 16500 through 18856 of 18977
    ##    ..working on genes 18857 through 18977 of 18977
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1    0.797 -1.62          0.747 1830.000  1.47e+03 4180.0
    ## 2      2    0.783 -1.67          0.721  580.000  3.57e+02 2400.0
    ## 3      3    0.716 -1.64          0.674  249.000  1.07e+02 1630.0
    ## 4      4    0.715 -1.58          0.729  129.000  3.74e+01 1200.0
    ## 5      5    0.741 -1.52          0.805   75.200  1.43e+01  920.0
    ## 6      6    0.740 -1.51          0.857   47.400  5.89e+00  724.0
    ## 7      7    0.742 -1.51          0.891   31.400  2.58e+00  579.0
    ## 8      8    0.743 -1.52          0.921   21.600  1.18e+00  470.0
    ## 9      9    0.724 -1.57          0.940   15.300  5.66e-01  385.0
    ## 10    10    0.745 -1.57          0.961   11.100  2.82e-01  318.0
    ## 11    12    0.773 -1.61          0.977    6.080  7.50e-02  222.0
    ## 12    14    0.797 -1.65          0.988    3.520  2.20e-02  159.0
    ## 13    16    0.822 -1.68          0.990    2.120  6.92e-03  116.0
    ## 14    18    0.844 -1.69          0.990    1.320  2.33e-03   85.5
    ## 15    20    0.861 -1.70          0.985    0.844  8.35e-04   64.1

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(muscle_sft$fitIndices[,1], 
     -sign(muscle_sft$fitIndices[,3])*muscle_sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),
     type="n",
     main = paste("Scale independence"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(muscle_sft$fitIndices[,1], 
     -sign(muscle_sft$fitIndices[,3])*muscle_sft$fitIndices[,2], 
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(muscle_sft$fitIndices[,1],
     muscle_sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n", 
     main = paste("Mean connectivity"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(muscle_sft$fitIndices[,1], 
     muscle_sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
# this line corresponds to a mean connectivty of 50
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
     line2user(line=5, side=3), 'Muscle', xpd=NA, cex=2.5, font=2)


tiff('Fig-S7-muscle.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5
plot(muscle_sft$fitIndices[,1], -sign(muscle_sft$fitIndices[,3])*muscle_sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),type="n",main = paste("Scale independence"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(muscle_sft$fitIndices[,1], -sign(muscle_sft$fitIndices[,3])*muscle_sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(muscle_sft$fitIndices[,1], muscle_sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(muscle_sft$fitIndices[,1], muscle_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),line2user(line=5, side=3), 'Muscle', xpd=NA, cex=2.5, font=2)
dev.off()
```

    ## pdf 
    ##   2

power: 6 (in this case most of the first 14 powers fail to reach 0.8 scale free topology fit) -&gt; mean connectivity around 50 is 6 which is in line with the choice based on the number of samples.

### Liver

#### Fig S6-Liver

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
liver_sft = pickSoftThreshold(liverDataExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
```

    ## pickSoftThreshold: will use block size 2219.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2219 of 20157
    ##    ..working on genes 2220 through 4438 of 20157
    ##    ..working on genes 4439 through 6657 of 20157
    ##    ..working on genes 6658 through 8876 of 20157
    ##    ..working on genes 8877 through 11095 of 20157
    ##    ..working on genes 11096 through 13314 of 20157
    ##    ..working on genes 13315 through 15533 of 20157
    ##    ..working on genes 15534 through 17752 of 20157
    ##    ..working on genes 17753 through 19971 of 20157
    ##    ..working on genes 19972 through 20157 of 20157
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1    0.531 -1.98          0.702 1860.000  1.62e+03 3890.0
    ## 2      2    0.755 -2.34          0.902  523.000  4.10e+02 1840.0
    ## 3      3    0.824 -2.43          0.958  188.000  1.29e+02 1040.0
    ## 4      4    0.849 -2.43          0.982   79.200  4.71e+01  633.0
    ## 5      5    0.860 -2.40          0.992   37.300  1.89e+01  407.0
    ## 6      6    0.871 -2.37          0.994   19.100  8.15e+00  271.0
    ## 7      7    0.883 -2.31          0.996   10.400  3.76e+00  186.0
    ## 8      8    0.885 -2.27          0.993    6.000  1.81e+00  130.0
    ## 9      9    0.888 -2.20          0.989    3.620  9.06e-01   93.0
    ## 10    10    0.888 -2.11          0.983    2.270  4.69e-01   67.7
    ## 11    12    0.953 -1.80          0.992    0.994  1.35e-01   37.4
    ## 12    14    0.976 -1.69          0.988    0.489  4.27e-02   27.2
    ## 13    16    0.980 -1.62          0.988    0.266  1.44e-02   21.3
    ## 14    18    0.977 -1.53          0.980    0.158  5.12e-03   16.9
    ## 15    20    0.979 -1.46          0.982    0.100  1.91e-03   13.6

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(liver_sft$fitIndices[,1], 
     -sign(liver_sft$fitIndices[,3])*liver_sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),
     type="n",
     main = paste("Scale independence"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(liver_sft$fitIndices[,1], 
     -sign(liver_sft$fitIndices[,3])*liver_sft$fitIndices[,2], 
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(liver_sft$fitIndices[,1],
     liver_sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n", 
     main = paste("Mean connectivity"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(liver_sft$fitIndices[,1], 
     liver_sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
# this line corresponds to a mean connectivty of 50
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
     line2user(line=5, side=3), 'Liver', xpd=NA, cex=2.5, font=2)


tiff('Fig-S7-liver.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5
plot(liver_sft$fitIndices[,1], -sign(liver_sft$fitIndices[,3])*liver_sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),type="n",main = paste("Scale independence"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(liver_sft$fitIndices[,1], -sign(liver_sft$fitIndices[,3])*liver_sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(liver_sft$fitIndices[,1], liver_sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(liver_sft$fitIndices[,1], liver_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),line2user(line=5, side=3), 'Liver', xpd=NA, cex=2.5, font=2)
dev.off()
```

    ## pdf 
    ##   2

power: 5 (mean connectivity between 30 and 50 and scale free tipology fit above 0.8)

### Pancreas

#### Fig S6-Pancreas

``` r
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
pancreas_sft = pickSoftThreshold(pancreasDataExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
```

    ## pickSoftThreshold: will use block size 2430.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2430 of 18411
    ##    ..working on genes 2431 through 4860 of 18411
    ##    ..working on genes 4861 through 7290 of 18411
    ##    ..working on genes 7291 through 9720 of 18411
    ##    ..working on genes 9721 through 12150 of 18411
    ##    ..working on genes 12151 through 14580 of 18411
    ##    ..working on genes 14581 through 17010 of 18411
    ##    ..working on genes 17011 through 18411 of 18411
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1    0.545 -1.62          0.613 2310.000  1.99e+03 5120.0
    ## 2      2    0.786 -1.54          0.860  802.000  5.99e+02 2670.0
    ## 3      3    0.808 -1.65          0.900  347.000  2.17e+02 1670.0
    ## 4      4    0.816 -1.72          0.922  173.000  8.60e+01 1140.0
    ## 5      5    0.829 -1.73          0.944   94.800  3.72e+01  816.0
    ## 6      6    0.823 -1.75          0.951   55.800  1.71e+01  604.0
    ## 7      7    0.831 -1.74          0.966   34.700  8.27e+00  458.0
    ## 8      8    0.839 -1.73          0.974   22.500  4.17e+00  354.0
    ## 9      9    0.854 -1.71          0.985   15.100  2.18e+00  278.0
    ## 10    10    0.852 -1.72          0.985   10.400  1.17e+00  221.0
    ## 11    12    0.862 -1.72          0.991    5.320  3.70e-01  144.0
    ## 12    14    0.863 -1.71          0.992    2.920  1.27e-01   97.3
    ## 13    16    0.868 -1.72          0.996    1.690  4.73e-02   67.9
    ## 14    18    0.877 -1.70          0.996    1.030  1.91e-02   48.5
    ## 15    20    0.881 -1.68          0.993    0.648  7.96e-03   35.4

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pancreas_sft$fitIndices[,1], 
     -sign(pancreas_sft$fitIndices[,3])*pancreas_sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
     ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),
     type="n",
     main = paste("Scale independence"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(pancreas_sft$fitIndices[,1], 
     -sign(pancreas_sft$fitIndices[,3])*pancreas_sft$fitIndices[,2], 
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pancreas_sft$fitIndices[,1],
     pancreas_sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n", 
     main = paste("Mean connectivity"), 
     cex.lab=1.6, 
     cex.axis=1.6, 
     cex.main=1.6)
text(pancreas_sft$fitIndices[,1], 
     pancreas_sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
# this line corresponds to a mean connectivty of 50
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
     line2user(line=5, side=3), 'Pancreas', xpd=NA, cex=2.5, font=2)


tiff('Fig-S7-pancreas.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
par(mfrow = c(1,2), mar=c(5,5,7,4))
cex1 = 1.5
plot(pancreas_sft$fitIndices[,1], -sign(pancreas_sft$fitIndices[,3])*pancreas_sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab=as.expression(bquote("Scale Free Topology Model Fit, signed" ~ R^2 ~ "")),type="n",main = paste("Scale independence"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(pancreas_sft$fitIndices[,1], -sign(pancreas_sft$fitIndices[,3])*pancreas_sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(pancreas_sft$fitIndices[,1], pancreas_sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"), cex.lab=1.6, cex.axis=1.6, cex.main=1.6)
text(pancreas_sft$fitIndices[,1], pancreas_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=50,col="red")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),line2user(line=5, side=3), 'Pancreas', xpd=NA, cex=2.5, font=2)
dev.off()
```

    ## pdf 
    ##   2

power: 7 (mean connectivity between 30 and 50; scale free topology fit above 0.8; in line with choice based on the sample number)

## 3.2. Block-wise network construction and module detection

Throughout this tutorial we work with a relatively small data set of 3600 measured probes. However, modern microarrays measure up to 50,000 probe expression levels at once. Constructing and analyzing networks with such large numbers of nodes is computationally challenging even on a large server. We now illustrate a method, implemented in the WGCNA package, that allows the user to perform a network analysis with such a large number of genes. Instead of actually using a very large data set, we will for simplicity pretend that hardware limitations restrict the number of genes that can be analyzed at once to 2000. The basic idea is to use a two-level clustering. First, we use a fast, computationally inexpensive and relatively crude clustering method to pre-cluster genes into blocks of size close to and not exceeding the maximum of 2000 genes. We then perform a full network analysis in each block separately. At the end, modules whose eigengenes are highly correlated are merged. The advantage of the block-wise approach is a much smaller memory footprint (which is the main problem with large data sets on standard desktop computers), and a significant speed-up of the calculations. The trade-off is that due to using a simpler clustering to obtain blocks, the blocks may not be optimal, causing some outlying genes to be assigned to a different module than they would be in a full network analysis. We will now pretend that even the relatively small number of genes, 3600, that we have been using here is too large, and the computer we run the analysis on is not capable of handling more than 2000 genes in one block.

The automatic network construction and module detection function blockwiseModules can handle the splitting into blocks automatically; the user just needs to specify the largest number of genes that can fit in a block:

We have chosen the soft thresholding power 6, a relatively large minimum module size of 30, and a medium sensitivity (deepSplit=2) to cluster splitting. The parameter `mergeCutHeight` is the threshold for merging of modules. We have also instructed the function to return numeric, rather than color, labels for modules, and to save the Topological Overlap Matrix. The output of the function may seem somewhat cryptic, but it is easy to use. For example, `bwnet$colors` contains the module assignment, and `bwnet$MEs` contains the module eigengenes of the modules.

A word of caution for the readers who would like to adapt this code for their own data. The function blockwiseModules has many parameters, and in this example most of them are left at their default value. We have attempted to provide reasonable default values, but they may not be appropriate for the particular data set the reader wishes to analyze. We encourage the user to read the help file provided within the package in the R environment and experiment with tweaking the network construction and module detection parameters. The potential reward is, of course, better (biologically more relevant) results of the analysis.

A second word of caution concerning block size. In particular, the parameter maxBlockSize tells the function how large the largest block can be that the reader’s computer can handle. In this example we have set the maximum block size to 2000 to illustrate the block-wise analysis and its results, but this value is needlessly small for most modern computers; the default is 5000 which is appropriate for most modern desktops. If the reader has access to a large workstation with more than 4 GB of memory, the parameter maxBlockSize can be increased. A 16GB workstation should handle up to 20000 probes; a 32GB workstation should handle perhaps 30000. A 4GB standard desktop or a laptop may handle up to 8000-10000 probes, depending on operating system and other running programs. In general it is preferable to analyze a data set in as few blocks as possible. Below we will compare the results of this analysis to the results of Section 2.a in which all genes were analyzed in a single block. To make the comparison easier, we relabel the block-wise module labels so that modules with a significant overlap with single-block modules have the same label:

### Brain

``` r
bwnet_brain = blockwiseModules(brainDataExpr, 
                               corType = "bicor",
                               maxPOutliers = 0.05,
                               power = 7, 
                               networkType = "signed hybrid", 
                               minModuleSize = 50,
                               numericLabels = TRUE, 
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "brainTOM-blockwise",
                               verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4    5    6    7 
    ## 4998 4987 4986 4978 4977 4643 4595 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file brainTOM-blockwise-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 2 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 3 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file brainTOM-blockwise-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 64 genes from module 1 because their KME is too low.
    ##  ..Working on block 3 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 3 into file brainTOM-blockwise-block.3.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 3 genes from module 1 because their KME is too low.
    ##      ..removing 3 genes from module 2 because their KME is too low.
    ##  ..Working on block 4 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 4 into file brainTOM-blockwise-block.4.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 13 genes from module 1 because their KME is too low.
    ##  ..Working on block 5 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 5 into file brainTOM-blockwise-block.5.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 9 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 6 because their KME is too low.
    ##      ..removing 1 genes from module 8 because their KME is too low.
    ##  ..Working on block 6 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 6 into file brainTOM-blockwise-block.6.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 7 genes from module 1 because their KME is too low.
    ##  ..Working on block 7 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 7 into file brainTOM-blockwise-block.7.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 1 genes from module 1 because their KME is too low.
    ##      ..removing 3 genes from module 2 because their KME is too low.
    ##   ..reassigning 5 genes from module 2 to modules with higher KME.
    ##   ..reassigning 4 genes from module 4 to modules with higher KME.
    ##   ..reassigning 16 genes from module 6 to modules with higher KME.
    ##   ..reassigning 6 genes from module 8 to modules with higher KME.
    ##   ..reassigning 3 genes from module 10 to modules with higher KME.
    ##   ..reassigning 3 genes from module 11 to modules with higher KME.
    ##   ..reassigning 1 genes from module 12 to modules with higher KME.
    ##   ..reassigning 3 genes from module 14 to modules with higher KME.
    ##   ..reassigning 1 genes from module 15 to modules with higher KME.
    ##   ..reassigning 2 genes from module 22 to modules with higher KME.
    ##   ..reassigning 1 genes from module 23 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.15
    ##        Calculating new MEs...

``` r
brain_moduleLabels = bwnet_brain$colors
brain_moduleColors = labels2colors(bwnet_brain$colors)
brain_MEs = bwnet_brain$MEs;
brain_geneTree1 = bwnet_brain$dendrograms[[1]];
brain_geneTree2 = bwnet_brain$dendrograms[[2]];
brain_geneTree3 = bwnet_brain$dendrograms[[3]];
brain_geneTree4 = bwnet_brain$dendrograms[[4]];
brain_geneTree5 = bwnet_brain$dendrograms[[5]];
brain_geneTree6 = bwnet_brain$dendrograms[[6]];
brain_geneTree7 = bwnet_brain$dendrograms[[7]];
brain_blockgenes1 = bwnet_brain$blockGenes[[1]];
brain_blockgenes2 = bwnet_brain$blockGenes[[2]];
brain_blockgenes3 = bwnet_brain$blockGenes[[3]];
brain_blockgenes4 = bwnet_brain$blockGenes[[4]];
brain_blockgenes5 = bwnet_brain$blockGenes[[5]];
brain_blockgenes6 = bwnet_brain$blockGenes[[6]];
brain_blockgenes7 = bwnet_brain$blockGenes[[7]];
save(brain_MEs, brain_moduleLabels, brain_moduleColors, brain_geneTree1, brain_geneTree2, brain_geneTree3, brain_geneTree4, brain_geneTree5, brain_geneTree6, brain_geneTree7, brain_blockgenes1, brain_blockgenes2,  brain_blockgenes3,  brain_blockgenes4,  brain_blockgenes5,  brain_blockgenes6, brain_blockgenes7,
     file = "Brain-networkConstruction-auto.RData")
```

### Heart

``` r
bwnet_heart = blockwiseModules(heartDataExpr, 
                               corType = "bicor",
                               maxPOutliers = 0.05,
                               power = 7, 
                               networkType = "signed hybrid", 
                               minModuleSize = 50,
                               numericLabels = TRUE, 
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "heartTOM-blockwise",
                               verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4    5    6 
    ## 5000 4993 4986 4980 4376 3738 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file heartTOM-blockwise-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 21 genes from module 1 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file heartTOM-blockwise-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 14 genes from module 1 because their KME is too low.
    ##  ..Working on block 3 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 3 into file heartTOM-blockwise-block.3.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 2 genes from module 1 because their KME is too low.
    ##  ..Working on block 4 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 4 into file heartTOM-blockwise-block.4.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 19 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 3 because their KME is too low.
    ##  ..Working on block 5 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 5 into file heartTOM-blockwise-block.5.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 16 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 2 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##      ..removing 1 genes from module 5 because their KME is too low.
    ##      ..removing 1 genes from module 6 because their KME is too low.
    ##  ..Working on block 6 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 6 into file heartTOM-blockwise-block.6.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 3 genes from module 2 because their KME is too low.
    ##      ..removing 2 genes from module 3 because their KME is too low.
    ##      ..removing 3 genes from module 4 because their KME is too low.
    ##   ..reassigning 2 genes from module 1 to modules with higher KME.
    ##   ..reassigning 1 genes from module 4 to modules with higher KME.
    ##   ..reassigning 2 genes from module 6 to modules with higher KME.
    ##   ..reassigning 2 genes from module 7 to modules with higher KME.
    ##   ..reassigning 4 genes from module 8 to modules with higher KME.
    ##   ..reassigning 1 genes from module 9 to modules with higher KME.
    ##   ..reassigning 1 genes from module 10 to modules with higher KME.
    ##   ..reassigning 1 genes from module 13 to modules with higher KME.
    ##   ..reassigning 1 genes from module 14 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.15
    ##        Calculating new MEs...

``` r
heart_moduleLabels = bwnet_heart$colors
heart_moduleColors = labels2colors(bwnet_heart$colors)
heart_MEs = bwnet_heart$MEs;
heart_geneTree1 = bwnet_heart$dendrograms[[1]];
heart_geneTree2 = bwnet_heart$dendrograms[[2]];
heart_geneTree3 = bwnet_heart$dendrograms[[3]];
heart_geneTree4 = bwnet_heart$dendrograms[[4]];
heart_geneTree5 = bwnet_heart$dendrograms[[5]];
heart_geneTree6 = bwnet_heart$dendrograms[[6]];
heart_blockgenes1 = bwnet_heart$blockGenes[[1]];
heart_blockgenes2 = bwnet_heart$blockGenes[[2]];
heart_blockgenes3 = bwnet_heart$blockGenes[[3]];
heart_blockgenes4 = bwnet_heart$blockGenes[[4]];
heart_blockgenes5 = bwnet_heart$blockGenes[[5]];
heart_blockgenes6 = bwnet_heart$blockGenes[[6]];
save(heart_MEs, heart_moduleLabels, heart_moduleColors, heart_geneTree1, heart_geneTree2, heart_geneTree3, heart_geneTree4, heart_geneTree5, heart_geneTree6, heart_blockgenes1, heart_blockgenes2, heart_blockgenes3, heart_blockgenes4, heart_blockgenes5, heart_blockgenes6,
     file = "Heart-networkConstruction-auto.RData")
```

### Muscle

``` r
bwnet_muscle = blockwiseModules(muscleDataExpr, 
                               corType = "bicor",
                               maxPOutliers = 0.05,
                               power = 6, 
                               networkType = "signed hybrid", 
                               minModuleSize = 50,
                               numericLabels = TRUE, 
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "muscleTOM-blockwise",
                               verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4 
    ## 5000 4859 4749 4369 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file muscleTOM-blockwise-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 297 genes from module 1 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file muscleTOM-blockwise-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 22 genes from module 1 because their KME is too low.
    ##      ..removing 25 genes from module 2 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##      ..removing 9 genes from module 5 because their KME is too low.
    ##      ..removing 1 genes from module 7 because their KME is too low.
    ##  ..Working on block 3 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 3 into file muscleTOM-blockwise-block.3.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 75 genes from module 1 because their KME is too low.
    ##      ..removing 53 genes from module 2 because their KME is too low.
    ##      ..removing 3 genes from module 3 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##  ..Working on block 4 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 4 into file muscleTOM-blockwise-block.4.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 32 genes from module 1 because their KME is too low.
    ##      ..removing 5 genes from module 2 because their KME is too low.
    ##      ..removing 2 genes from module 4 because their KME is too low.
    ##      ..removing 1 genes from module 5 because their KME is too low.
    ##      ..removing 2 genes from module 6 because their KME is too low.
    ##   ..reassigning 1 genes from module 4 to modules with higher KME.
    ##   ..reassigning 1 genes from module 11 to modules with higher KME.
    ##   ..reassigning 2 genes from module 15 to modules with higher KME.
    ##   ..reassigning 1 genes from module 16 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.15
    ##        Calculating new MEs...

``` r
muscle_moduleLabels = bwnet_muscle$colors
muscle_moduleColors = labels2colors(bwnet_muscle$colors)
muscle_MEs = bwnet_muscle$MEs;
muscle_geneTree1 = bwnet_muscle$dendrograms[[1]];
muscle_geneTree2 = bwnet_muscle$dendrograms[[2]];
muscle_geneTree3 = bwnet_muscle$dendrograms[[3]];
muscle_geneTree4 = bwnet_muscle$dendrograms[[4]];
muscle_blockgenes1 = bwnet_muscle$blockGenes[[1]];
muscle_blockgenes2 = bwnet_muscle$blockGenes[[2]];
muscle_blockgenes3 = bwnet_muscle$blockGenes[[3]];
muscle_blockgenes4 = bwnet_muscle$blockGenes[[4]];
save(muscle_MEs, muscle_moduleLabels, muscle_moduleColors, muscle_geneTree1, muscle_geneTree2, muscle_geneTree3, muscle_geneTree4, muscle_blockgenes1, muscle_blockgenes2,muscle_blockgenes3, muscle_blockgenes4,
     file = "Muscle-networkConstruction-auto.RData")
```

### Liver

``` r
bwnet_liver = blockwiseModules(liverDataExpr, 
                               corType = "bicor",
                               maxPOutliers = 0.05,
                               power = 5, 
                               networkType = "signed hybrid", 
                               minModuleSize = 50,
                               numericLabels = TRUE, 
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "liverTOM-blockwise",
                               verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4    5 
    ## 4997 4993 4719 3651 1797 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file liverTOM-blockwise-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 21 genes from module 1 because their KME is too low.
    ##      ..removing 17 genes from module 2 because their KME is too low.
    ##      ..removing 6 genes from module 3 because their KME is too low.
    ##      ..removing 7 genes from module 4 because their KME is too low.
    ##      ..removing 2 genes from module 5 because their KME is too low.
    ##      ..removing 5 genes from module 6 because their KME is too low.
    ##      ..removing 3 genes from module 7 because their KME is too low.
    ##      ..removing 2 genes from module 8 because their KME is too low.
    ##      ..removing 1 genes from module 9 because their KME is too low.
    ##      ..removing 1 genes from module 11 because their KME is too low.
    ##      ..removing 1 genes from module 12 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file liverTOM-blockwise-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 288 genes from module 1 because their KME is too low.
    ##  ..Working on block 3 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 3 into file liverTOM-blockwise-block.3.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 159 genes from module 1 because their KME is too low.
    ##      ..removing 117 genes from module 2 because their KME is too low.
    ##      ..removing 65 genes from module 3 because their KME is too low.
    ##      ..removing 14 genes from module 4 because their KME is too low.
    ##      ..removing 13 genes from module 5 because their KME is too low.
    ##      ..removing 4 genes from module 6 because their KME is too low.
    ##      ..removing 3 genes from module 7 because their KME is too low.
    ##      ..removing 7 genes from module 8 because their KME is too low.
    ##      ..removing 1 genes from module 9 because their KME is too low.
    ##      ..removing 1 genes from module 10 because their KME is too low.
    ##  ..Working on block 4 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 4 into file liverTOM-blockwise-block.4.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 29 genes from module 1 because their KME is too low.
    ##      ..removing 23 genes from module 2 because their KME is too low.
    ##      ..removing 18 genes from module 3 because their KME is too low.
    ##      ..removing 7 genes from module 4 because their KME is too low.
    ##      ..removing 11 genes from module 5 because their KME is too low.
    ##      ..removing 6 genes from module 6 because their KME is too low.
    ##      ..removing 8 genes from module 7 because their KME is too low.
    ##      ..removing 1 genes from module 8 because their KME is too low.
    ##      ..removing 3 genes from module 9 because their KME is too low.
    ##      ..removing 5 genes from module 10 because their KME is too low.
    ##      ..removing 3 genes from module 11 because their KME is too low.
    ##      ..removing 2 genes from module 12 because their KME is too low.
    ##  ..Working on block 5 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 5 into file liverTOM-blockwise-block.5.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 107 genes from module 1 because their KME is too low.
    ##      ..removing 136 genes from module 2 because their KME is too low.
    ##   ..reassigning 1 genes from module 1 to modules with higher KME.
    ##   ..reassigning 9 genes from module 14 to modules with higher KME.
    ##   ..reassigning 3 genes from module 15 to modules with higher KME.
    ##   ..reassigning 5 genes from module 16 to modules with higher KME.
    ##   ..reassigning 3 genes from module 18 to modules with higher KME.
    ##   ..reassigning 1 genes from module 21 to modules with higher KME.
    ##   ..reassigning 2 genes from module 29 to modules with higher KME.
    ##   ..reassigning 1 genes from module 31 to modules with higher KME.
    ##   ..reassigning 1 genes from module 34 to modules with higher KME.
    ##   ..reassigning 8 genes from module 37 to modules with higher KME.
    ##   ..reassigning 14 genes from module 38 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.15
    ##        Calculating new MEs...

``` r
liver_moduleLabels = bwnet_liver$colors
liver_moduleColors = labels2colors(bwnet_liver$colors)
liver_MEs = bwnet_liver$MEs;
liver_geneTree1 = bwnet_liver$dendrograms[[1]];
liver_geneTree2 = bwnet_liver$dendrograms[[2]];
liver_geneTree3 = bwnet_liver$dendrograms[[3]];
liver_geneTree4 = bwnet_liver$dendrograms[[4]];
liver_geneTree5 = bwnet_liver$dendrograms[[5]];
liver_blockgenes1 = bwnet_liver$blockGenes[[1]];
liver_blockgenes2 = bwnet_liver$blockGenes[[2]];
liver_blockgenes3 = bwnet_liver$blockGenes[[3]];
liver_blockgenes4 = bwnet_liver$blockGenes[[4]];
liver_blockgenes5 = bwnet_liver$blockGenes[[5]];
save(liver_MEs, liver_moduleLabels, liver_moduleColors, liver_geneTree1, liver_geneTree2, liver_geneTree3, liver_geneTree4, liver_geneTree5, liver_blockgenes1, liver_blockgenes2,liver_blockgenes3, liver_blockgenes4,liver_blockgenes5, 
     file = "Liver-networkConstruction-auto.RData")
```

### Pancreas

``` r
bwnet_pancreas = blockwiseModules(pancreasDataExpr, 
                               corType = "bicor",
                               maxPOutliers = 0.05,
                               power = 7, 
                               networkType = "signed hybrid", 
                               minModuleSize = 50,
                               numericLabels = TRUE, 
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "pancreasTOM-blockwise",
                               verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4 
    ## 4994 4990 4832 3595 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file pancreasTOM-blockwise-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 69 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 2 because their KME is too low.
    ##      ..removing 3 genes from module 3 because their KME is too low.
    ##      ..removing 5 genes from module 4 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file pancreasTOM-blockwise-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 77 genes from module 1 because their KME is too low.
    ##  ..Working on block 3 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 3 into file pancreasTOM-blockwise-block.3.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##  ..Working on block 4 .
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 4 into file pancreasTOM-blockwise-block.4.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 9 genes from module 1 because their KME is too low.
    ##      ..removing 1 genes from module 3 because their KME is too low.
    ##   ..reassigning 4 genes from module 1 to modules with higher KME.
    ##   ..reassigning 2 genes from module 7 to modules with higher KME.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.15
    ##        Calculating new MEs...

``` r
pancreas_moduleLabels = bwnet_pancreas$colors
pancreas_moduleColors = labels2colors(bwnet_pancreas$colors)
pancreas_MEs = bwnet_pancreas$MEs;
pancreas_geneTree1 = bwnet_pancreas$dendrograms[[1]];
pancreas_geneTree2 = bwnet_pancreas$dendrograms[[2]];
pancreas_geneTree3 = bwnet_pancreas$dendrograms[[3]];
pancreas_geneTree4 = bwnet_pancreas$dendrograms[[4]];
pancreas_blockgenes1 = bwnet_pancreas$blockGenes[[1]];
pancreas_blockgenes2 = bwnet_pancreas$blockGenes[[2]];
pancreas_blockgenes3 = bwnet_pancreas$blockGenes[[3]];
pancreas_blockgenes4 = bwnet_pancreas$blockGenes[[4]];
save(pancreas_MEs, pancreas_moduleLabels, pancreas_moduleColors, pancreas_geneTree1, pancreas_geneTree2,pancreas_geneTree3, pancreas_geneTree4, pancreas_blockgenes1, pancreas_blockgenes2, pancreas_blockgenes3, pancreas_blockgenes4,
     file = "Pancreas-networkConstruction-auto.RData")
```

Langfelder, Peter, and Steve Horvath. 2008. “WGCNA: an R package for weighted correlation network analysis.” doi:[10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559).
