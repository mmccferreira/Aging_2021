Integration of segmented regression analysis with weighted gene correlation network analysis identifies genes whose expression is remodeled throughout physiological aging in mouse tissues
================
Margarida Ferreira
Last updated July 05, 2021

-   [1. Initial setup](#initial-setup)
-   [2. Data input](#data-input)
    -   [2.1. Removal of low coverage samples](#removal-of-low-coverage-samples)
    -   [2.2. Annotation](#annotation)
    -   [2.3. DESeqDataSet object](#deseqdataset-object)
    -   [2.3. DESeqDataSet subsetting \#1](#deseqdataset-subsetting-1)
-   [3. Exploratory analysis and visualization](#exploratory-analysis-and-visualization)
    -   [3.1. Principal component analysis - All studied tissues together](#principal-component-analysis---all-studied-tissues-together)
    -   [3.2. DESeqDataSet subsetting \#2](#deseqdataset-subsetting-2)
    -   [3.2. Filtering of low expressed genes](#filtering-of-low-expressed-genes)
    -   [3.3. Identification and removal of outlier samples](#identification-and-removal-of-outlier-samples)

## 1. Initial setup

Loading the required packages:

``` r
library(DESeq2)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(grid)
library(tidyverse)
library(plyr)
library(cowplot)
library(ggpubr)
```

## 2. Data input

We used DESeq2 (Love, Huber, and Anders 2014) to perform the filtering and normalization steps. As input, the DESeq2 package expects count data in the form of a matrix of integer values, with samples in columns and features (i.e. genes/transcripts) in rows. We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata. The .csv files we are uploading were downloaded from the [GSE132040](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132040) web page and extracted to the folder in which this script is located.

``` r
##read counts
cts <- as.matrix(read.csv("GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv", header=TRUE, row.names = "gene"))
#remove ".gencode.vM19" expression from samples names in the count matrix to match the sample names provided in the metadata file
colnames(cts) = gsub(pattern = ".gencode.vM19", replacement = "", x = colnames(cts))
#remove the last 5 rows because they don't contain read counts
cts <- cts[-c(54353:54357),]
dim(cts)
```

    ## [1] 54352   947

``` r
##sample information
coldata <- read.csv("GSE132040_MACA_Bulk_metadata.csv", header=TRUE)#, row.names=1)

coldata[c(2,4,6,8:13)] <- NULL
colnames(coldata) <- c("Sample","Tissue","Age","Sex")
rownames(coldata) <- coldata$Sample

#removal of "_" for simplifcation
coldata$Tissue = gsub(pattern = "Limb_Muscle_.*", replacement = "Muscle", x = coldata$Tissue)
coldata$Tissue = gsub(pattern = "Small_Intestine_.*", replacement = "Intestine", x = coldata$Tissue)
coldata$Tissue = gsub(pattern = "_.*", replacement = "", x = coldata$Tissue)

#replacement of "f" and "m by "Female" and "Male" for ease of understanding
coldata$Sex = gsub(pattern = "m", replacement = "Male", x = coldata$Sex)
coldata$Sex = gsub(pattern = "f", replacement = "Female", x = coldata$Sex)


dim(coldata)
```

    ## [1] 947   4

### 2.1. Removal of low coverage samples

This step was performed as indicated by the authors of the original dataset (Schaum et al. 2020).

``` r
#identify which samples have a sequencing coverage of less than 4 million reads
low_coverage <- colnames(cts[,colSums(cts) < 4000000])
length(low_coverage)
```

    ## [1] 115

``` r
#keep track of the removed samples by saving the removed samples' names in a .txt file
write.table(low_coverage, "low_coverage_samples.txt", quote=FALSE)

#remove those samples from the count matrix
cts <- cts[,!colnames(cts) %in% low_coverage]
ncol(cts)
```

    ## [1] 832

``` r
#remove those samples from the sample information table
coldata <- coldata[!rownames(coldata) %in% low_coverage,]
nrow(coldata)
```

    ## [1] 832

### 2.2. Annotation

The next step uses the R package biomaRt (Durinck et al. 2005) for cross-referencing gene symbols with Ensembl biotype annotations.

These are the details of the annotation:

``` r
mouseMart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
datasets <- listDatasets(mouseMart)
#check the current version of the annotation dataset
datasets[datasets$dataset=="mmusculus_gene_ensembl",]
```

    ##                    dataset          description version
    ## 107 mmusculus_gene_ensembl Mouse genes (GRCm39)  GRCm39

``` r
# #cross-reference gene symbols and gene biotypes
# annotation <- getBM(attributes = c("external_gene_name","gene_biotype"), filters="external_gene_name", values=rownames(cts), mart=mouseMart)
# colnames(annotation) <- c("SYMBOL","BIOTYPE")
# dim(annotation)
```

``` r
# #duplicate symbol removal
# dup_symbols <- duplicated(annotation$SYMBOL)
# table(dup_symbols)
# #12 duplicate gene names were found, comprising a total of 24 entries to remove
# dupsToRemove <- annotation[dup_symbols,]$SYMBOL
# annotation <- annotation[!annotation$SYMBOL %in% dupsToRemove,]
# dim(annotation)
# 
# write.table(annotation, "annotation.txt", row.names = F, quote=F)

#annotation downloaded in May 17th 2021
annotation <- read.table("annotation.txt", header = T)
```

Now we will remove from the read count matrix the genes without biotype annotation.

``` r
cts <- as.matrix(cts[rownames(cts) %in% annotation$SYMBOL,])
dim(cts)
```

    ## [1] 49512   832

### 2.3. DESeqDataSet object

We will use the R package DESeq2 to normalize the data so that the gene expression values can be compared between samples.

First, we will make sure the count matrix and column data to see if they are consistent in terms of sample order.

As stated in [DESeq2's vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html):

*"It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order."*

``` r
cts <- cts[,order(colnames(cts), decreasing=FALSE)]

coldata <- coldata[order(rownames(coldata), decreasing=FALSE),]

all(rownames(coldata) %in% colnames(cts))
```

    ## [1] TRUE

``` r
#if TRUE, then rownames(coldata) and colnames(cts) are the same.
all(rownames(coldata) == colnames(cts))
```

    ## [1] TRUE

``` r
#if TRUE, then rownames(coldata) and colnames(cts) are in the same order.
```

As we can see, the sample order is the same in the `cts` matrix and in the `coldata` table.

Now we can construct the `DESeqDataSet` object:

``` r
dds <- DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=~Sex + Age)
dds
```

    ## class: DESeqDataSet 
    ## dim: 49512 832 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(49512): 0610005C13Rik 0610006L08Rik ... n-R5s85 n-R5s88
    ## rowData names(0):
    ## colnames(832): A1_384Bulk_Plate1_S1 A1_384Bulk_Plate3_S1 ...
    ##   P9_384Bulk_Plate2_S369 P9_384Bulk_Plate3_S369
    ## colData names(4): Sample Tissue Age Sex

### 2.3. DESeqDataSet subsetting \#1

For this study, we excluded the 1-month-old samples to avoid the influence of developmental genes and selected the brain, heart, muscle, liver, and pancreas for further analyses.

``` r
dim(dds) #dimensions of count matrix before the removal of 1-month samples
```

    ## [1] 49512   832

``` r
dds <- subset(dds, select=dds$Age!="1") #remove 1-month samples
dds <- subset(dds, select=dds$Tissue %in% c("Brain","Heart","Liver","Muscle","Pancreas")) #only keep samples from brain, heart, muscle, liver and pancreas
dim(dds) #dimensions of count matrix after the removal of 1-month samples
```

    ## [1] 49512   222

``` r
#Remove the 1-month samples from the sample information table
coldata <- coldata[rownames(coldata) %in% colnames(dds) ,]
dim(coldata)
```

    ## [1] 222   4

``` r
#Drop the unused levels from the factors included in the design
dds$Age <- droplevels(dds$Age)
dds$Age <- factor(dds$Age, levels=c(3,6,9,12,15,18,21,24,27)) #reorder the levels
dds$Sex <- droplevels(dds$Sex)
dim(dds)
```

    ## [1] 49512   222

## 3. Exploratory analysis and visualization

The next steps involve the transformation of the counts in order to stabilize the variances across different mean values and visually explore sample relationships.

We used the *variance stabilizing transformation* implemented in DESeq2's `vst` function.

``` r
vst <- vst(dds, blind = FALSE)
```

#### 3.1. Principal component analysis - All studied tissues together

##### Fig 1A

``` r
object=vst

intgroup = "Tissue"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar1 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d1 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d1, "percentVar") <- percentVar1[1:2]
        return(d1)
    }

intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar1 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d2 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d2, "percentVar") <- percentVar1[1:2]
        return(d2)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar2 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d3 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d3, "percentVar") <- percentVar2[1:2]
        return(d3)
    }


d1$Effect <- "Tissue"
d2$Effect <- "Age"
d3$Effect <- "Sex"


name <- c(d1$name,d2$name,d3$name)
PC1 <- c(d1$PC1,d2$PC1,d3$PC1)
PC2 <- c(d1$PC2,d2$PC2,d3$PC2)
group <- c(as.character(d1$group),as.character(d2$group),as.character(d3$group))
effect <- c(d1$Effect,d2$Effect,d3$Effect)

d = data.frame(name,PC1,PC2,group,effect)


d$group <- factor(d$group, levels = c("Brain","Heart","Liver","Muscle","Pancreas","3","6","9","12","15","18","21","24","27","Female","Male"))
d$effect <- factor(d$effect, levels = c("Tissue","Age","Sex"))

mycolors <- c('#77AADD', '#EE8866','#EEDD88','#FFAABB', '#44BB99','#e6194B', '#f58231', '#ffe119', '#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#F8766D','#619CFF')


a <- ggplot(data = d, aes_string(x = "PC1", y = "PC2",
      color = "group")) + geom_point(size = 3) + xlab(paste0("PC1: ",
        round(percentVar1[1] * 100), "% variance")) + ylab(paste0("PC2: ",
        round(percentVar1[2] * 100), "% variance")) +
      theme_light() + scale_color_manual(values = mycolors) +
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
            axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            legend.title =element_blank(), legend.text = element_text(size=20),
            axis.text.x = element_text(hjust = 0.5, face="bold"),
            strip.text = element_text(size=25, face="bold")) +
  facet_grid(rows=vars(effect), scales = "free") + guides(color=guide_legend(ncol=1))

a
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
tiff('Fig-1A.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
a
dev.off()
```

    ## png 
    ##   2

Because the samples segregate mainly by tissue type, further analyses will be performed on each tissue independently.

### 3.2. DESeqDataSet subsetting \#2

First we re-subset the `dds` object for each tissue.

##### Brain

``` r
#Subset the count matrix to include only the brain samples
dds_brain <- subset(dds, select=dds$Tissue=="Brain")
dim(dds_brain)
```

    ## [1] 49512    48

``` r
#Create a sample information table comprising only the information of the brain samples
coldata_brain <- coldata[rownames(coldata) %in% colnames(dds_brain) ,]
dim(coldata_brain)
```

    ## [1] 48  4

``` r
#Check the sample order
all(rownames(coldata_brain) %in% colnames(dds_brain))
```

    ## [1] TRUE

``` r
all(rownames(coldata_brain) == colnames(dds_brain))
```

    ## [1] TRUE

``` r
#Drop the unused levels
dds_brain$Age <- droplevels(dds_brain$Age)
dds_brain$Age <- factor(dds_brain$Age, levels=c(3,6,9,12,15,18,21,24,27))
dds_brain$Sex <- droplevels(dds_brain$Sex)
dim(dds_brain)
```

    ## [1] 49512    48

##### Heart

``` r
#Subset the count matrix to include only the heart samples
dds_heart <- subset(dds, select=dds$Tissue=="Heart")
dim(dds_heart)
```

    ## [1] 49512    47

``` r
#Create a sample information table comprising only the information of the heart samples
coldata_heart <- coldata[rownames(coldata) %in% colnames(dds_heart) ,]
dim(coldata_heart)
```

    ## [1] 47  4

``` r
#Check the sample order
all(rownames(coldata_heart) %in% colnames(dds_heart))
```

    ## [1] TRUE

``` r
all(rownames(coldata_heart) == colnames(dds_heart))
```

    ## [1] TRUE

``` r
#Drop the unused levels
dds_heart$Age <- droplevels(dds_heart$Age)
dds_heart$Age <- factor(dds_heart$Age, levels=c(3,6,9,12,15,18,21,24,27))
dds_heart$Sex <- droplevels(dds_heart$Sex)
dim(dds_heart)
```

    ## [1] 49512    47

##### Muscle

``` r
#Subset the count matrix to include only the muscle samples
dds_muscle <- subset(dds, select=dds$Tissue=="Muscle")
dim(dds_muscle)
```

    ## [1] 49512    47

``` r
#Create a sample information table comprising only the information of the muscle samples
coldata_muscle <- coldata[rownames(coldata) %in% colnames(dds_muscle) ,]
dim(coldata_muscle)
```

    ## [1] 47  4

``` r
#Check the sample order
all(rownames(coldata_muscle) %in% colnames(dds_muscle))
```

    ## [1] TRUE

``` r
all(rownames(coldata_muscle) == colnames(dds_muscle))
```

    ## [1] TRUE

``` r
#Drop the unused levels
dds_muscle$Age <- droplevels(dds_muscle$Age)
dds_muscle$Age <- factor(dds_muscle$Age, levels=c(3,6,9,12,15,18,21,24,27))
dds_muscle$Sex <- droplevels(dds_muscle$Sex)
dim(dds_muscle)
```

    ## [1] 49512    47

##### Liver

``` r
#Subset the count matrix to include only the liver samples
dds_liver <- subset(dds, select=dds$Tissue=="Liver")
dim(dds_liver)
```

    ## [1] 49512    46

``` r
#Create a sample information table comprising only the information of the liver samples
coldata_liver <- coldata[rownames(coldata) %in% colnames(dds_liver) ,]
dim(coldata_liver)
```

    ## [1] 46  4

``` r
#Check the sample order
all(rownames(coldata_liver) %in% colnames(dds_liver))
```

    ## [1] TRUE

``` r
all(rownames(coldata_liver) == colnames(dds_liver))
```

    ## [1] TRUE

``` r
#Drop the unused levels
dds_liver$Age <- droplevels(dds_liver$Age)
dds_liver$Age <- factor(dds_liver$Age, levels=c(3,6,9,12,15,18,21,24,27))
dds_liver$Sex <- droplevels(dds_liver$Sex)
dim(dds_liver)
```

    ## [1] 49512    46

##### Pancreas

``` r
#Subset the count matrix to include only the pancreas samples
dds_pancreas <- subset(dds, select=dds$Tissue=="Pancreas")
dim(dds_pancreas)
```

    ## [1] 49512    34

``` r
#Create a sample information table comprising only the information of the pancreas samples
coldata_pancreas <- coldata[rownames(coldata) %in% colnames(dds_pancreas) ,]
dim(coldata_pancreas)
```

    ## [1] 34  4

``` r
#Check the sample order
all(rownames(coldata_pancreas) %in% colnames(dds_pancreas))
```

    ## [1] TRUE

``` r
all(rownames(coldata_pancreas) == colnames(dds_pancreas))
```

    ## [1] TRUE

``` r
#Drop the unused levels
dds_pancreas$Age <- droplevels(dds_pancreas$Age)
dds_pancreas$Age <- factor(dds_pancreas$Age, levels=c(3,6,9,12,15,18,21,24,27))
dds_pancreas$Sex <- droplevels(dds_pancreas$Sex)
dim(dds_pancreas)
```

    ## [1] 49512    34

### 3.2. Filtering of low expressed genes

##### Brain

``` r
#check minimum age group size
table(dds_brain$Age) #3 samples in the 27-month timepoint 
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  6  6  6  6  6  5  4  3

``` r
#remove genes with less than 3 samples (minimum group size) with counts greater or equal to 10
keep_brain <- rowSums(counts(dds_brain) >= 10 ) >= 3
dds_brain <- dds_brain[keep_brain,]
dim(dds_brain)
```

    ## [1] 34164    48

``` r
#vst transformation
vst_brain <- vst(dds_brain, blind = FALSE)
```

##### Heart

``` r
#check minimum age group size
table(dds_heart$Age) #4 samples in the 24- and 27-month timepoint 
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  5  5  6  6  6  5  6  4  4

``` r
#remove genes with less than 4 samples (minimum group size) with counts greater or equal to 10
keep_heart <- rowSums(counts(dds_heart) >= 10 ) >= 4
dds_heart <- dds_heart[keep_heart,]
dim(dds_heart)
```

    ## [1] 28073    47

``` r
#vst transformation
vst_heart <- vst(dds_heart, blind = FALSE)
```

##### Muscle

``` r
#check minimum age group size
table(dds_muscle$Age) #3 samples in the 24-month timepoint 
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  5  6  6  6  6  5  3  4

``` r
#remove genes with less than 3 samples (minimum group size) with counts greater or equal to 10
keep_muscle <- rowSums(counts(dds_muscle) >= 10 ) >= 3
dds_muscle <- dds_muscle[keep_muscle,]
dim(dds_muscle)
```

    ## [1] 18978    47

``` r
#vst transformation
vst_muscle <- vst(dds_muscle, blind = FALSE)
```

##### Liver

``` r
#check minimum age group size
table(dds_liver$Age) #3 samples in the 24-month timepoint 
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  6  4  6  6  5  6  3  4

``` r
#remove genes with less than 3 samples (minimum group size) with counts greater or equal to 10
keep_liver <- rowSums(counts(dds_liver) >= 10 ) >= 3
dds_liver <- dds_liver[keep_liver,]
dim(dds_liver)
```

    ## [1] 20157    46

``` r
#vst transformation
vst_liver <- vst(dds_liver, blind = FALSE)
```

##### Pancreas

``` r
#check minimum age group size
table(dds_pancreas$Age) #2 samples in the 24-month timepoint 
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  3  4  4  5  5  4  3  2  4

``` r
#remove genes with less than 3 samples (minimum group size) with counts greater or equal to 10
keep_pancreas <- rowSums(counts(dds_pancreas) >= 10 ) >= 2
dds_pancreas <- dds_pancreas[keep_pancreas,]
dim(dds_pancreas)
```

    ## [1] 18414    34

``` r
#vst transformation
vst_pancreas <- vst(dds_pancreas, blind = FALSE)
```

### 3.3. Identification and removal of outlier samples

Next, we used the sample network approach proposed by Oldham et al. (Oldham, Langfelder, and Horvath 2012) to identify outlier samples.

##### Brain

###### Fig S6-Brain

``` r
# sample network based on biweight midcorrelation (bicor)
A = adjacency(assay(vst_brain), type = "distance", corFnc="bicor")
# whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized whole network connectivity
Z.k = scale(k)
# identify samples as outliers if their Z.k values are below the threshold, in this case, more than two standard deviations away from the mean standardized connectivity
thresholdZ.k = -2 


# calculate the cluster tree
sampleTree = hclust(as.dist(1 - A), method = "average")

# define colors for outlier and non-outlier samples
outlierColor = ifelse(Z.k < thresholdZ.k, "darkred", "grey50")
datColors = data.frame(Outlier = outlierColor)

# Plot the sample dendrogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Brain",autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
tiff('Fig-S6-brain.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Brain", autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
dev.off()
```

    ## png 
    ##   2

``` r
#brain outlier samples
outliers_brain <- colnames(vst_brain)[(Z.k < thresholdZ.k)]

dds_brain <- subset(dds_brain, select=!(dds_brain$Sample %in% outliers_brain)) 
dim(dds_brain)
```

    ## [1] 34164    46

``` r
#update sample information
coldata_brain <- coldata_brain[!row.names(coldata_brain) %in% outliers_brain,]
dim(coldata_brain)
```

    ## [1] 46  4

``` r
#Save the sample information table for downstream analyses
write.table(coldata_brain, file="coldata_brain.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#retrieve the normalized and vst-transformed counts for downstream analyses (outside DESeq2, such as WGCNA or segmented regression analysis)
dds_brain <- estimateSizeFactors(dds_brain)
norm_brain_counts <- counts(dds_brain, normalized=TRUE)
write.table(norm_brain_counts, file="norm_brain_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
brain_counts <- counts(dds_brain, normalized=FALSE)
write.table(brain_counts, file="brain_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

vst_brain_counts <- assay(vst(dds_brain, blind = FALSE)) 
write.table(vst_brain_counts, file="vst_brain_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

##### Heart

###### Fig S6-Heart

``` r
# sample network based on biweight midcorrelation (bicor)
A = adjacency(assay(vst_heart), type = "distance", corFnc="bicor")
# whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized whole network connectivity
Z.k = scale(k)
# identify samples as outliers if their Z.k values are below the threshold, in this case, more than two standard deviations away from the mean standardized connectivity
thresholdZ.k = -2 


# calculate the cluster tree
sampleTree = hclust(as.dist(1 - A), method = "average")

# define colors for outlier and non-outlier samples
outlierColor = ifelse(Z.k < thresholdZ.k, "darkred", "grey50")
datColors = data.frame(Outlier = outlierColor)

# Plot the sample dendrogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Heart",autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
tiff('Fig-S6-heart.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Heart", autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
dev.off()
```

    ## png 
    ##   2

``` r
#heart outlier samples
outliers_heart <- colnames(vst_heart)[(Z.k < thresholdZ.k)]

dds_heart <- subset(dds_heart, select=!(dds_heart$Sample %in% outliers_heart)) 
dim(dds_heart)
```

    ## [1] 28073    45

``` r
#update sample information
coldata_heart <- coldata_heart[!row.names(coldata_heart) %in% outliers_heart,]
dim(coldata_heart)
```

    ## [1] 45  4

``` r
#Save the sample information table for downstream analyses
write.table(coldata_heart, file="coldata_heart.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#retrieve the normalized and vst-transformed counts for downstream analyses (outside DESeq2, such as WGCNA or segmented regression analysis)
dds_heart <- estimateSizeFactors(dds_heart)
norm_heart_counts <- counts(dds_heart, normalized=TRUE)
write.table(norm_heart_counts, file="norm_heart_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
heart_counts <- counts(dds_heart, normalized=FALSE)
write.table(heart_counts, file="heart_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

vst_heart_counts <- assay(vst(dds_heart, blind = FALSE)) 
write.table(vst_heart_counts, file="vst_heart_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

##### Muscle

###### Fig S6-Muscle

``` r
# sample network based on biweight midcorrelation (bicor)
A = adjacency(assay(vst_muscle), type = "distance", corFnc="bicor")
# whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized whole network connectivity
Z.k = scale(k)
# identify samples as outliers if their Z.k values are below the threshold, in this case, more than two standard deviations away from the mean standardized connectivity
thresholdZ.k = -2 


# calculate the cluster tree
sampleTree = hclust(as.dist(1 - A), method = "average")

# define colors for outlier and non-outlier samples
outlierColor = ifelse(Z.k < thresholdZ.k, "darkred", "grey50")
datColors = data.frame(Outlier = outlierColor)

# Plot the sample dendrogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Muscle",autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
tiff('Fig-S6-muscle.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Muscle", autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
dev.off()
```

    ## png 
    ##   2

``` r
#muscle outlier samples
outliers_muscle <- colnames(vst_muscle)[(Z.k < thresholdZ.k)]

dds_muscle <- subset(dds_muscle, select=!(dds_muscle$Sample %in% outliers_muscle)) 
dim(dds_muscle)
```

    ## [1] 18978    43

``` r
#update sample information
coldata_muscle <- coldata_muscle[!row.names(coldata_muscle) %in% outliers_muscle,]
dim(coldata_muscle)
```

    ## [1] 43  4

``` r
#Save the sample information table for downstream analyses
write.table(coldata_muscle, file="coldata_muscle.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#retrieve the normalized and vst-transformed counts for downstream analyses (outside DESeq2, such as WGCNA or segmented regression analysis)
dds_muscle <- estimateSizeFactors(dds_muscle)
norm_muscle_counts <- counts(dds_muscle, normalized=TRUE)
write.table(norm_muscle_counts, file="norm_muscle_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
muscle_counts <- counts(dds_muscle, normalized=FALSE)
write.table(muscle_counts, file="muscle_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

vst_muscle_counts <- assay(vst(dds_muscle, blind = FALSE)) 
write.table(vst_muscle_counts, file="vst_muscle_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

##### Liver

###### Fig S6-Liver

``` r
# sample network based on biweight midcorrelation (bicor)
A = adjacency(assay(vst_liver), type = "distance", corFnc="bicor")
# whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized whole network connectivity
Z.k = scale(k)
# identify samples as outliers if their Z.k values are below the threshold, in this case, more than two standard deviations away from the mean standardized connectivity
thresholdZ.k = -2 


# calculate the cluster tree
sampleTree = hclust(as.dist(1 - A), method = "average")

# define colors for outlier and non-outlier samples
outlierColor = ifelse(Z.k < thresholdZ.k, "darkred", "grey50")
datColors = data.frame(Outlier = outlierColor)

# Plot the sample dendrogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Liver",autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
tiff('Fig-S6-liver.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Liver", autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
dev.off()
```

    ## png 
    ##   2

``` r
#liver outlier samples
outliers_liver <- colnames(vst_liver)[(Z.k < thresholdZ.k)]

dds_liver <- subset(dds_liver, select=!(dds_liver$Sample %in% outliers_liver)) 
dim(dds_liver)
```

    ## [1] 20157    45

``` r
#update sample information
coldata_liver <- coldata_liver[!row.names(coldata_liver) %in% outliers_liver,]
dim(coldata_liver)
```

    ## [1] 45  4

``` r
#Save the sample information table for downstream analyses
write.table(coldata_liver, file="coldata_liver.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#retrieve the normalized and vst-transformed counts for downstream analyses (outside DESeq2, such as WGCNA or segmented regression analysis)
dds_liver <- estimateSizeFactors(dds_liver)
norm_liver_counts <- counts(dds_liver, normalized=TRUE)
write.table(norm_liver_counts, file="norm_liver_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
liver_counts <- counts(dds_liver, normalized=FALSE)
write.table(liver_counts, file="liver_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

vst_liver_counts <- assay(vst(dds_liver, blind = FALSE)) 
write.table(vst_liver_counts, file="vst_liver_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

##### Pancreas

###### Fig S6-Pancreas

``` r
# sample network based on biweight midcorrelation (bicor)
A = adjacency(assay(vst_pancreas), type = "distance", corFnc="bicor")
# whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized whole network connectivity
Z.k = scale(k)
# identify samples as outliers if their Z.k values are below the threshold, in this case, more than two standard deviations away from the mean standardized connectivity
thresholdZ.k = -2 


# calculate the cluster tree
sampleTree = hclust(as.dist(1 - A), method = "average")

# define colors for outlier and non-outlier samples
outlierColor = ifelse(Z.k < thresholdZ.k, "darkred", "grey50")
datColors = data.frame(Outlier = outlierColor)

# Plot the sample dendrogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Pancreas",autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
tiff('Fig-S6-pancreas.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,
    main = "Pancreas", autoColorHeight = FALSE, colorHeight = 0.1, cex.main=2.3, cex.lab=1.8, cex.colorLabels = 1.8, cex.dendroLabels = 1.3)
dev.off()
```

    ## png 
    ##   2

``` r
#pancreas outlier samples
outliers_pancreas <- colnames(vst_pancreas)[(Z.k < thresholdZ.k)]

dds_pancreas <- subset(dds_pancreas, select=!(dds_pancreas$Sample %in% outliers_pancreas)) 
dim(dds_pancreas)
```

    ## [1] 18414    32

``` r
#update sample information
coldata_pancreas <- coldata_pancreas[!row.names(coldata_pancreas) %in% outliers_pancreas,]
dim(coldata_pancreas)
```

    ## [1] 32  4

``` r
#Save the sample information table for downstream analyses
write.table(coldata_pancreas, file="coldata_pancreas.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#retrieve the normalized and vst-transformed counts for downstream analyses (outside DESeq2, such as WGCNA or segmented regression analysis)
dds_pancreas <- estimateSizeFactors(dds_pancreas)
norm_pancreas_counts <- counts(dds_pancreas, normalized=TRUE)
write.table(norm_pancreas_counts, file="norm_pancreas_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
pancreas_counts <- counts(dds_pancreas, normalized=FALSE)
write.table(pancreas_counts, file="pancreas_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

vst_pancreas_counts <- assay(vst(dds_pancreas, blind = FALSE)) 
write.table(vst_pancreas_counts, file="vst_pancreas_counts.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

#### 3.4. Principal component analysis of each tissue

##### Fig S1

``` r
#brain
object=vst(dds_brain, blind = FALSE)
intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar1 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d1 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d1, "percentVar") <- percentVar1[1:2]
        return(d1)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar2 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d2 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d2, "percentVar") <- percentVar2[1:2]
        return(d2)
    }

#heart
object=vst(dds_heart, blind = FALSE)
intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar3 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d3 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d3, "percentVar") <- percentVar3[1:2]
        return(d3)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar4 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d4 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d4, "percentVar") <- percentVar4[1:2]
        return(d4)
    }

#liver
object=vst(dds_liver, blind = FALSE)
intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar5 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d5 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d5, "percentVar") <- percentVar5[1:2]
        return(d5)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar6 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d6 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d6, "percentVar") <- percentVar6[1:2]
        return(d6)
    }

#muscle
object=vst(dds_muscle, blind = FALSE)
intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar7 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d7 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d7, "percentVar") <- percentVar7[1:2]
        return(d7)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar8 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d8 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d8, "percentVar") <- percentVar8[1:2]
        return(d8)
    }

#pancreas
object=vst(dds_pancreas, blind = FALSE)
intgroup = "Age"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar9 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d9 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d9, "percentVar") <- percentVar9[1:2]
        return(d9)
    }

intgroup = "Sex"
ntop = 500
returnData = FALSE

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar10 <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

d10 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d10, "percentVar") <- percentVar10[1:2]
        return(d10)
    }


d1$Effect <- "Age"
d2$Effect <- "Sex"
d3$Effect <- "Age"
d4$Effect <- "Sex"
d5$Effect <- "Age"
d6$Effect <- "Sex"
d7$Effect <- "Age"
d8$Effect <- "Sex"
d9$Effect <- "Age"
d10$Effect <- "Sex"

d1$Tissue <- "Brain"
d2$Tissue <- "Brain"
d3$Tissue <- "Heart"
d4$Tissue <- "Heart"
d5$Tissue <- "Liver"
d6$Tissue <- "Liver"
d7$Tissue <- "Muscle"
d8$Tissue <- "Muscle"
d9$Tissue <- "Pancreas"
d10$Tissue <- "Pancreas"

name <- c(d1$name,d2$name,d3$name,d4$name,d5$name,d6$name,d7$name,d8$name,d9$name,d10$name)
PC1 <- c(d1$PC1,d2$PC1,d3$PC1,d4$PC1,d5$PC1,d6$PC1,d7$PC1,d8$PC1,d9$PC1,d10$PC1)
PC2 <- c(d1$PC2,d2$PC2,d3$PC2,d4$PC2,d5$PC2,d6$PC2,d7$PC2,d8$PC2,d9$PC2,d10$PC2)
group <- c(as.character(d1$group),as.character(d2$group),as.character(d3$group),as.character(d4$group),as.character(d5$group),as.character(d6$group),as.character(d7$group),as.character(d8$group),as.character(d9$group),as.character(d10$group))
tissue <- c(d1$Tissue,d2$Tissue,d3$Tissue,d4$Tissue,d5$Tissue,d6$Tissue,d7$Tissue,d8$Tissue,d9$Tissue,d10$Tissue)
effect <- c(d1$Effect,d2$Effect,d3$Effect,d4$Effect,d5$Effect,d6$Effect,d7$Effect,d8$Effect,d9$Effect,d10$Effect)

d = data.frame(name,PC1,PC2,group,tissue,effect)
d$group <- factor(d$group, levels = c("3","6","9","12","15","18","21","24","27","Female","Male"))


mycolors <- c('#e6194B', '#f58231', '#ffe119',
      '#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#F8766D','#619CFF')


a <- ggplot(data = d, aes_string(x = "PC1", y = "PC2",
      color = "group")) + geom_point(size = 3) + xlab("PC1") + ylab("PC2") +
      theme_light() + scale_color_manual(values = mycolors) +
      theme(axis.text=element_text(size=15, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
            axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            legend.title =element_blank(), legend.text = element_text(size=20),
            axis.text.x = element_text(hjust = 0.5, face="bold"),
            strip.text = element_text(size=25, face="bold")) +
      facet_grid(rows=vars(effect), cols=vars(tissue), scales = "free", labeller=labeller(tissue = as_labeller(c(
        Brain = paste0("Brain\nPC1: ", round(percentVar1[1] * 100), "%\nPC2: ",round(percentVar1[2] * 100), "%"),
        Heart = paste0("Heart\nPC1: ", round(percentVar3[1] * 100), "%\nPC2: ",round(percentVar3[2] * 100), "%"),
        Liver = paste0("Liver\nPC1: ", round(percentVar5[1] * 100), "%\nPC2: ",round(percentVar5[2] * 100), "%"),
        Muscle = paste0("Muscle\nPC1: ", round(percentVar7[1] * 100), "%\nPC2: ",round(percentVar7[2] * 100), "%"),
        Pancreas = paste0("Pancreas\nPC1: ", round(percentVar9[1] * 100), "%\nPC2: ",round(percentVar9[2] * 100), "%")) )))


a

g <- ggplot_gtable(ggplot_build(a))
stript <- which(grepl('strip-t', g$layout$name))
fills <- c("#77AADD","#EE8866","#EEDD88","#FFAABB","#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g)
```

![](Normalization_files/figure-markdown_github/unnamed-chunk-32-1.png)

``` r
tiff('Fig-S1.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
grid.draw(g)
dev.off()
```

    ## png 
    ##   2

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 15063)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggpubr_0.4.0                cowplot_1.1.0              
    ##  [3] plyr_1.8.6                  forcats_0.5.0              
    ##  [5] stringr_1.4.0               purrr_0.3.4                
    ##  [7] readr_1.4.0                 tidyr_1.1.2                
    ##  [9] tibble_3.0.4                tidyverse_1.3.0            
    ## [11] WGCNA_1.69                  fastcluster_1.1.25         
    ## [13] dynamicTreeCut_1.63-1       ggplot2_3.3.2              
    ## [15] dplyr_1.0.2                 biomaRt_2.44.4             
    ## [17] DESeq2_1.28.1               SummarizedExperiment_1.18.2
    ## [19] DelayedArray_0.14.1         matrixStats_0.57.0         
    ## [21] Biobase_2.50.0              GenomicRanges_1.40.0       
    ## [23] GenomeInfoDb_1.24.2         IRanges_2.24.0             
    ## [25] S4Vectors_0.28.0            BiocGenerics_0.36.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1           backports_1.2.0        Hmisc_4.4-1           
    ##   [4] BiocFileCache_1.12.1   splines_4.0.3          BiocParallel_1.24.1   
    ##   [7] digest_0.6.27          foreach_1.5.1          htmltools_0.5.0       
    ##  [10] GO.db_3.12.1           magrittr_2.0.1         checkmate_2.0.0       
    ##  [13] memoise_1.1.0          cluster_2.1.0          doParallel_1.0.16     
    ##  [16] openxlsx_4.2.3         annotate_1.66.0        modelr_0.1.8          
    ##  [19] askpass_1.1            prettyunits_1.1.1      jpeg_0.1-8.1          
    ##  [22] colorspace_2.0-0       blob_1.2.1             rvest_0.3.6           
    ##  [25] rappdirs_0.3.1         haven_2.3.1            xfun_0.19             
    ##  [28] crayon_1.3.4           RCurl_1.98-1.2         jsonlite_1.7.1        
    ##  [31] genefilter_1.70.0      impute_1.62.0          survival_3.2-7        
    ##  [34] iterators_1.0.13       glue_1.4.2             gtable_0.3.0          
    ##  [37] zlibbioc_1.34.0        XVector_0.28.0         car_3.0-10            
    ##  [40] abind_1.4-5            scales_1.1.1           DBI_1.1.0             
    ##  [43] rstatix_0.6.0          Rcpp_1.0.5             xtable_1.8-4          
    ##  [46] progress_1.2.2         htmlTable_2.1.0        foreign_0.8-80        
    ##  [49] bit_4.0.4              preprocessCore_1.50.0  Formula_1.2-4         
    ##  [52] htmlwidgets_1.5.2      httr_1.4.2             RColorBrewer_1.1-2    
    ##  [55] ellipsis_0.3.1         farver_2.0.3           pkgconfig_2.0.3       
    ##  [58] XML_3.99-0.5           nnet_7.3-14            dbplyr_1.4.4          
    ##  [61] locfit_1.5-9.4         tidyselect_1.1.0       labeling_0.4.2        
    ##  [64] rlang_0.4.9            AnnotationDbi_1.52.0   munsell_0.5.0         
    ##  [67] cellranger_1.1.0       tools_4.0.3            cli_3.0.0             
    ##  [70] generics_0.1.0         RSQLite_2.2.1          broom_0.7.2           
    ##  [73] evaluate_0.14          yaml_2.2.1             knitr_1.30            
    ##  [76] bit64_4.0.5            fs_1.5.0               zip_2.1.1             
    ##  [79] xml2_1.3.2             compiler_4.0.3         rstudioapi_0.13       
    ##  [82] curl_4.3               png_0.1-7              ggsignif_0.6.0        
    ##  [85] reprex_2.0.0           geneplotter_1.66.0     stringi_1.5.3         
    ##  [88] lattice_0.20-41        Matrix_1.2-18          vctrs_0.3.5           
    ##  [91] pillar_1.4.7           lifecycle_0.2.0        data.table_1.13.4     
    ##  [94] bitops_1.0-6           R6_2.5.0               latticeExtra_0.6-29   
    ##  [97] gridExtra_2.3          rio_0.5.16             codetools_0.2-16      
    ## [100] assertthat_0.2.1       openssl_1.4.3          withr_2.3.0           
    ## [103] GenomeInfoDbData_1.2.3 hms_0.5.3              rpart_4.1-15          
    ## [106] rmarkdown_2.5          carData_3.0-4          lubridate_1.7.9       
    ## [109] base64enc_0.1-3

Durinck, Steffen, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma, and Wolfgang Huber. 2005. “BioMart and Bioconductor: A powerful link between biological databases and microarray data analysis.” *Bioinformatics* 21 (16): 3439–40. doi:[10.1093/bioinformatics/bti525](https://doi.org/10.1093/bioinformatics/bti525).

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” *Genome Biology* 15 (12): 550. doi:[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).

Oldham, Michael C., Peter Langfelder, and Steve Horvath. 2012. “Network methods for describing sample relationships in genomic datasets: application to Huntington’s disease.” *BMC Systems Biology* 6 (1). BioMed Central: 63. doi:[10.1186/1752-0509-6-63](https://doi.org/10.1186/1752-0509-6-63).

Schaum, Nicholas, Benoit Lehallier, Oliver Hahn, Róbert Pálovics, Shayan Hosseinzadeh, Song E Lee, Rene Sit, et al. 2020. “Ageing hallmarks exhibit organ-specific temporal signatures.” *Nature* 583: 596–602. doi:[10.1038/s41586-020-2499-y](https://doi.org/10.1038/s41586-020-2499-y).
