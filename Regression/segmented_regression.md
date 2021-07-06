Segmented regression analysis
================
Margarida Ferreira
Last updated July 2021

-   [1. Initial setup](#initial-setup)
-   [2. Data input](#data-input)
    -   [2.1 Normalized data](#normalized-data)
    -   [2.2 Time vector and gender co-variate](#time-vector-and-gender-co-variate)
-   [3. Regression analysis](#regression-analysis)
    -   [3.1 Altered functions](#altered-functions)
    -   [3.2 Run method](#run-method)
    -   [3.3 Choice of threshold for dynamic Trendy fit](#choice-of-threshold-for-dynamic-trendy-fit)
    -   [3.4 Top Dynamic Genes](#top-dynamic-genes)
        -   [3.4.1 Summary tables](#summary-tables)
        -   [Fig 1B](#fig-1b)
    -   [3.5 Test differences in median breakpoint distributions](#test-differences-in-median-breakpoint-distributions)
        -   [Fig 1C](#fig-1c)
        -   [Fig 1D](#fig-1d)
    -   [3.6 Test biotype enrichment](#test-biotype-enrichment)
        -   [Fig 1E](#fig-1e)
-   [Session info](#session-info)

The regression analysis was carried out using the Bioconductor's package Trendy (Bacher et al. 2018).

# 1. Initial setup

Loading the required packages:

``` r
library(Trendy)
library(BiocParallel)
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(cowplot)
library(ggpubr)
library(plyr)
library(FSA)
```

# 2. Data input

## 2.1 Normalized data

We use DESeq2's median of ratios normalization (performed in the `Normalization.md`).

``` r
brain_counts <- as.matrix(read.table("norm_brain_counts.txt", header=TRUE))

heart_counts <- as.matrix(read.table("norm_heart_counts.txt", header=TRUE))

muscle_counts <- as.matrix(read.table("norm_muscle_counts.txt", header=TRUE))

liver_counts <- as.matrix(read.table("norm_liver_counts.txt", header=TRUE))

pancreas_counts <- as.matrix(read.table("norm_pancreas_counts.txt", header=TRUE))

dim(brain_counts)
```

    ## [1] 34164    46

``` r
dim(heart_counts)
```

    ## [1] 28073    45

``` r
dim(muscle_counts)
```

    ## [1] 18978    43

``` r
dim(liver_counts)
```

    ## [1] 20157    45

``` r
dim(pancreas_counts)
```

    ## [1] 18414    32

It is important to note that the [*"samples should be sorted in the time course order."*](http://www.bioconductor.org/packages/release/bioc/vignettes/Trendy/inst/doc/Trendy_vignette.pdf)

``` r
###BRAIN
#read in the sample information table
coldata_brain <- read.table("coldata_brain.txt", header=TRUE)
#order by increasing age
coldata_brain <- coldata_brain[order(coldata_brain$Age, decreasing=F),]
#re-order samples in counts table based on increasing age
brain_counts <- brain_counts[,match(coldata_brain$Sample, colnames(brain_counts))]

###HEART
coldata_heart <- read.table("coldata_heart.txt", header=TRUE)
coldata_heart <- coldata_heart[order(coldata_heart$Age, decreasing=F),]
heart_counts <- heart_counts[,match(coldata_heart$Sample, colnames(heart_counts))]

###MUSCLE
coldata_muscle <- read.table("coldata_muscle.txt", header=TRUE)
coldata_muscle <- coldata_muscle[order(coldata_muscle$Age, decreasing=F),]
muscle_counts <- muscle_counts[,match(coldata_muscle$Sample, colnames(muscle_counts))]

###LIVER
coldata_liver <- read.table("coldata_liver.txt", header=TRUE)
coldata_liver <- coldata_liver[order(coldata_liver$Age, decreasing=F),]
liver_counts <- liver_counts[,match(coldata_liver$Sample, colnames(liver_counts))]

###PANCREAS
coldata_pancreas <- read.table("coldata_pancreas.txt", header=TRUE)
coldata_pancreas <- coldata_pancreas[order(coldata_pancreas$Age, decreasing=F),]
pancreas_counts <- pancreas_counts[,match(coldata_pancreas$Sample, colnames(pancreas_counts))]
```

## 2.2 Time vector and gender co-variate

``` r
#BRAIN
table(coldata_brain$Age)
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  6  5  6  5  6  5  4  3

``` r
coldata_brain$Gender[coldata_brain$Sex == "Female"] <- "0"
coldata_brain$Gender[coldata_brain$Sex == "Male"] <- "1"

#replace column names in the counts matrix for simplicity
colnames(brain_counts) <- c(paste0("3m_s", 1:6),paste0("6m_s", 1:6),paste0("9m_s", 1:5),paste0("12m_s", 1:6),paste0("15m_s", 1:5),paste0("18m_s", 1:6),paste0("21m_s", 1:5),paste0("24m_s", 1:4),paste0("27m_s", 1:3))

#HEART
table(coldata_heart$Age)
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  5  5  6  5  6  4  6  4  4

``` r
coldata_heart$Gender[coldata_heart$Sex == "Female"] <- "0"
coldata_heart$Gender[coldata_heart$Sex == "Male"] <- "1"

#replace column names in the counts matrix for simplicity
colnames(heart_counts) <- c(paste0("3m_s", 1:5),paste0("6m_s", 1:5),paste0("9m_s", 1:6),paste0("12m_s", 1:5),paste0("15m_s", 1:6),paste0("18m_s", 1:4),paste0("21m_s", 1:6),paste0("24m_s", 1:4),paste0("27m_s", 1:4))

#MUSCLE
table(coldata_muscle$Age)
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  4  6  5  6  6  5  2  3

``` r
coldata_muscle$Gender[coldata_muscle$Sex == "Female"] <- "0"
coldata_muscle$Gender[coldata_muscle$Sex == "Male"] <- "1"

#replace column names in the counts matrix for simplicity
colnames(muscle_counts) <- c(paste0("3m_s", 1:6),paste0("6m_s", 1:4),paste0("9m_s", 1:6),paste0("12m_s", 1:5),paste0("15m_s", 1:6),paste0("18m_s", 1:6),paste0("21m_s", 1:5),paste0("24m_s", 1:2),paste0("27m_s", 1:3))

#LIVER
table(coldata_liver$Age)
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  6  6  4  6  6  5  5  3  4

``` r
coldata_liver$Gender[coldata_liver$Sex == "Female"] <- "0"
coldata_liver$Gender[coldata_liver$Sex == "Male"] <- "1"

#replace column names in the counts matrix for simplicity
colnames(liver_counts) <- c(paste0("3m_s", 1:6),paste0("6m_s", 1:6),paste0("9m_s", 1:4),paste0("12m_s", 1:6),paste0("15m_s", 1:6),paste0("18m_s", 1:5),paste0("21m_s", 1:5),paste0("24m_s", 1:3),paste0("27m_s", 1:4))

#PANCREAS
table(coldata_pancreas$Age)
```

    ## 
    ##  3  6  9 12 15 18 21 24 27 
    ##  3  4  4  5  4  4  3  2  3

``` r
coldata_pancreas$Gender[coldata_pancreas$Sex == "Female"] <- "0"
coldata_pancreas$Gender[coldata_pancreas$Sex == "Male"] <- "1"

#replace column names in the counts matrix for simplicity
colnames(pancreas_counts) <- c(paste0("3m_s", 1:3),paste0("6m_s", 1:4),paste0("9m_s", 1:4),paste0("12m_s", 1:5),paste0("15m_s", 1:4),paste0("18m_s", 1:4),paste0("21m_s", 1:3),paste0("24m_s", 1:2),paste0("27m_s", 1:3))
```

# 3. Regression analysis

## 3.1 Altered functions

To accomodate adjustment for co-variates, some functions were altered to include a 'covar' parameter.

``` r
fitSegBIC <- function(Data, maxK = 2, tVectIn = NULL, covar=NULL,
                      minNumInSeg = 5, pvalCut = .1,
                      numTry = 5, keepFit = FALSE)

{

    whichFit <- seq_len(maxK)

        # If any replicates, jitter a small amount to help with the segmented fitting:
        if (length(unique(tVectIn)) < length(tVectIn)) {tVectIn <- jitter(tVectIn, .1)}

    # Start with lm without any breaks
    lmLinear <- lm(Data ~ tVectIn+covar)
    lm.radj <- summary(lmLinear)$adj.r.squared
    lm.rsq <- summary(lmLinear)$r.squared
    lm.slp <- coef(lmLinear)[2]
    names(lm.slp) <- paste0("Segment", seq_len(length(lm.slp)), ".Slope")
    lm.fit <- fitted.values(lmLinear)
    lm.pval <- coef(summary(lmLinear))[2,4]
    lm.sign <- ifelse(lm.slp > 0, 1, -1)

    lm.sign[which(lm.pval > pvalCut)] <- 0
    names(lm.sign) <- paste0("Segment", seq_len(length(lm.sign)), ".Trend")

    lm.id.sign <- rep(lm.sign, length(tVectIn))
    names(lm.id.sign) <- paste0(names(tVectIn), ".Trend")
    names(lm.fit) <- paste0(names(tVectIn), ".Fitted")
    bic.lm <- BIC(lmLinear)

    # Fit all possible breakpoints now:
    fit.bp.all <- lapply(whichFit, breakpointFit, tVectIn, lmLinear, numTry)

    isna <- which(vapply(fit.bp.all, function(i) {
                (class(i)[1] == "character")
              }, logical(1)))

    if (length(isna) == maxK) { #No bp models valid, return linear results
        OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
            Segment.Trends = lm.sign,
            Segment.Pvalues = lm.pval, Breakpoints = NA,
            Fitted.Values = lm.fit,
            AdjustedR2 = lm.radj, Fit = lmLinear)
            if (keepFit == FALSE) { OUT <- OUT[seq_len(7)] }
                return(OUT)
            # If it is not solved in 100 trys then return lm
            # results (if numTry=100)
    }

    fit.keep <- fit.bp.all
    if (length(isna) > 0) { # if one of whichFit cant be fitted then remove it.
        fit.keep <- fit.keep[-isna]
        whichFit <- whichFit[-isna]
    }

    # Get info for each fit:
    slp.l <- lapply(fit.keep, function(i) {segmented::slope(i)})
    radj <- vapply(fit.keep, function(i) {summary(i)$adj.r.squared},numeric(1))
    brk.l <- lapply(fit.keep, function(i) {i$psi[,2]})
    id.l <- lapply(fit.keep, function(i) {i$id.group})
    allBIC <- vapply(fit.keep, function(i) {BIC(i)}, numeric(1))

    if (length(whichFit) >= 1) {
        bic.whichmin <- which.min(allBIC)
        if (length(bic.whichmin) > 0) {
            r.choose <- max(bic.whichmin)
        }
        # compare here using bic.
        if (allBIC[r.choose] >= bic.lm) {
            # if linear has smaller BIC then take the linear, return values
            OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
                Segment.Trends = lm.sign,
                Segment.Pvalues = lm.pval, Breakpoints = NA,
                Fitted.Values = lm.fit,
                AdjustedR2 = lm.radj, Fit = lmLinear)
                if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                    return(OUT)
            }

            # now make sure that best fit model satisfies having minimum
            # number of segments, if it does not then decrease breakpoints.
            while((min(table(id.l[[r.choose]])) < minNumInSeg | any(brk.l[[r.choose]] < tVectIn[minNumInSeg])) & r.choose > 1) {
                r.choose <- r.choose - 1
            }
                        checkMinSeg = 0
                        while (any(checkMinSeg < minNumInSeg)) {
                            checkMinSeg <- c()
              breaks.temp <- brk.l[[r.choose]]
              if (length(breaks.temp) == 1) {
                checkMinSeg <- c(length(tVectIn[tVectIn < breaks.temp]), length(tVectIn[tVectIn > breaks.temp]))
              }
              if (length(breaks.temp) >= 2) {
                breaks.temp <- c(tVectIn[1], breaks.temp, tVectIn[length(tVectIn)])
                for (i in 1:(length(breaks.temp)-1)) {
                                checkMinSeg <- c(checkMinSeg, length(tVectIn[tVectIn < breaks.temp[i+1] & tVectIn > breaks.temp[i]]))
                              }
              }
                            if (any(checkMinSeg < minNumInSeg)) {r.choose <- r.choose - 1}
                            if (r.choose == 0) {break}
                        }
                if (r.choose == 0) {
                    # take the linear, return values
                    OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
                        Segment.Trends = lm.sign,
                        Segment.Pvalues = lm.pval, Breakpoints = NA,
                        Fitted.Values = lm.fit,
                        AdjustedR2 = lm.radj, Fit = lmLinear)
                        if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                            return(OUT)
                    }
            if (r.choose == 1 & (min(table(id.l[[r.choose]])) < minNumInSeg | any(brk.l[[r.choose]] < tVectIn[minNumInSeg]))){
                # if 1 bp gives too small segment, then take the linear
                OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
                    Segment.Trends = lm.sign,
                    Segment.Pvalues = lm.pval, Breakpoints = NA,
                    Fitted.Values = lm.fit,
                    AdjustedR2 = lm.radj, Fit = lmLinear)
                    if(keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
                        return(OUT)
            }
    }

    # Finally decide if best BP model is better than linear,
    # If not return linear
    if (allBIC[r.choose] >= bic.lm) {
        OUT <- list(Trends = lm.id.sign, Segment.Slopes = lm.slp,
            Segment.Trends = lm.sign,
            Segment.Pvalues = lm.pval, Breakpoints = NA,
            Fitted.Values = lm.fit,
            AdjustedR2 = lm.radj, Fit = lmLinear)
        if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
        return(OUT)
        # If lm is better
    }

    # # Actual k (before remove the NA fitting)
    r.choose.ori <- whichFit[r.choose]

    fit.choose <- fit.keep[[r.choose]]
    fv.choose <- fitted.values(fit.choose)
    names(fv.choose) <- paste0(names(tVectIn), ".Fitted")

    bp.choose <- brk.l[[r.choose]]
    if (length(bp.choose) >= 1) {
        names(bp.choose) <- paste0("Breakpoint", seq_along(bp.choose))}
    slp.choose <- slp.l[[r.choose]][[1]][,1]
    names(slp.choose) <- paste0("Segment",seq_len(length(slp.choose)),".Slope")
  

    slp.t <- slp.l[[r.choose]][[1]][,3]
    slp.pval <- pt(-abs(slp.t), 1)
    names(slp.pval) <- paste0("Segment", seq_len(length(slp.pval)), ".Pvalue")

    slp.sign <- ifelse(slp.t > 0, 1, -1)
    slp.sign[which(slp.pval > pvalCut)] <- 0
    names(slp.sign) <- paste0("Segment", seq_len(length(slp.sign)), ".Trend")

    id.choose <- id.l[[r.choose]]
    id.sign <- slp.sign[id.choose + 1]
    names(id.sign) <- paste0(names(tVectIn), ".Trend")

    OUT = list(Trends = id.sign, Segment.Slopes = slp.choose,
        Segment.Trends = slp.sign,
        Segment.Pvalues = slp.pval, Breakpoints = bp.choose,
        Fitted.Values = fv.choose,
        AdjustedR2 = radj[r.choose], Fit = fit.choose)
    if (keepFit == FALSE) {OUT <- OUT[seq_len(7)]}
  
    return(OUT)
}
```

``` r
trendy <-
    function(Data = NULL, tVectIn = NULL, covar= NULL, saveObject = FALSE, fileName = NULL,
            meanCut = 10, maxK = 3, minNumInSeg = 5, pvalCut = .1,
            numTry = 5, keepFit = FALSE, NCores = NULL, featureNames = NULL)
{
    # Checks
    if (methods::is(Data, "SummarizedExperiment")) {
        if (is.null(SummarizedExperiment::assayNames(Data)) ||
        SummarizedExperiment::assayNames(Data)[1] != "Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"
            if (is.null(colnames(Trendy::getCounts(Data)))) {
                stop("Must supply sample/cell names!")
            }
        }
    }
    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
        Data<-SummarizedExperiment::SummarizedExperiment(assays=
            list("Counts"=Data))
    }
    if (anyNA(Trendy::getCounts(Data))) {stop("Data contains at least one
        value of NA. Unsure how to proceed.")}
    if (is.null(rownames(Trendy::getCounts(Data)))) {
        stop("Must supply feature/gene/row names!")
    }
    if (is.null(colnames(Trendy::getCounts(Data)))) {
        stop("Must supply sample/column names!")
    }
    NSample <- ncol(Trendy::getCounts(Data))
    if (is.null(tVectIn)) {
        warning(paste0("No values for parameter tVectIn were given.
        Trendy will assume data goes from 1:",NSample))
        tVectIn <- seq_len(NSample)
        names(tVectIn) <- colnames(Trendy::getCounts(Data))
    }
    if (is.null(names(tVectIn))) {
        names(tVectIn) <- colnames(Trendy::getCounts(Data))
    }
    if (is.null(NCores)) {NCores <- max(1, parallel::detectCores() - 1)}
    if (.Platform$OS.type == "windows") {
        param = SnowParam(workers=NCores)
    }

    param = BiocParallel::MulticoreParam(workers=NCores)

    BiocParallel::register(BPPARAM = param)

    if (is.null(featureNames)) {
        featureNames <- rownames(Data)
    }

    Data <- Data[rownames(Data) %in% featureNames, ]

    if (length(featureNames) == 1) {

        if (mean(Trendy::getCounts(Data)) >= meanCut) {
            Data.MeanFiltered <- t(data.matrix(Trendy::getCounts(Data)))
            row.names(Data.MeanFiltered) <- featureNames
        } else {
            stop("Gene does not pass the mean cutoff filter!")
        }
    } else {
        toKeep <- which(rowMeans(Trendy::getCounts(Data)) >= meanCut)
        Data.MeanFiltered<-
        Trendy::getCounts(Data)[toKeep,]
        if (sum(rowMeans(Trendy::getCounts(Data)) >= meanCut) == 0)  {
        stop("No genes pass the mean cutoff filter!")
    }
    }

    if (NSample < (maxK + 1) * minNumInSeg) {
        maxK <- floor(NSample / minNumInSeg) - 1
        message("Number of samples (", NSample, ") is less than
        [# segments] * [min number of samples in a segment]. maxK has been
        set to", maxK)
    }
    if (length(unique(tVectIn)) <= (maxK + 1)) {
        maxK <- length(unique(tVectIn)) - 2
        message("Number of unique times (", length(unique(tVectIn)), ") is less than
        setting of maxK. Trendy has automatically set maxK to ", maxK)
    }
    if (maxK < 1) {
        stop("Invalid value for maxK. Adjust minNumInSeg setting
        in order to run Trendy.")
    }

    segAll <- BiocParallel::bplapply(X = seq_len(nrow(Data.MeanFiltered)),
      function(X) {
        inGene = Data.MeanFiltered[X,]
        fitSegBIC(Data = inGene,
            tVectIn = tVectIn,
            covar = covar,
            maxK = maxK,
            minNumInSeg = minNumInSeg,
            pvalCut = pvalCut,
            numTry = numTry,
            keepFit = keepFit)
        })
    names(segAll) <- rownames(Data.MeanFiltered)

    if (saveObject == TRUE) {
        if (is.null(fileName)){
            fileName <- "trendyForShiny.RData"
        } else {
            fileName <- paste0(fileName, "_trendyForShiny.RData")
        }
        origData <- Trendy::getCounts(Data)
        trendyOut <- segAll
        tVectIn <- tVectIn

        save(trendyOut, origData, tVectIn, file = fileName)}

        S4Vectors::metadata(Data)[["TrendyFits"]] <- segAll
        return(Data)
}
```

## 3.2 Run method

``` r
# #BRAIN
# res_brain_cov <- trendy(Data = brain_counts, tVectIn = coldata_brain$Age, covar=coldata_brain$Gender, minNumInSeg=3, pvalCut = 0.1, maxK = 7, meanCut = 0)
# res_brain_cov <- results(res_brain_cov)
# 
# #HEART
# res_heart_cov <- trendy(Data = heart_counts, tVectIn = coldata_heart$Age, covar=coldata_heart$Gender, minNumInSeg=4, pvalCut = 0.1, maxK = 7, meanCut = 0)
# res_heart_cov <- results(res_heart_cov)
# 
# #MUSCLE
# res_muscle_cov <- trendy(Data = muscle_counts, tVectIn = coldata_muscle$Age, covar = coldata_muscle$Gender, minNumInSeg=2, pvalCut = 0.1, maxK = 7, meanCut = 0)
# res_muscle_cov <- results(res_muscle_cov)
# 
# #LIVER
# res_liver_cov <- trendy(Data = liver_counts, tVectIn = coldata_liver$Age, covar = coldata_liver$Gender, minNumInSeg=3, pvalCut = 0.1, maxK = 7, meanCut = 0)
# res_liver_cov <- results(res_liver_cov)
# 
# #PANCREAS
# res_pancreas_cov <- trendy(Data = pancreas_counts, tVectIn = coldata_pancreas$Age, covar = coldata_pancreas$Gender, minNumInSeg=2, pvalCut = 0.1, maxK = 7, meanCut = 0)
# res_pancreas_cov <- results(res_pancreas_cov)

#save R objects to load next time
# saveRDS(res_brain_cov, file = "res_brain_cov.rds")
# 
# saveRDS(res_heart_cov, file = "res_heart_cov.rds")
# 
# saveRDS(res_muscle_cov, file = "res_muscle_cov.rds")
# 
# saveRDS(res_liver_cov, file = "res_liver_cov.rds")
# 
# saveRDS(res_pancreas_cov, file = "res_pancreas_cov.rds")

#restore the objects
res_brain_cov <- readRDS(file = "res_brain_cov.rds")

res_heart_cov <- readRDS(file = "res_heart_cov.rds")

res_muscle_cov <- readRDS(file = "res_muscle_cov.rds")

res_liver_cov <- readRDS(file = "res_liver_cov.rds")

res_pancreas_cov <- readRDS(file = "res_pancreas_cov.rds")
```

## 3.3 Choice of threshold for dynamic Trendy fit

As suggested in the section 4.4 of the Trendy [vignette](https://rdrr.io/bioc/Trendy/f/inst/doc/Trendy_vignette.pdf) and as performed in [Barry et al, 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007543#sec017).

Because the goodness of fit parameter (adjusted R2) is experiment-dependent, to determine the appropriate cutoffs for each tissue we performed a permutation procedure on each dataset. We shuffled the order of the samples 100 times and for each shuffle, and for 100 randomly sampled genes, we ran Trendy using the same parameters as on the original datasets. We then evaluated the adjusted R2 distribution of the Trendy fits on the permuted data.

``` r
#BRAIN
# res_brain.r2 <- c()
# 
# #shuffle time vector and get corresponding gender info
# df_brain<-coldata_brain[,c(3,5)]
# 
# sample_brain <- sample(nrow(df_brain))
# 
# df_sample_brain <- df_brain[sample_brain,]
# 
# 
# for(i in 1:100) { # permute 100 times at least
# BiocParallel::register(BiocParallel::SerialParam())
# seg.shuffle_brain <- trendy(brain_counts[sample(1:nrow(brain_counts), 100),], #sample genes each time
# tVectIn = df_sample_brain$Age, # shuffled time vector
# covar = df_sample_brain$Gender, #shuffled gender vector
# minNumInSeg=3, pvalCut = 0.1, maxK = 8, meanCut = 0,
# saveObject=FALSE)
# res_brain_cov.r2 <- results(seg.shuffle_brain)
# res_brain.r2 <- c(res_brain.r2, sapply(res_brain_cov.r2, function(x) x$AdjustedR2))
# }
# 
# saveRDS(res_brain.r2, file = "res_brain.r2.rds")
res_brain.r2 <- readRDS(file = "res_brain.r2.rds")

# Say you want to use the value such that less than 1% of permutations reach:
round(sort(res_brain.r2, decreasing=T)[round(.01 * length(res_brain.r2))],1)
```

    ## Gm20257 
    ##     0.2

``` r
#We chose an adjusted R2 of 0.2, as less than 1% of the genes in the brain permutations achieved this threshold.

# Histogram of all R^2
pdf("res_brain.r2.pdf")
hist(res_brain.r2, ylim=c(0,1000), xlim=c(0,1), xlab=expression(paste("Adjusted R"^"2")),main="Brain")
abline(v=round(sort(res_brain.r2, decreasing=T)[round(.01 * length(res_brain.r2))],1), lty=2)
dev.off()
```

    ## png 
    ##   2

``` r
#HEART
# res_heart.r2 <- c()
# 
# #shuffle time vector and get corresponding gender info
# df_heart<-coldata_heart[,c(3,5)]
# 
# sample_heart <- sample(nrow(df_heart))
# 
# df_sample_heart <- df_heart[sample_heart,]
# 
# 
# for(i in 1:100) { # permute 100 times at least
# BiocParallel::register(BiocParallel::SerialParam())
# seg.shuffle_heart <- trendy(heart_counts[sample(1:nrow(heart_counts), 100),], #sample genes each time
# tVectIn = df_sample_heart$Age, # shuffled time vector
# covar = df_sample_heart$Gender, #shuffled gender vector
# minNumInSeg=4, pvalCut = 0.1, maxK = 8, meanCut = 0,
# saveObject=FALSE)
# res_heart_cov.r2 <- results(seg.shuffle_heart)
# res_heart.r2 <- c(res_heart.r2, sapply(res_heart_cov.r2, function(x) x$AdjustedR2))
# }
# 
# saveRDS(res_heart.r2, file = "res_heart.r2.rds")
res_heart.r2 <- readRDS(file = "res_heart.r2.rds")

# Say you want to use the value such that less than 1% of permutations reach:
round(sort(res_heart.r2, decreasing=T)[round(.01 * length(res_heart.r2))],1)
```

    ## Spag11b 
    ##     0.2

``` r
#We chose an adjusted R2 of 0.2, as less than 1% of the genes in the heart permutations achieved this threshold.

# Histogram of all R^2
pdf("res_heart.r2.pdf")
hist(res_heart.r2, ylim=c(0,1000), xlim=c(0,1), xlab=expression(paste("Adjusted R"^"2")),main="Heart")
abline(v=round(sort(res_heart.r2, decreasing=T)[round(.01 * length(res_heart.r2))],1), lty=2)
dev.off()
```

    ## png 
    ##   2

``` r
#MUSCLE
# res_muscle.r2 <- c()
# 
# #shuffle time vector and get corresponding gender info
# df_muscle<-coldata_muscle[,c(3,5)]
# 
# sample_muscle <- sample(nrow(df_muscle))
# 
# df_sample_muscle <- df_muscle[sample_muscle,]
# 
# 
# for(i in 1:100) { # permute 100 times at least
# BiocParallel::register(BiocParallel::SerialParam())
# seg.shuffle_muscle <- trendy(muscle_counts[sample(1:nrow(muscle_counts), 100),], #sample genes each time
# tVectIn = df_sample_muscle$Age, # shuffled time vector
# covar = df_sample_muscle$Gender, #shuffled gender vector
# minNumInSeg=2, pvalCut = 0.1, maxK = 8, meanCut = 0,
# saveObject=FALSE)
# res_muscle_cov.r2 <- results(seg.shuffle_muscle)
# res_muscle.r2 <- c(res_muscle.r2, sapply(res_muscle_cov.r2, function(x) x$AdjustedR2))
# }
# 
# saveRDS(res_muscle.r2, file = "res_muscle.r2.rds")
res_muscle.r2 <- readRDS(file = "res_muscle.r2.rds")

# Say you want to use the value such that less than 1% of permutations reach:
round(sort(res_muscle.r2, decreasing=T)[round(.01 * length(res_muscle.r2))],1)
```

    ## A430078I02Rik 
    ##           0.3

``` r
#We chose an adjusted R2 of 0.3, as less than 1% of the genes in the muscle permutations achieved this threshold.

# Histogram of all R^2
pdf("res_muscle.r2.pdf")
hist(res_muscle.r2, ylim=c(0,1000), xlim=c(0,1), xlab=expression(paste("Adjusted R"^"2")),main="Muscle")
abline(v=round(sort(res_muscle.r2, decreasing=T)[round(.01 * length(res_muscle.r2))],1), lty=2)
dev.off()
```

    ## png 
    ##   2

``` r
#LIVER
# res_liver.r2 <- c()
# 
# #shuffle time vector and get corresponding gender info
# df_liver<-coldata_liver[,c(3,5)]
# 
# sample_liver <- sample(nrow(df_liver))
# 
# df_sample_liver <- df_liver[sample_liver,]
# 
# 
# for(i in 1:100) { # permute 100 times at least
# BiocParallel::register(BiocParallel::SerialParam())
# seg.shuffle_liver <- trendy(liver_counts[sample(1:nrow(liver_counts), 100),], #sample genes each time
# tVectIn = df_sample_liver$Age, # shuffled time vector
# covar = df_sample_liver$Gender, #shuffled gender vector
# minNumInSeg=3, pvalCut = 0.1, maxK = 8, meanCut = 0,
# saveObject=FALSE)
# res_liver_cov.r2 <- results(seg.shuffle_liver)
# res_liver.r2 <- c(res_liver.r2, sapply(res_liver_cov.r2, function(x) x$AdjustedR2))
# }
# 
# saveRDS(res_liver.r2, file = "res_liver.r2.rds")
res_liver.r2 <- readRDS(file = "res_liver.r2.rds")

# Say you want to use the value such that less than 1% of permutations reach:
round(sort(res_liver.r2, decreasing=T)[round(.01 * length(res_liver.r2))],1)
```

    ## Rapgef5 
    ##     0.1

``` r
#We chose an adjusted R2 of 0.1, as less than 1% of the genes in the liver permutations achieved this threshold.

# Histogram of all R^2
pdf("res_liver.r2.pdf")
hist(res_liver.r2, ylim=c(0,1000), xlim=c(0,1), xlab=expression(paste("Adjusted R"^"2")),main="Liver")
abline(v=round(sort(res_liver.r2, decreasing=T)[round(.01 * length(res_liver.r2))],1), lty=2)
dev.off()
```

    ## png 
    ##   2

``` r
#PANCREAS
# res_pancreas.r2 <- c()
# 
# #shuffle time vector and get corresponding gender info
# df_pancreas<-coldata_pancreas[,c(3,5)]
# 
# sample_pancreas <- sample(nrow(df_pancreas))
# 
# df_sample_pancreas <- df_pancreas[sample_pancreas,]
# 
# 
# for(i in 1:100) { # permute 100 times at least
# BiocParallel::register(BiocParallel::SerialParam())
# seg.shuffle_pancreas <- trendy(pancreas_counts[sample(1:nrow(pancreas_counts), 100),], #sample genes each time
# tVectIn = df_sample_pancreas$Age, # shuffled time vector
# covar = df_sample_pancreas$Gender, #shuffled gender vector
# minNumInSeg=2, pvalCut = 0.1, maxK = 8, meanCut = 0,
# saveObject=FALSE)
# res_pancreas_cov.r2 <- results(seg.shuffle_pancreas)
# res_pancreas.r2 <- c(res_pancreas.r2, sapply(res_pancreas_cov.r2, function(x) x$AdjustedR2))
# }
# 
# saveRDS(res_pancreas.r2, file = "res_pancreas.r2.rds")
res_pancreas.r2 <- readRDS(file = "res_pancreas.r2.rds")

# Say you want to use the value such that less than 1% of permutations reach:
round(sort(res_pancreas.r2, decreasing=T)[round(.01 * length(res_pancreas.r2))],1)
```

    ## Crocc2 
    ##    0.3

``` r
#We chose an adjusted R2 of 0.3, as less than 1% of the genes in the pancreas permutations achieved this threshold.

# Histogram of all R^2
pdf("res_pancreas.r2.pdf")
hist(res_pancreas.r2, ylim=c(0,1000), xlim=c(0,1), xlab=expression(paste("Adjusted R"^"2")),main="Pancreas")
abline(v=round(sort(res_pancreas.r2, decreasing=T)[round(.01 * length(res_pancreas.r2))],1), lty=2)
dev.off()
```

    ## png 
    ##   2

## 3.4 Top Dynamic Genes

``` r
# #BRAIN
# res.top_brain_cov <- topTrendy(res_brain_cov, adjR2Cut = 0.2)
# 
# #HEART
# res.top_heart_cov <- topTrendy(res_heart_cov, adjR2Cut = 0.2)
# 
# #MUSCLE
# res.top_muscle_cov <- topTrendy(res_muscle_cov, adjR2Cut = 0.3)
# 
# #LIVER
# res.top_liver_cov <- topTrendy(res_liver_cov, adjR2Cut = 0.1)
# 
# #PANCREAS
# res.top_pancreas_cov <- topTrendy(res_pancreas_cov, adjR2Cut = 0.3)

#save top dynamic genes
# write.table(res.top_brain_cov$AdjustedR2 , "res.top_brain_cov.txt", quote=F, row.names=T)
# write.table(res.top_heart_cov$AdjustedR2 , "res.top_heart_cov.txt", quote=F, row.names=T)
# write.table(res.top_muscle_cov$AdjustedR2 , "res.top_muscle_cov.txt", quote=F, row.names=T)
# write.table(res.top_liver_cov$AdjustedR2 , "res.top_liver_cov.txt", quote=F, row.names=T)
# write.table(res.top_pancreas_cov$AdjustedR2 , "res.top_pancreas_cov.txt", quote=F, row.names=T)

#save R objects to load next time
# saveRDS(res.top_brain_cov, file = "res.top_brain_cov.rds")
# 
# saveRDS(res.top_heart_cov, file = "res.top_heart_cov.rds")
# 
# saveRDS(res.top_muscle_cov, file = "res.top_muscle_cov.rds")
# 
# saveRDS(res.top_liver_cov, file = "res.top_liver_cov.rds")
# 
# saveRDS(res.top_pancreas_cov, file = "res.top_pancreas_cov.rds")

#restore the objects
res.top_brain_cov <- readRDS(file = "res.top_brain_cov.rds")

res.top_heart_cov <- readRDS(file = "res.top_heart_cov.rds")

res.top_muscle_cov <- readRDS(file = "res.top_muscle_cov.rds")

res.top_liver_cov <- readRDS(file = "res.top_liver_cov.rds")

res.top_pancreas_cov <- readRDS(file = "res.top_pancreas_cov.rds")
```

### 3.4.1 Summary tables

``` r
#Brain
brain_table <- formatResults(res.top_brain_cov)
length(which(grepl("Pvalue",names(brain_table)))) # no. of segments found in the results
```

    ## [1] 5

``` r
brain_table_sig <- brain_table[which(brain_table$Segment1.Pvalue<0.1 | brain_table$Segment2.Pvalue<0.1 | brain_table$Segment3.Pvalue<0.1 | brain_table$Segment4.Pvalue<0.1 | brain_table$Segment5.Pvalue<0.1), ]

write.csv(brain_table , "brain_table.csv", quote=F, row.names=T)
write.csv(brain_table_sig , "brain_table_sig.csv", quote=F, row.names=T)

#Heart
heart_table <- formatResults(res.top_heart_cov)
length(which(grepl("Pvalue",names(heart_table)))) # no. of segments found in the results
```

    ## [1] 6

``` r
heart_table_sig <- heart_table[which(heart_table$Segment1.Pvalue<0.1 | heart_table$Segment2.Pvalue<0.1 | heart_table$Segment3.Pvalue<0.1 | heart_table$Segment4.Pvalue<0.1 | heart_table$Segment5.Pvalue<0.1 | heart_table$Segment6.Pvalue<0.1), ]

write.csv(heart_table , "heart_table.csv", quote=F, row.names=T)
write.csv(heart_table_sig , "heart_table_sig.csv", quote=F, row.names=T)

#Muscle
muscle_table <- formatResults(res.top_muscle_cov)
length(which(grepl("Pvalue",names(muscle_table)))) # no. of segments found in the results
```

    ## [1] 7

``` r
muscle_table_sig <- muscle_table[which(muscle_table$Segment1.Pvalue<0.1 | muscle_table$Segment2.Pvalue<0.1 | muscle_table$Segment3.Pvalue<0.1 | muscle_table$Segment4.Pvalue<0.1 | muscle_table$Segment5.Pvalue<0.1 | muscle_table$Segment6.Pvalue<0.1 | muscle_table$Segment7.Pvalue<0.1), ]

write.csv(muscle_table , "muscle_table.csv", quote=F, row.names=T)
write.csv(muscle_table_sig , "muscle_table_sig.csv", quote=F, row.names=T)

#Liver
liver_table <- formatResults(res.top_liver_cov)
length(which(grepl("Pvalue",names(liver_table)))) # no. of segments found in the results
```

    ## [1] 8

``` r
liver_table_sig <- liver_table[which(liver_table$Segment1.Pvalue<0.1 | liver_table$Segment2.Pvalue<0.1 | liver_table$Segment3.Pvalue<0.1 | liver_table$Segment4.Pvalue<0.1 | liver_table$Segment5.Pvalue<0.1 | liver_table$Segment6.Pvalue<0.1 | liver_table$Segment7.Pvalue<0.1 | liver_table$Segment8.Pvalue<0.1), ]

write.csv(liver_table , "liver_table.csv", quote=F, row.names=T)
write.csv(liver_table_sig , "liver_table_sig.csv", quote=F, row.names=T)

#Pancreas
pancreas_table <- formatResults(res.top_pancreas_cov)
length(which(grepl("Pvalue",names(pancreas_table)))) # no. of segments found in the results
```

    ## [1] 8

``` r
pancreas_table_sig <- pancreas_table[which(pancreas_table$Segment1.Pvalue<0.1 | pancreas_table$Segment2.Pvalue<0.1 | pancreas_table$Segment3.Pvalue<0.1 | pancreas_table$Segment4.Pvalue<0.1 | pancreas_table$Segment5.Pvalue<0.1 | pancreas_table$Segment6.Pvalue<0.1 | pancreas_table$Segment7.Pvalue<0.1 | pancreas_table$Segment8.Pvalue<0.1), ]

write.csv(pancreas_table , "pancreas_table.csv", quote=F, row.names=T)
write.csv(pancreas_table_sig , "pancreas_table_sig.csv", quote=F, row.names=T)
```

### Fig 1B

``` r
#calculate the percentage of top dynamic genes relative to the total of expressed genes in each tissue (obtained after filtering of low expressed genes - section 3.2 of Mouse_Aging_DESeq2.Rmd)
brain=round(((nrow(brain_table_sig)/length(res_brain_cov))*100),1)
heart=round(((nrow(heart_table_sig)/length(res_heart_cov))*100),1)
muscle=round(((nrow(muscle_table_sig)/length(res_muscle_cov))*100),1)
liver=round(((nrow(liver_table_sig)/length(res_liver_cov))*100),1)
pancreas=round(((nrow(pancreas_table_sig)/length(res_pancreas_cov))*100),1)

#prepare columns of data frame to input to ggplot
tissue <- c('Brain','Heart','Muscle','Liver','Pancreas')
degsPercent <- c(brain,heart,muscle,liver,pancreas)

#create data frame to input to ggplot
data <- data.frame(tissue, degsPercent)

gg_dot <- data %>%
  ggplot() +
  geom_point(aes(x = degsPercent, y = tissue), size =25, col="grey50",position="identity") + labs(x="Top dynamic genes (%)") +
  geom_text(aes(x = degsPercent, y = tissue, label = degsPercent), col = "white", size=8) +
  theme_light() +
  theme(axis.title=element_text(size=20,face="bold"), 
        axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title =element_text(size=7,face="bold"), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="bottom", 
        axis.text = element_blank(), 
        axis.ticks=element_blank(), 
        axis.title.y=element_blank()) + 
  facet_grid(rows=vars(tissue), scales="free") +
  scale_x_continuous(limits=c(0,20))

gg_dot

g <- ggplot_gtable(ggplot_build(gg_dot))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#77AADD","#EE8866","#EEDD88","#FFAABB","#44BB99")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
```

![](segmented_regression_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
tiff('Fig-1B.tiff', units="cm", width=25, height=20, res=300, compression = 'lzw')
grid.draw(g)
dev.off()
```

    ## png 
    ##   2

## 3.5 Test differences in median breakpoint distributions

``` r
allTimepoints <- as.character(seq(3, 27, 1))

#BRAIN
#prepare data for plotting
which(grepl("Breakpoint",names(brain_table_sig)))
```

    ## [1] 17 18 19 20

``` r
res.bp_brain <- table(round(unlist(brain_table_sig[,17:20]), 0)) 
res.bp_brain <- rep(0, length(allTimepoints))
names(res.bp_brain) <- allTimepoints
res.bp_brain[names(res.bp_brain)] <- res.bp_brain

#prepare data for statistical testing
res.bp_brain2 <- c(rep(3,as.vector(res.bp_brain[1])),rep(4,as.vector(res.bp_brain[2])),rep(5,as.vector(res.bp_brain[3])),rep(6,as.vector(res.bp_brain[4])),rep(7,as.vector(res.bp_brain[5])),rep(8,as.vector(res.bp_brain[6])),rep(9,as.vector(res.bp_brain[7])),rep(10,as.vector(res.bp_brain[8])),rep(11,as.vector(res.bp_brain[9])),rep(12,as.vector(res.bp_brain[10])),rep(13,as.vector(res.bp_brain[11])),rep(14,as.vector(res.bp_brain[12])),rep(15,as.vector(res.bp_brain[13])),rep(16,as.vector(res.bp_brain[14])),rep(17,as.vector(res.bp_brain[15])),rep(18,as.vector(res.bp_brain[16])),rep(19,as.vector(res.bp_brain[17])),rep(20,as.vector(res.bp_brain[18])),rep(21,as.vector(res.bp_brain[9])),rep(22,as.vector(res.bp_brain[20])),rep(23,as.vector(res.bp_brain[21])),rep(24,as.vector(res.bp_brain[22])),rep(25,as.vector(res.bp_brain[23])),rep(26,as.vector(res.bp_brain[24])),rep(27,as.vector(res.bp_brain[25])))
length(res.bp_brain2)
```

    ## [1] 0

``` r
#HEART
#ploting
which(grepl("Breakpoint",names(heart_table_sig)))
```

    ## [1] 20 21 22 23 24

``` r
heart_bp <- table(round(unlist(heart_table_sig[,20:24]), 0))
res.bp_heart <- rep(0, length(allTimepoints))
names(res.bp_heart) <- allTimepoints
res.bp_heart[names(heart_bp)] <- heart_bp

#statistical testing
res.bp_heart2 <- c(rep(3,as.vector(res.bp_heart[1])),rep(4,as.vector(res.bp_heart[2])),rep(5,as.vector(res.bp_heart[3])),rep(6,as.vector(res.bp_heart[4])),rep(7,as.vector(res.bp_heart[5])),rep(8,as.vector(res.bp_heart[6])),rep(9,as.vector(res.bp_heart[7])),rep(10,as.vector(res.bp_heart[8])),rep(11,as.vector(res.bp_heart[9])),rep(12,as.vector(res.bp_heart[10])),rep(13,as.vector(res.bp_heart[11])),rep(14,as.vector(res.bp_heart[12])),rep(15,as.vector(res.bp_heart[13])),rep(16,as.vector(res.bp_heart[14])),rep(17,as.vector(res.bp_heart[15])),rep(18,as.vector(res.bp_heart[16])),rep(19,as.vector(res.bp_heart[17])),rep(20,as.vector(res.bp_heart[18])),rep(21,as.vector(res.bp_heart[9])),rep(22,as.vector(res.bp_heart[20])),rep(23,as.vector(res.bp_heart[21])),rep(24,as.vector(res.bp_heart[22])),rep(25,as.vector(res.bp_heart[23])),rep(26,as.vector(res.bp_heart[24])),rep(27,as.vector(res.bp_heart[25])))
length(res.bp_heart2)
```

    ## [1] 535

``` r
#LIVER
#plotting
which(grepl("Breakpoint",names(liver_table_sig)))
```

    ## [1] 26 27 28 29 30 31 32

``` r
liver_bp <- table(round(unlist(liver_table_sig[,26:32]), 0))
res.bp_liver <- rep(0, length(allTimepoints))
names(res.bp_liver) <- allTimepoints
res.bp_liver[names(liver_bp)] <- liver_bp

#statistical testing
res.bp_liver2 <- c(rep(3,as.vector(res.bp_liver[1])),rep(4,as.vector(res.bp_liver[2])),rep(5,as.vector(res.bp_liver[3])),rep(6,as.vector(res.bp_liver[4])),rep(7,as.vector(res.bp_liver[5])),rep(8,as.vector(res.bp_liver[6])),rep(9,as.vector(res.bp_liver[7])),rep(10,as.vector(res.bp_liver[8])),rep(11,as.vector(res.bp_liver[9])),rep(12,as.vector(res.bp_liver[10])),rep(13,as.vector(res.bp_liver[11])),rep(14,as.vector(res.bp_liver[12])),rep(15,as.vector(res.bp_liver[13])),rep(16,as.vector(res.bp_liver[14])),rep(17,as.vector(res.bp_liver[15])),rep(18,as.vector(res.bp_liver[16])),rep(19,as.vector(res.bp_liver[17])),rep(20,as.vector(res.bp_liver[18])),rep(21,as.vector(res.bp_liver[9])),rep(22,as.vector(res.bp_liver[20])),rep(23,as.vector(res.bp_liver[21])),rep(24,as.vector(res.bp_liver[22])),rep(25,as.vector(res.bp_liver[23])),rep(26,as.vector(res.bp_liver[24])),rep(27,as.vector(res.bp_liver[25])))
length(res.bp_liver2)
```

    ## [1] 651

``` r
#MUSCLE
#plotting
which(grepl("Breakpoint",names(muscle_table_sig)))
```

    ## [1] 23 24 25 26 27 28

``` r
muscle_bp <- table(round(unlist(muscle_table_sig[,23:28]), 0))
res.bp_muscle <- rep(0, length(allTimepoints))
names(res.bp_muscle) <- allTimepoints
res.bp_muscle[names(muscle_bp)] <- muscle_bp

#statistical testing
res.bp_muscle2 <- c(rep(3,as.vector(res.bp_muscle[1])),rep(4,as.vector(res.bp_muscle[2])),rep(5,as.vector(res.bp_muscle[3])),rep(6,as.vector(res.bp_muscle[4])),rep(7,as.vector(res.bp_muscle[5])),rep(8,as.vector(res.bp_muscle[6])),rep(9,as.vector(res.bp_muscle[7])),rep(10,as.vector(res.bp_muscle[8])),rep(11,as.vector(res.bp_muscle[9])),rep(12,as.vector(res.bp_muscle[10])),rep(13,as.vector(res.bp_muscle[11])),rep(14,as.vector(res.bp_muscle[12])),rep(15,as.vector(res.bp_muscle[13])),rep(16,as.vector(res.bp_muscle[14])),rep(17,as.vector(res.bp_muscle[15])),rep(18,as.vector(res.bp_muscle[16])),rep(19,as.vector(res.bp_muscle[17])),rep(20,as.vector(res.bp_muscle[18])),rep(21,as.vector(res.bp_muscle[9])),rep(22,as.vector(res.bp_muscle[20])),rep(23,as.vector(res.bp_muscle[21])),rep(24,as.vector(res.bp_muscle[22])),rep(25,as.vector(res.bp_muscle[23])),rep(26,as.vector(res.bp_muscle[24])),rep(27,as.vector(res.bp_muscle[25])))
length(res.bp_muscle2)
```

    ## [1] 769

``` r
#PANCREAS
#plotting
which(grepl("Breakpoint",names(pancreas_table_sig)))
```

    ## [1] 26 27 28 29 30 31 32

``` r
pancreas_bp <- table(round(unlist(pancreas_table_sig[,26:32]), 0))
res.bp_pancreas <- rep(0, length(allTimepoints))
names(res.bp_pancreas) <- allTimepoints
res.bp_pancreas[names(pancreas_bp)] <- pancreas_bp

#statistical testing
res.bp_pancreas2 <- c(rep(3,as.vector(res.bp_pancreas[1])),rep(4,as.vector(res.bp_pancreas[2])),rep(5,as.vector(res.bp_pancreas[3])),rep(6,as.vector(res.bp_pancreas[4])),rep(7,as.vector(res.bp_pancreas[5])),rep(8,as.vector(res.bp_pancreas[6])),rep(9,as.vector(res.bp_pancreas[7])),rep(10,as.vector(res.bp_pancreas[8])),rep(11,as.vector(res.bp_pancreas[9])),rep(12,as.vector(res.bp_pancreas[10])),rep(13,as.vector(res.bp_pancreas[11])),rep(14,as.vector(res.bp_pancreas[12])),rep(15,as.vector(res.bp_pancreas[13])),rep(16,as.vector(res.bp_pancreas[14])),rep(17,as.vector(res.bp_pancreas[15])),rep(18,as.vector(res.bp_pancreas[16])),rep(19,as.vector(res.bp_pancreas[17])),rep(20,as.vector(res.bp_pancreas[18])),rep(21,as.vector(res.bp_pancreas[9])),rep(22,as.vector(res.bp_pancreas[20])),rep(23,as.vector(res.bp_pancreas[21])),rep(24,as.vector(res.bp_pancreas[22])),rep(25,as.vector(res.bp_pancreas[23])),rep(26,as.vector(res.bp_pancreas[24])),rep(27,as.vector(res.bp_pancreas[25])))
length(res.bp_pancreas2)
```

    ## [1] 652

``` r
#identify group differences in breakpoint distribution
tissue <- c(rep("Brain",length(res.bp_brain2)),rep("Heart",length(res.bp_heart2)),rep("Liver",length(res.bp_liver2)),rep("Muscle",length(res.bp_muscle2)),rep("Pancreas",length(res.bp_pancreas2)))

data <- data.frame(tissue,c(res.bp_brain2,res.bp_heart2,res.bp_liver2,res.bp_muscle2,res.bp_pancreas2))
names(data) <- c("Tissue","Breakpoints")

kruskal.test(data$Breakpoints~data$Tissue, data=data)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  data$Breakpoints by data$Tissue
    ## Kruskal-Wallis chi-squared = 178.22, df = 3, p-value < 2.2e-16

``` r
dunnTest(data$Breakpoints~data$Tissue, data=data, method='bh')
```

    ##          Comparison           Z      P.unadj        P.adj
    ## 1     Heart - Liver  -0.9360464 3.492493e-01 3.492493e-01
    ## 2    Heart - Muscle -11.1321526 8.749667e-29 5.249800e-28
    ## 3    Liver - Muscle -10.7419769 6.464511e-27 1.939353e-26
    ## 4  Heart - Pancreas  -1.9417637 5.216572e-02 7.824857e-02
    ## 5  Liver - Pancreas  -1.0585290 2.898144e-01 3.477772e-01
    ## 6 Muscle - Pancreas   9.6447731 5.172543e-22 1.034509e-21

``` r
#identify median breakpoint distribution
median(data$Breakpoints[data$Tissue=="Brain"])
```

    ## [1] NA

``` r
median(data$Breakpoints[data$Tissue=="Heart"])
```

    ## [1] 14

``` r
median(data$Breakpoints[data$Tissue=="Liver"])
```

    ## [1] 12

``` r
median(data$Breakpoints[data$Tissue=="Muscle"])
```

    ## [1] 18

``` r
median(data$Breakpoints[data$Tissue=="Pancreas"])
```

    ## [1] 15

``` r
#calculate the %of breakpoints until median
tissue <- c(rep("Brain",25),rep("Heart",25),rep("Liver",25),rep("Muscle",25),rep("Pancreas",25))
data2=data.frame(tissue,c(names(res.bp_brain),names(res.bp_heart),names(res.bp_liver),names(res.bp_muscle),names(res.bp_pancreas)),c(as.vector(res.bp_brain),as.vector(res.bp_heart),as.vector(res.bp_liver),as.vector(res.bp_muscle),as.vector(res.bp_pancreas)))
names(data2) <- c("Tissue","Breakpoint","Frequency")
data2 <- ddply(data2, .(Tissue), transform, relFreq = (Frequency/sum(Frequency)*100))

#cumulative percentage to midlle age
brain_percentToMiddle <- sum(data2$relFreq[c(1:13)])
heart_percentToMiddle <- sum(data2$relFreq[c(26:38)])
liver_percentToMiddle <- sum(data2$relFreq[c(51:63)])
muscle_percentToMiddle <- sum(data2$relFreq[c(76:88)])
pancreas_percentToMiddle <- sum(data2$relFreq[c(101:113)])

brain_percentToMiddle
```

    ## [1] NaN

``` r
heart_percentToMiddle
```

    ## [1] 68.39187

``` r
liver_percentToMiddle
```

    ## [1] 65.14032

``` r
muscle_percentToMiddle
```

    ## [1] 41.875

``` r
pancreas_percentToMiddle
```

    ## [1] 79.09774

### Fig 1C

``` r
#Breakpoint distribution histogram (based on Trendy's original function breakpointDist (https://github.com/rhondabacher/Trendy/blob/RELEASE_3_10/R/breakpointDist.R))

#plot all tissues together
data <- data.frame("timepoint"=c(names(res.bp_brain),names(res.bp_heart),names(res.bp_muscle),names(res.bp_liver),names(res.bp_pancreas)),"breakpoints"=c(res.bp_brain,res.bp_heart,res.bp_muscle,res.bp_liver,res.bp_pancreas),"tissue"=rep(c("Brain","Heart","Muscle","Liver","Pancreas"),each=25))

data$timepoint <- factor(data$timepoint, levels=c(3:27))

data1 <- data[data$tissue=="Brain",]
data2 <- data[data$tissue=="Heart",]
data3 <- data[data$tissue=="Liver",]
data4 <- data[data$tissue=="Muscle",]
data5 <- data[data$tissue=="Pancreas",]

hist1 <- ggplot(data1, aes(timepoint,breakpoints)) +
  geom_bar(position="identity", stat="identity") +
  geom_vline(xintercept = 10, linetype="dotted", color = "#77AADD", size=1.5) + #because the x axis starts in time point 3, the xintercept need to be the median values - 2
  annotate("segment", x = 0, xend = 13, y = 35, yend = 35,colour = "grey40", size=1.2)+  annotate("text", x=(0 + 13)/2, y = 31, label = "64.5% breakpoints", colour="grey40",size=6) +
  facet_grid(cols =vars(tissue)) +
  theme_light() +
  labs(x=" ", y=" \n ") +
  theme(title=element_text(size=25,face="bold"),axis.text=element_text(size=15, face="bold"),axis.title=element_text(size=25,face="bold"),        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
legend.title =element_blank(), legend.text = element_text(size=20),
axis.text.x = element_text(hjust = 0.5, face="bold"),
strip.text = element_text(size=25, face="bold"))

#rebuild the plot without the legend
hist1 <- hist1 + theme(legend.position = "none")

hist1
```

![](segmented_regression_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
g1 <- ggplot_gtable(ggplot_build(hist1))
stript <- which(grepl('strip-t', g1$layout$name))
fills <- c("#77AADD")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

hist2 <- ggplot(data2, aes(timepoint,breakpoints)) +
  geom_bar(position="identity", stat="identity") +
  geom_vline(xintercept = 12, linetype="dotted", color = "#EE8866", size=1.5) +
  annotate("segment", x = 0, xend = 13, y = 90, yend = 90,colour = "grey40", size=1.2)+  annotate("text", x=(0 + 13)/2, y = 80, label = "68.4% breakpoints", colour="grey40",size=6) +
  facet_grid(cols =vars(tissue)) +
  theme_light() +
  labs(x=" ", y=" \n ") +
  theme(title=element_text(size=25,face="bold"),axis.text=element_text(size=15, face="bold"),axis.title=element_text(size=25,face="bold"),        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
legend.title =element_blank(), legend.text = element_text(size=20),
axis.text.x = element_text(hjust = 0.5, face="bold"),
strip.text = element_text(size=25, face="bold"))

g2 <- ggplot_gtable(ggplot_build(hist2))
stript <- which(grepl('strip-t', g2$layout$name))
fills <- c("#EE8866")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

hist3 <- ggplot(data3, aes(timepoint,breakpoints)) +
  geom_bar(position="identity", stat="identity") +
  geom_vline(xintercept = 10, linetype="dotted", color = "#EEDD88", size=1.5) +
  annotate("segment", x = 0, xend = 13, y = 90, yend = 90,colour = "grey40", size=1.2)+  annotate("text", x=(0 + 13)/2, y = 80, label = "65.1% breakpoints", colour="grey40",size=6) +
  facet_grid(cols =vars(tissue)) +
  theme_light() +
  labs(x=" ", y=" \n ") +
  theme(title=element_text(size=25,face="bold"),axis.text=element_text(size=15, face="bold"),axis.title=element_text(size=25,face="bold"),        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
legend.title =element_blank(), legend.text = element_text(size=20),
axis.text.x = element_text(hjust = 0.5, face="bold"),
strip.text = element_text(size=25, face="bold"))

g3 <- ggplot_gtable(ggplot_build(hist3))
stript <- which(grepl('strip-t', g3$layout$name))
fills <- c("#EEDD88")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g3$grobs[[i]]$grobs[[1]]$childrenOrder))
  g3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

hist4 <- ggplot(data4, aes(timepoint,breakpoints)) +
  geom_bar(position="identity", stat="identity") +
  geom_vline(xintercept =16, linetype="dotted", color = "#FFAABB", size=1.5) +
  annotate("segment", x = 0, xend = 13, y = 250, yend = 250,colour = "grey40", size=1.2)+  annotate("text", x=(0 + 13)/2, y = 220, label = "41.9% breakpoints", colour="grey40",size=6) +
  facet_grid(cols =vars(tissue)) +
  theme_light() +
  labs(x=" ", y=" \n ") +
  theme(title=element_text(size=25,face="bold"),axis.text=element_text(size=15, face="bold"),axis.title=element_text(size=25,face="bold"),        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
legend.title =element_blank(), legend.text = element_text(size=20),
axis.text.x = element_text(hjust = 0.5, face="bold"),
strip.text = element_text(size=25, face="bold"))

g4 <- ggplot_gtable(ggplot_build(hist4))
stript <- which(grepl('strip-t', g4$layout$name))
fills <- c("#FFAABB")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g4$grobs[[i]]$grobs[[1]]$childrenOrder))
  g4$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

hist5 <- ggplot(data5, aes(timepoint,breakpoints)) +
  geom_bar(position="identity", stat="identity") +
  geom_vline(xintercept = 13, linetype="dotted", color = "#44BB99", size=1.5) +
  annotate("segment", x = 0, xend = 13, y = 225, yend = 225,colour = "grey40", size=1.2)+  annotate("text", x=(0 + 13)/2, y = 195, label = "79.1% breakpoints", colour="grey40",size=6) +
  facet_grid(cols =vars(tissue)) +
  theme_light() +
  labs(x=" ", y=" \n ") +
  theme(title=element_text(size=25,face="bold"),axis.text=element_text(size=15, face="bold"),axis.title=element_text(size=25,face="bold"),        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
legend.title =element_blank(), legend.text = element_text(size=20),
axis.text.x = element_text(hjust = 0.5, face="bold"),
strip.text = element_text(size=25, face="bold"))

g5 <- ggplot_gtable(ggplot_build(hist5))
stript <- which(grepl('strip-t', g5$layout$name))
fills <- c("#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g5$grobs[[i]]$grobs[[1]]$childrenOrder))
  g5$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff('Fig-1C.tiff', units="cm", width=30, height=40, res=300, compression = 'lzw')
grid::grid.newpage()
ggarrange(g1,g2,g3,g4,g5, nrow=5)
grid.text(label = "Lifespan time point", x = unit(0.55, "npc"), y = unit(0.01, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "No. of genes  with breakpoints in gene expression", x = unit(0.05, "npc"), y = unit(0.5, "npc"),just = "centre",rot=90, gp=gpar(fontsize=20,fontface="bold"))
dev.off()
```

    ## png 
    ##   2

### Fig 1D

``` r
#Calculate the number of significant genes with min (0) to max(7) breakpoints
brain_7 <- sum(!is.na(brain_table_sig$Breakpoint7))
brain_6 <- sum(!is.na(brain_table_sig$Breakpoint6))
brain_5 <- sum(!is.na(brain_table_sig$Breakpoint5))
brain_4 <- sum(!is.na(brain_table_sig$Breakpoint4))
brain_3 <- sum(!is.na(brain_table_sig$Breakpoint3))
brain_2 <- sum(!is.na(brain_table_sig$Breakpoint2))
brain_1 <- sum(!is.na(brain_table_sig$Breakpoint1))
brain_0 <- sum(is.na(brain_table_sig$Breakpoint1)) #if a gene is increasingly up- or down-regulated it will not exhibit breakpoints thus it will appear as NA already in the first breakpoint

heart_7 <- sum(!is.na(heart_table_sig$Breakpoint7))
heart_6 <- sum(!is.na(heart_table_sig$Breakpoint6))
heart_5 <- sum(!is.na(heart_table_sig$Breakpoint5))
heart_4 <- sum(!is.na(heart_table_sig$Breakpoint4))
heart_3 <- sum(!is.na(heart_table_sig$Breakpoint3))
heart_2 <- sum(!is.na(heart_table_sig$Breakpoint2))
heart_1 <- sum(!is.na(heart_table_sig$Breakpoint1))
heart_0 <- sum(is.na(heart_table_sig$Breakpoint1)) #if a gene is increasingly up- or down-regulated it will not exhibit breakpoints thus it will appear as NA already in the first breakpoint

muscle_7 <- sum(!is.na(muscle_table_sig$Breakpoint7))
muscle_6 <- sum(!is.na(muscle_table_sig$Breakpoint6))
muscle_5 <- sum(!is.na(muscle_table_sig$Breakpoint5))
muscle_4 <- sum(!is.na(muscle_table_sig$Breakpoint4))
muscle_3 <- sum(!is.na(muscle_table_sig$Breakpoint3))
muscle_2 <- sum(!is.na(muscle_table_sig$Breakpoint2))
muscle_1 <- sum(!is.na(muscle_table_sig$Breakpoint1))
muscle_0 <- sum(is.na(muscle_table_sig$Breakpoint1)) #if a gene is increasingly up- or down-regulated it will not exhibit breakpoints thus it will appear as NA already in the first breakpoint

liver_7 <- sum(!is.na(liver_table_sig$Breakpoint7))
liver_6 <- sum(!is.na(liver_table_sig$Breakpoint6))
liver_5 <- sum(!is.na(liver_table_sig$Breakpoint5))
liver_4 <- sum(!is.na(liver_table_sig$Breakpoint4))
liver_3 <- sum(!is.na(liver_table_sig$Breakpoint3))
liver_2 <- sum(!is.na(liver_table_sig$Breakpoint2))
liver_1 <- sum(!is.na(liver_table_sig$Breakpoint1))
liver_0 <- sum(is.na(liver_table_sig$Breakpoint1)) #if a gene is increasingly up- or down-regulated it will not exhibit breakpoints thus it will appear as NA already in the first breakpoint

pancreas_7 <- sum(!is.na(pancreas_table_sig$Breakpoint7))
pancreas_6 <- sum(!is.na(pancreas_table_sig$Breakpoint6))
pancreas_5 <- sum(!is.na(pancreas_table_sig$Breakpoint5))
pancreas_4 <- sum(!is.na(pancreas_table_sig$Breakpoint4))
pancreas_3 <- sum(!is.na(pancreas_table_sig$Breakpoint3))
pancreas_2 <- sum(!is.na(pancreas_table_sig$Breakpoint2))
pancreas_1 <- sum(!is.na(pancreas_table_sig$Breakpoint1))
pancreas_0 <- sum(is.na(pancreas_table_sig$Breakpoint1)) #if a gene is increasingly up- or down-regulated it will not exhibit breakpoints thus it will appear as NA already in the first breakpoint

#prepare columns of data frame to input to ggplot
tissue=c(rep("Brain",8),rep("Heart",8),rep("Muscle",8),rep("Liver",8),rep("Pancreas",8)) 
breakpoints=rep(c("k=7","k=6","k=5","k=4", "k=3", "k=2", "k=1", "k=0") , 5)
frequency=c(brain_7,brain_6,brain_5,brain_4,brain_3,brain_2,brain_1,brain_0,
            heart_7,heart_6,heart_5,heart_4,heart_3,heart_2,heart_1,heart_0,
            muscle_7,muscle_6,muscle_5,muscle_4,muscle_3,muscle_2,muscle_1,muscle_0,
            liver_7,liver_6,liver_5,liver_4,liver_3,liver_2,liver_1,liver_0,
            pancreas_7,pancreas_6,pancreas_5,pancreas_4,pancreas_3,pancreas_2,pancreas_1,pancreas_0)

#create data frame to input to ggplot
data=data.frame(tissue,breakpoints,frequency)
data <- ddply(data, .(tissue), transform, relFreq = (frequency/sum(frequency)*100))

data1 <- data[data$tissue=="Brain",]
data2 <- data[data$tissue=="Heart",]
data3 <- data[data$tissue=="Liver",]
data4 <- data[data$tissue=="Muscle",]
data5 <- data[data$tissue=="Pancreas",]


bk1 <- ggplot(data1, aes(x="", y=relFreq, fill=breakpoints)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_text(size=15,face="bold"), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="right") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("% of genes exhibiting", values = c("#f9d5e5", "#eeac99", "#e06377", "#c83349", "#5b9aa0", "#d6d4e0", "#b8a9c9", "#622569")) 

#retrieve legend from first plot to serve as common legend in the final plot
mylegend<-get_legend(bk1)

#rebuild the plot without the legend
bk1 <- bk1 + theme(legend.position = "none")

#create pie
bk_pie1 <- bk1 + coord_polar("y", start=0) 

bk_pie1
```

![](segmented_regression_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
g6 <- ggplot_gtable(ggplot_build(bk_pie1))
stript <- which(grepl('strip-t', g6$layout$name))
fills <- c("#77AADD")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g6$grobs[[i]]$grobs[[1]]$childrenOrder))
  g6$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk2 <- ggplot(data2, aes(x="", y=relFreq, fill=breakpoints)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("No. of genes exhibiting:", values = c("#f9d5e5", "#eeac99", "#e06377", "#c83349", "#5b9aa0", "#d6d4e0", "#b8a9c9", "#622569")) 

bk_pie2 <- bk2 + coord_polar("y", start=0) 

g7 <- ggplot_gtable(ggplot_build(bk_pie2))
stript <- which(grepl('strip-t', g7$layout$name))
fills <- c("#EE8866")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g7$grobs[[i]]$grobs[[1]]$childrenOrder))
  g7$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk3 <- ggplot(data3, aes(x="", y=relFreq, fill=breakpoints)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("No. of genes exhibiting:", values = c("#f9d5e5", "#eeac99", "#e06377", "#c83349", "#5b9aa0", "#d6d4e0", "#b8a9c9", "#622569")) 

bk_pie3 <- bk3 + coord_polar("y", start=0) 

g8 <- ggplot_gtable(ggplot_build(bk_pie3))
stript <- which(grepl('strip-t', g8$layout$name))
fills <- c("#EEDD88")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g8$grobs[[i]]$grobs[[1]]$childrenOrder))
  g8$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk4 <- ggplot(data4, aes(x="", y=relFreq, fill=breakpoints)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("No. of genes exhibiting:", values = c("#f9d5e5", "#eeac99", "#e06377", "#c83349", "#5b9aa0", "#d6d4e0", "#b8a9c9", "#622569"))  

bk_pie4 <- bk4 + coord_polar("y", start=0) 

g9 <- ggplot_gtable(ggplot_build(bk_pie4))
stript <- which(grepl('strip-t', g9$layout$name))
fills <- c("#FFAABB")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g9$grobs[[i]]$grobs[[1]]$childrenOrder))
  g9$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk5 <- ggplot(data5, aes(x="", y=relFreq, fill=breakpoints)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("No. of genes exhibiting:", values = c("#f9d5e5", "#eeac99", "#e06377", "#c83349", "#5b9aa0", "#d6d4e0", "#b8a9c9", "#622569"))  

bk_pie5 <- bk5 + coord_polar("y", start=0) 

g10 <- ggplot_gtable(ggplot_build(bk_pie5))
stript <- which(grepl('strip-t', g10$layout$name))
fills <- c("#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g10$grobs[[i]]$grobs[[1]]$childrenOrder))
  g10$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff('Fig-1D.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
grid::grid.newpage()
ggarrange(g6,g7,g8,g9,g10,mylegend, ncol=3, nrow=2)
dev.off()
```

    ## png 
    ##   2

## 3.6 Test biotype enrichment

``` r
annotation <- read.table("annotation.txt", header = T) #input biotype annotation
names(annotation) <- c("Feature","Biotype") # update colnames to match top dynamic tables feature columns

brain_table_sig <- merge(brain_table_sig, annotation, by="Feature", all=FALSE)
nrow(brain_table_sig)
```

    ## [1] 747

``` r
heart_table_sig <- merge(heart_table_sig, annotation, by="Feature", all=FALSE)
nrow(heart_table_sig)
```

    ## [1] 1799

``` r
liver_table_sig <- merge(liver_table_sig, annotation, by="Feature", all=FALSE)
nrow(liver_table_sig)
```

    ## [1] 3690

``` r
muscle_table_sig <- merge(muscle_table_sig, annotation, by="Feature", all=FALSE)
nrow(muscle_table_sig)
```

    ## [1] 1048

``` r
pancreas_table_sig <- merge(pancreas_table_sig, annotation, by="Feature", all=FALSE)
nrow(pancreas_table_sig)
```

    ## [1] 336

``` r
nrow(annotation) #total no. of annotated genes in the reference genome
```

    ## [1] 49512

``` r
#retrieve the proportion of protein coding genes, lncRNAs and other biotypes from the reference genome (GRCm39)
annotation_pc <- sum(annotation$Biotype=="protein_coding")
annotation_pc #total no. of annotated protein coding genes in the reference genome
```

    ## [1] 21177

``` r
annotation_lncRNA <- sum(annotation$Biotype=="lncRNA")
annotation_lncRNA #total no. of annotated lncRNAs in the reference genome
```

    ## [1] 9082

``` r
#Calculate the number of significant genes in each biotype
brain_pc <- sum(brain_table_sig$Biotype=="protein_coding")
brain_lncRNA <- sum(brain_table_sig$Biotype=="lncRNA")
brain_all <- length(brain_table_sig$Biotype)
brain_other <- brain_all-brain_pc-brain_lncRNA

heart_pc <- sum(heart_table_sig$Biotype=="protein_coding")
heart_lncRNA <- sum(heart_table_sig$Biotype=="lncRNA")
heart_all <- length(heart_table_sig$Biotype)
heart_other <- heart_all-heart_pc-heart_lncRNA

liver_pc <- sum(liver_table_sig$Biotype=="protein_coding")
liver_lncRNA <- sum(liver_table_sig$Biotype=="lncRNA")
liver_all <- length(liver_table_sig$Biotype)
liver_other <- liver_all-liver_pc-liver_lncRNA

muscle_pc <- sum(muscle_table_sig$Biotype=="protein_coding")
muscle_lncRNA <- sum(muscle_table_sig$Biotype=="lncRNA")
muscle_all <- length(muscle_table_sig$Biotype)
muscle_other <- muscle_all-muscle_pc-muscle_lncRNA

pancreas_pc <- sum(pancreas_table_sig$Biotype=="protein_coding")
pancreas_lncRNA <- sum(pancreas_table_sig$Biotype=="lncRNA")
pancreas_all <- length(pancreas_table_sig$Biotype)
pancreas_other <- pancreas_all-pancreas_pc-pancreas_lncRNA

#create data frame with tissue-specific annotation information
tissue <- c("Brain","Heart","Liver","Muscle","Pancreas")
protein_coding <- c(brain_pc,heart_pc,liver_pc,muscle_pc,pancreas_pc)
lncRNA <- c(brain_lncRNA,heart_lncRNA,liver_lncRNA,muscle_lncRNA,pancreas_lncRNA)
total <- c(brain_all, heart_all, liver_all, muscle_all, pancreas_all)

annot_data <- data.frame(protein_coding,lncRNA,total)
rownames(annot_data) <- tissue
```

``` r
#Fisher's Exact Test

#PROTEIN CODING GENES
fisher_pc <- c()

for (i in 1:nrow(annot_data)) {
  M <- as.table(rbind(c(annot_data[i,1],annotation_pc), c(annot_data[i,3]-annot_data[i,1],nrow(annotation)-annotation_pc)))
  X <- chisq.test(M)
  fisher_pc[[i]] <- data.frame(pval=X$p.value)
}

fisher_pc_sig <- do.call(rbind.data.frame,fisher_pc)
fisher_pc_sig$tissue <- rownames(annot_data)

#lncRNAs
fisher_lncRNA <- c()

for (i in 1:nrow(annot_data)) {
  M <- as.table(rbind(c(annot_data[i,2],annotation_lncRNA), c(annot_data[i,3]-annot_data[i,2],nrow(annotation)-annotation_lncRNA)))
  X <- chisq.test(M)
  fisher_lncRNA[[i]] <- data.frame(pval=X$p.value)
}

fisher_lncRNA_sig <- do.call(rbind.data.frame,fisher_lncRNA)
fisher_lncRNA_sig$tissue <- rownames(annot_data)
```

### Fig 1E

``` r
#prepare columns of data frame to input to ggplot
tissue=c(rep("Brain",3),rep("Heart",3),rep("Muscle",3),rep("Liver",3),rep("Pancreas",3)) 
biotype=rep(c("Protein Coding","lncRNA","Other") , 5)
frequency=c(brain_pc,brain_lncRNA,brain_other,heart_pc,heart_lncRNA,heart_other,liver_pc,liver_lncRNA,liver_other,muscle_pc,muscle_lncRNA,muscle_other,pancreas_pc,pancreas_lncRNA,pancreas_other)

#create data frame to input to ggplot
data=data.frame(tissue,biotype,frequency)
data <- ddply(data, .(tissue), transform, relFreq = (frequency/sum(frequency)*100))

data1 <- data[data$tissue=="Brain",]
data2 <- data[data$tissue=="Heart",]
data3 <- data[data$tissue=="Liver",]
data4 <- data[data$tissue=="Muscle",]
data5 <- data[data$tissue=="Pancreas",]


bk1 <- ggplot(data1, aes(x="", y=relFreq, fill=biotype)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_text(size=15,face="bold"), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="right") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual("% of genes", values = c("#74A089", "#FDDDA0", "#F8AFA8")) 

#retrieve legend from first plot to serve as common legend in the final plot
mylegend<-get_legend(bk1)

#rebuild the plot without the legend
bk1 <- bk1 + theme(legend.position = "none")

#create pie
bk_pie1 <- bk1 + coord_polar("y", start=0) 

bk_pie1
```

![](segmented_regression_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
g6 <- ggplot_gtable(ggplot_build(bk_pie1))
stript <- which(grepl('strip-t', g6$layout$name))
fills <- c("#77AADD")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g6$grobs[[i]]$grobs[[1]]$childrenOrder))
  g6$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk2 <- ggplot(data2, aes(x="", y=relFreq, fill=biotype)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual(" ", values = c("#74A089", "#FDDDA0", "#F8AFA8")) 

bk_pie2 <- bk2 + coord_polar("y", start=0) 

g7 <- ggplot_gtable(ggplot_build(bk_pie2))
stript <- which(grepl('strip-t', g7$layout$name))
fills <- c("#EE8866")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g7$grobs[[i]]$grobs[[1]]$childrenOrder))
  g7$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk3 <- ggplot(data3, aes(x="", y=relFreq, fill=biotype)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual(" ", values = c("#74A089", "#FDDDA0", "#F8AFA8")) 

bk_pie3 <- bk3 + coord_polar("y", start=0) 

g8 <- ggplot_gtable(ggplot_build(bk_pie3))
stript <- which(grepl('strip-t', g8$layout$name))
fills <- c("#EEDD88")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g8$grobs[[i]]$grobs[[1]]$childrenOrder))
  g8$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk4 <- ggplot(data4, aes(x="", y=relFreq, fill=biotype)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual(" ", values = c("#74A089", "#FDDDA0", "#F8AFA8"))  

bk_pie4 <- bk4 + coord_polar("y", start=0) 

g9 <- ggplot_gtable(ggplot_build(bk_pie4))
stript <- which(grepl('strip-t', g9$layout$name))
fills <- c("#FFAABB")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g9$grobs[[i]]$grobs[[1]]$childrenOrder))
  g9$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

bk5 <- ggplot(data5, aes(x="", y=relFreq, fill=biotype)) +
  geom_bar(width=1, stat="identity")+
  labs(x="No. of breakpoints", y="Top dynamic genes") + theme_light() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title=element_blank(),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        strip.text = element_text(size=20, face="bold"), 
        legend.position="none") + facet_grid(cols=vars(tissue), scales="fixed") +  scale_fill_manual(" ", values = c("#74A089", "#FDDDA0", "#F8AFA8"))  

bk_pie5 <- bk5 + coord_polar("y", start=0) 

g10 <- ggplot_gtable(ggplot_build(bk_pie5))
stript <- which(grepl('strip-t', g10$layout$name))
fills <- c("#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g10$grobs[[i]]$grobs[[1]]$childrenOrder))
  g10$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff('Fig-1E.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
grid::grid.newpage()
ggarrange(g6,g7,g8,g9,g10,mylegend, ncol=3, nrow=2)
dev.off()
```

    ## png 
    ##   2

# Session info

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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] FSA_0.8.32           plyr_1.8.6           ggpubr_0.4.0        
    ##  [4] cowplot_1.1.0        ComplexHeatmap_2.4.3 ggplot2_3.3.2       
    ##  [7] dplyr_1.0.2          data.table_1.13.4    BiocParallel_1.24.1 
    ## [10] Trendy_1.10.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] segmented_1.3-4             bitops_1.0-6               
    ##  [3] matrixStats_0.57.0          fs_1.5.0                   
    ##  [5] RColorBrewer_1.1-2          GenomeInfoDb_1.24.2        
    ##  [7] tools_4.0.3                 backports_1.2.0            
    ##  [9] R6_2.5.0                    KernSmooth_2.23-18         
    ## [11] BiocGenerics_0.36.0         colorspace_2.0-0           
    ## [13] GetoptLong_1.0.4            withr_2.3.0                
    ## [15] tidyselect_1.1.0            curl_4.3                   
    ## [17] compiler_4.0.3              Biobase_2.50.0             
    ## [19] DelayedArray_0.14.1         labeling_0.4.2             
    ## [21] caTools_1.18.0              scales_1.1.1               
    ## [23] stringr_1.4.0               digest_0.6.27              
    ## [25] foreign_0.8-80              rmarkdown_2.5              
    ## [27] rio_0.5.16                  XVector_0.28.0             
    ## [29] pkgconfig_2.0.3             htmltools_0.5.0            
    ## [31] dunn.test_1.3.5             fastmap_1.0.1              
    ## [33] readxl_1.3.1                rlang_0.4.9                
    ## [35] GlobalOptions_0.1.2         shiny_1.5.0                
    ## [37] farver_2.0.3                shape_1.4.5                
    ## [39] generics_0.1.0              jsonlite_1.7.1             
    ## [41] gtools_3.8.2                zip_2.1.1                  
    ## [43] car_3.0-10                  RCurl_1.98-1.2             
    ## [45] magrittr_2.0.1              GenomeInfoDbData_1.2.3     
    ## [47] Matrix_1.2-18               Rcpp_1.0.5                 
    ## [49] munsell_0.5.0               S4Vectors_0.28.0           
    ## [51] abind_1.4-5                 lifecycle_0.2.0            
    ## [53] stringi_1.5.3               yaml_2.2.1                 
    ## [55] carData_3.0-4               SummarizedExperiment_1.18.2
    ## [57] zlibbioc_1.34.0             gplots_3.1.1               
    ## [59] shinyFiles_0.9.0            parallel_4.0.3             
    ## [61] promises_1.1.1              forcats_0.5.0              
    ## [63] crayon_1.3.4                lattice_0.20-41            
    ## [65] haven_2.3.1                 splines_4.0.3              
    ## [67] hms_0.5.3                   circlize_0.4.11            
    ## [69] knitr_1.30                  pillar_1.4.7               
    ## [71] GenomicRanges_1.40.0        rjson_0.2.20               
    ## [73] ggsignif_0.6.0              stats4_4.0.3               
    ## [75] glue_1.4.2                  evaluate_0.14              
    ## [77] png_0.1-7                   vctrs_0.3.5                
    ## [79] httpuv_1.5.4                cellranger_1.1.0           
    ## [81] gtable_0.3.0                purrr_0.3.4                
    ## [83] tidyr_1.1.2                 clue_0.3-57                
    ## [85] openxlsx_4.2.3              xfun_0.19                  
    ## [87] mime_0.9                    xtable_1.8-4               
    ## [89] broom_0.7.2                 rstatix_0.6.0              
    ## [91] later_1.1.0.1               tibble_3.0.4               
    ## [93] IRanges_2.24.0              cluster_2.1.0              
    ## [95] ellipsis_0.3.1

Bacher, Rhonda, Ning Leng, Li Fang Chu, Zijian Ni, James A. Thomson, Christina Kendziorski, and Ron Stewart. 2018. “Trendy: Segmented regression analysis of expression dynamics in high-throughput ordered profiling experiments.” *BMC Bioinformatics* 19 (1). BioMed Central Ltd.: 380. doi:[10.1186/s12859-018-2405-x](https://doi.org/10.1186/s12859-018-2405-x).
