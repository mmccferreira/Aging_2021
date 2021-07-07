Functional Analysis
================
Margarida Ferreira
Last updated July 2021

-   [1. Initial setup](#initial-setup)
-   [2. Age-associated modules](#age-associated-modules)
    -   [2.1 DEGs-Hub genes' overlap data input](#degs-hub-genes-overlap-data-input)
    -   [2.2. Functional analysis](#functional-analysis)
        -   [2.2.1 Conversion of gene symbols to EGIDs](#conversion-of-gene-symbols-to-egids)
            -   [Brain](#brain)
            -   [Heart](#heart)
            -   [Muscle](#muscle)
            -   [Liver](#liver)
        -   [2.2.2 Gene Ontology over-representation analysis](#gene-ontology-over-representation-analysis)
            -   [Brain](#brain-1)
            -   [Heart](#heart-1)
            -   [Muscle](#muscle-1)
            -   [Liver](#liver-1)
-   [3. Sex-associated modules](#sex-associated-modules)
    -   [3.1 DEGs-Hub genes' overlap data input](#degs-hub-genes-overlap-data-input-1)
    -   [3.2 Functional analysis](#functional-analysis-1)
        -   [3.2.1 Conversion of gene symbols to EGIDs](#conversion-of-gene-symbols-to-egids-1)
            -   [Brain](#brain-2)
            -   [Muscle](#muscle-2)
            -   [Liver](#liver-2)
        -   [3.2.2 Gene Ontology over-representation analysis](#gene-ontology-over-representation-analysis-1)
            -   [Brain](#brain-3)
            -   [Muscle](#muscle-3)
            -   [Liver](#liver-3)
-   [4. Upset plots](#upset-plots)
    -   [Fig 5 - right](#fig-5---right)

The functional analysis was carried out using the Bioconductor's package clusterProfiler (Yu et al. 2012).

# 1. Initial setup

Loading the required packages:

``` r
library(clusterProfiler)
library(org.Mm.eg.db)
library(UpSetR)
library(ggplot2)
library(grid)
library(ReactomePA)
library(reactome.db)
```

# 2. Age-associated modules

## 2.1 DEGs-Hub genes' overlap data input

``` r
### BRAIN
#tan
upsetList_brain_tan = read.table(file = "upsetList_brain_trendyTan.txt", header = T, sep="\t",check.names = F)
brain_hubTan = na.omit(upsetList_brain_tan$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(brain_hubTan)
```

    ## [1] 34

``` r
brain_moduleTan = na.omit(unique(c(upsetList_brain_tan$`Trendy genes:Module genes`,upsetList_brain_tan$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(brain_moduleTan)
```

    ## [1] 181

``` r
### HEART
#tan
upsetList_heart_tan = read.table(file = "upsetList_heart_trendyTan.txt", header = T, sep="\t",check.names = F)
heart_hubTan = na.omit(upsetList_heart_tan$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(heart_hubTan)
```

    ## [1] 11

``` r
heart_moduleTan = na.omit(unique(c(upsetList_heart_tan$`Trendy genes:Module genes`,upsetList_heart_tan$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(heart_moduleTan)
```

    ## [1] 88

``` r
#blue
upsetList_heart_blue = read.table(file = "upsetList_heart_trendyBlue.txt", header = T, sep="\t",check.names = F)
heart_hubBlue = na.omit(upsetList_heart_blue$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(heart_hubBlue)
```

    ## [1] 94

``` r
heart_moduleBlue = na.omit(unique(c(upsetList_heart_blue$`Trendy genes:Module genes`,upsetList_heart_blue$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(heart_moduleBlue)
```

    ## [1] 374

``` r
### MUSCLE
#magenta
upsetList_muscle_magenta = read.table(file = "upsetList_muscle_trendyMagenta.txt", header = T, sep="\t",check.names = F)
muscle_hubMagenta = na.omit(upsetList_muscle_magenta$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubMagenta)
```

    ## [1] 9

``` r
muscle_moduleMagenta = na.omit(unique(c(upsetList_muscle_magenta$`Trendy genes:Module genes`,upsetList_muscle_magenta$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_moduleMagenta)
```

    ## [1] 70

``` r
#brown
upsetList_muscle_brown = read.table(file = "upsetList_muscle_trendyBrown.txt", header = T, sep="\t",check.names = F)
muscle_hubBrown = na.omit(upsetList_muscle_brown$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubBrown)
```

    ## [1] 31

``` r
muscle_moduleBrown = na.omit(unique(c(upsetList_muscle_brown$`Trendy genes:Module genes`,upsetList_muscle_brown$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_moduleBrown)
```

    ## [1] 146

``` r
### LIVER
#salmon
upsetList_liver_salmon = read.table(file = "upsetList_liver_trendySalmon.txt", header = T, sep="\t",check.names = F)
liver_hubSalmon = na.omit(upsetList_liver_salmon$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubSalmon)
```

    ## [1] 5

``` r
liver_moduleSalmon = na.omit(unique(c(upsetList_liver_salmon$`Trendy genes:Module genes`,upsetList_liver_salmon$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleSalmon)
```

    ## [1] 157

``` r
#darkturquoise
upsetList_liver_darkturquoise = read.table(file = "upsetList_liver_trendyDarkturquoise.txt", header = T, sep="\t",check.names = F)
liver_hubDarkturquoise = na.omit(upsetList_liver_darkturquoise$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubDarkturquoise)
```

    ## [1] 9

``` r
liver_moduleDarkturquoise = na.omit(unique(c(upsetList_liver_darkturquoise$`Trendy genes:Module genes`,upsetList_liver_darkturquoise$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleDarkturquoise)
```

    ## [1] 42

## 2.2. Functional analysis

### 2.2.1 Conversion of gene symbols to EGIDs

This package requires as input a vector of Entrez Gene IDs (EGIDs). For that reason, we used the Bioconductor's R package [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) to perform this conversion.

#### Brain

``` r
### TAN module
length(brain_hubTan) # number of gene symbols
```

    ## [1] 34

``` r
brain_hubTan_conv <- bitr(brain_hubTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(brain_hubTan_conv)
```

    ## [1] 34

``` r
brain_hubTan_egid <- as.vector(brain_hubTan_conv[,2])
length(brain_hubTan_egid)
```

    ## [1] 34

``` r
length(brain_moduleTan) # number of gene symbols
```

    ## [1] 181

``` r
brain_moduleTan_conv <- bitr(brain_moduleTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(brain_moduleTan_conv)
```

    ## [1] 180

``` r
setdiff(brain_moduleTan,brain_moduleTan_conv$SYMBOL)
```

    ## [1] "Gm13905"

``` r
#[1] "Gm13905" has no entrez gene id
brain_moduleTan_egid <- as.vector(brain_moduleTan_conv[,2])
length(brain_moduleTan_egid)
```

    ## [1] 180

#### Heart

``` r
### TAN module
length(heart_hubTan) # number of gene symbols
```

    ## [1] 11

``` r
heart_hubTan_conv <- bitr(heart_hubTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(heart_hubTan_conv)
```

    ## [1] 11

``` r
heart_hubTan_egid <- as.vector(heart_hubTan_conv[,2])
length(heart_hubTan_egid)
```

    ## [1] 11

``` r
length(heart_moduleTan) # number of gene symbols
```

    ## [1] 88

``` r
heart_moduleTan_conv <- bitr(heart_moduleTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(heart_moduleTan_conv)
```

    ## [1] 86

``` r
setdiff(heart_moduleTan,heart_moduleTan_conv$SYMBOL)
```

    ## [1] "Gm44608" "Gm42517"

``` r
#[1] "Gm44608" "Gm42517" (Gm44608" has no entrez gene id)
heart_moduleTan_conv[nrow(heart_moduleTan_conv)+1,] <- c("Gm42517","115490131") #updated recently
heart_moduleTan_egid <- as.vector(heart_moduleTan_conv[,2])
length(heart_moduleTan_egid)
```

    ## [1] 87

``` r
### BLUE module
length(heart_hubBlue) # number of gene symbols
```

    ## [1] 94

``` r
heart_hubBlue_conv <- bitr(heart_hubBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(heart_hubBlue_conv)
```

    ## [1] 94

``` r
heart_hubBlue_egid <- as.vector(heart_hubBlue_conv[,2])
length(heart_hubBlue_egid)
```

    ## [1] 94

``` r
length(heart_moduleBlue) # number of gene symbols
```

    ## [1] 374

``` r
heart_moduleBlue_conv <- bitr(heart_moduleBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(heart_moduleBlue_conv)
```

    ## [1] 374

``` r
heart_moduleBlue_egid <- as.vector(heart_moduleBlue_conv[,2])
length(heart_moduleBlue_egid)
```

    ## [1] 374

#### Muscle

``` r
### MAGENTA module
length(muscle_hubMagenta) # number of gene symbols
```

    ## [1] 9

``` r
muscle_hubMagenta_conv <- bitr(muscle_hubMagenta, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubMagenta_conv)
```

    ## [1] 9

``` r
muscle_hubMagenta_egid <- as.vector(muscle_hubMagenta_conv[,2])
length(muscle_hubMagenta_egid)
```

    ## [1] 9

``` r
length(muscle_moduleMagenta) # number of gene symbols
```

    ## [1] 70

``` r
muscle_moduleMagenta_conv <- bitr(muscle_moduleMagenta, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_moduleMagenta_conv)
```

    ## [1] 68

``` r
setdiff(muscle_moduleMagenta,muscle_moduleMagenta_conv$SYMBOL)
```

    ## [1] "mt-Nd5" "mt-Nd1"

``` r
#[1] "mt-Nd5" "mt-Nd1" 
muscle_moduleMagenta_conv[nrow(muscle_moduleMagenta_conv)+1,] <- c("mt-Nd5","17721")
muscle_moduleMagenta_conv[nrow(muscle_moduleMagenta_conv)+1,] <- c("mt-Nd1","17716")
muscle_moduleMagenta_egid <- as.vector(muscle_moduleMagenta_conv[,2])
length(muscle_moduleMagenta_egid)
```

    ## [1] 70

``` r
### BROWN module
length(muscle_hubBrown) # number of gene symbols
```

    ## [1] 31

``` r
muscle_hubBrown_conv <- bitr(muscle_hubBrown, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubBrown_conv)
```

    ## [1] 31

``` r
muscle_hubBrown_egid <- as.vector(muscle_hubBrown_conv[,2])
length(muscle_hubBrown_egid)
```

    ## [1] 31

``` r
length(muscle_moduleBrown) # number of gene symbols
```

    ## [1] 146

``` r
muscle_moduleBrown_conv <- bitr(muscle_moduleBrown, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_moduleBrown_conv)
```

    ## [1] 146

``` r
muscle_moduleBrown_egid <- as.vector(muscle_moduleBrown_conv[,2])
length(muscle_moduleBrown_egid)
```

    ## [1] 146

#### Liver

``` r
### SALMON module
length(liver_hubSalmon) # number of gene symbols
```

    ## [1] 5

``` r
liver_hubSalmon_conv <- bitr(liver_hubSalmon, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubSalmon_conv)
```

    ## [1] 5

``` r
liver_hubSalmon_egid <- as.vector(liver_hubSalmon_conv[,2])
length(liver_hubSalmon_egid)
```

    ## [1] 5

``` r
length(liver_moduleSalmon) # number of gene symbols
```

    ## [1] 157

``` r
liver_moduleSalmon_conv <- bitr(liver_moduleSalmon, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleSalmon_conv)
```

    ## [1] 156

``` r
setdiff(liver_moduleSalmon,liver_moduleSalmon_conv$SYMBOL)
```

    ## [1] "Gm49410"

``` r
#[1] "Gm49410" has no entrez gene id
liver_moduleSalmon_egid <- as.vector(liver_moduleSalmon_conv[,2])
length(liver_moduleSalmon_egid)
```

    ## [1] 156

``` r
### DARKTURQUOISE module
length(liver_hubDarkturquoise) # number of gene symbols
```

    ## [1] 9

``` r
liver_hubDarkturquoise_conv <- bitr(liver_hubDarkturquoise, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubDarkturquoise_conv)
```

    ## [1] 9

``` r
liver_hubDarkturquoise_egid <- as.vector(liver_hubDarkturquoise_conv[,2])
length(liver_hubDarkturquoise_egid)
```

    ## [1] 9

``` r
length(liver_moduleDarkturquoise) # number of gene symbols
```

    ## [1] 42

``` r
liver_moduleDarkturquoise_conv <- bitr(liver_moduleDarkturquoise, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleDarkturquoise_conv)
```

    ## [1] 42

``` r
liver_moduleDarkturquoise_egid <- as.vector(liver_moduleDarkturquoise_conv[,2])
length(liver_moduleDarkturquoise_egid)
```

    ## [1] 42

### 2.2.2 Gene Ontology over-representation analysis

``` r
#set the universe for the gene ontology analysis as the total number of expressed genes in each tissue
brain_universe <- as.vector(rownames(read.table("vst_brain_counts.txt")))
heart_universe <- as.vector(rownames(read.table("vst_heart_counts.txt")))
muscle_universe <- as.vector(rownames(read.table("vst_muscle_counts.txt")))
liver_universe <- as.vector(rownames(read.table("vst_liver_counts.txt")))

#convert the gene symbols to entrez gene ids
#BRAIN
length(brain_universe) # number of gene symbols
```

    ## [1] 34164

``` r
brain_universe_conv <- bitr(brain_universe, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(brain_universe_conv)
```

    ## [1] 28803

``` r
length(brain_universe_conv$SYMBOL) # number of genes with corresponding entrez gene id
```

    ## [1] 28803

``` r
brain_universe_egid <- as.vector(brain_universe_conv[,2])
length(brain_universe_egid)
```

    ## [1] 28803

``` r
length(brain_universe)-nrow(brain_universe_conv) 
```

    ## [1] 5361

``` r
#[1] 5361 genes lost in the conversion

#HEART
length(heart_universe) # number of gene symbols
```

    ## [1] 28073

``` r
heart_universe_conv <- bitr(heart_universe, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(heart_universe_conv)
```

    ## [1] 24528

``` r
#add previously added egids
heart_universe_conv[nrow(heart_universe_conv)+1,] <- c("Gm42517","115490131")
heart_universe_egid <- as.vector(heart_universe_conv[,2])
length(heart_universe_egid)
```

    ## [1] 24529

``` r
length(heart_universe)-nrow(heart_universe_conv) 
```

    ## [1] 3544

``` r
#[1] 3544 genes lost in the conversion

#MUSCLE
length(muscle_universe) # number of gene symbols
```

    ## [1] 18978

``` r
muscle_universe_conv <- bitr(muscle_universe, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_universe_conv) 
```

    ## [1] 17644

``` r
#add previously added egids
muscle_universe_conv[nrow(muscle_universe_conv)+1,] <- c("mt-Nd5","17721")
muscle_universe_conv[nrow(muscle_universe_conv)+1,] <- c("mt-Nd1","17716")
muscle_universe_egid <- as.vector(muscle_universe_conv[,2])
length(muscle_universe_egid)
```

    ## [1] 17646

``` r
length(muscle_universe)-nrow(muscle_universe_conv) 
```

    ## [1] 1332

``` r
#[1] 1332 genes lost in the conversion

#LIVER
length(liver_universe) # number of gene symbols
```

    ## [1] 20157

``` r
liver_universe_conv <- bitr(liver_universe, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_universe_conv)
```

    ## [1] 18434

``` r
liver_universe_egid <- as.vector(liver_universe_conv[,2])
length(liver_universe_egid)
```

    ## [1] 18434

``` r
length(liver_universe)-nrow(liver_universe_conv) 
```

    ## [1] 1723

``` r
#[1] 1723 genes lost in the conversion
```

Note that the enrichGO function will remove from the universe genes that have no annotation (in this case, genes not present in the GO annotation).

#### Brain

``` r
### TAN module
#hubTan
brain_hubTan_BP <- enrichGO(gene          = brain_hubTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = brain_universe_egid,
                readable = TRUE)
write.table(brain_hubTan_BP@result[brain_hubTan_BP@result$p.adjust <0.05,], "brain_hubTan_BP.txt", quote=FALSE, sep="\t")

brain_hubTan_BP_simp <- clusterProfiler::simplify(brain_hubTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(brain_hubTan_BP_simp@result[brain_hubTan_BP_simp@result$p.adjust <0.05,], "brain_hubTan_BP_simp.txt", quote=FALSE, sep="\t")

brain_hubTan_kegg <- enrichKEGG(gene          = brain_hubTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = brain_universe_egid,
              pvalueCutoff  = 0.05)
write.table(brain_hubTan_kegg@result[brain_hubTan_kegg@result$p.adjust <0.05,], file="brain_hubTan_kegg.txt", quote=FALSE, sep="\t")

brain_hubTan_react <- enrichPathway(gene          = brain_hubTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = brain_universe_egid,
              readable = TRUE)
write.table(brain_hubTan_react@result[brain_hubTan_react@result$p.adjust <0.05,], file="brain_hubTan_react.txt", quote=FALSE, sep="\t")

#moduleTan
brain_moduleTan_BP <- enrichGO(gene          = brain_moduleTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = brain_universe_egid,
                readable = TRUE)
write.table(brain_moduleTan_BP@result[brain_moduleTan_BP@result$p.adjust <0.05,], "brain_moduleTan_BP.txt", quote=FALSE, sep="\t")

brain_moduleTan_BP_simp <- clusterProfiler::simplify(brain_moduleTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(brain_moduleTan_BP_simp@result[brain_moduleTan_BP_simp@result$p.adjust <0.05,], "brain_moduleTan_BP_simp.txt", quote=FALSE, sep="\t")

brain_moduleTan_kegg <- enrichKEGG(gene          = brain_moduleTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = brain_universe_egid,
              pvalueCutoff  = 0.05)
write.table(brain_moduleTan_kegg@result[brain_moduleTan_kegg@result$p.adjust <0.05,], file="brain_moduleTan_kegg.txt", quote=FALSE, sep="\t")

brain_moduleTan_react <- enrichPathway(gene          = brain_moduleTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = brain_universe_egid,
              readable = TRUE)
write.table(brain_moduleTan_react@result[brain_moduleTan_react@result$p.adjust <0.05,], file="brain_moduleTan_react.txt", quote=FALSE, sep="\t")
```

#### Heart

``` r
### TAN module
#hubTan
heart_hubTan_BP <- enrichGO(gene          = heart_hubTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = heart_universe_egid,
                readable = TRUE)
write.table(heart_hubTan_BP@result[heart_hubTan_BP@result$p.adjust <0.05,], "heart_hubTan_BP.txt", quote=FALSE, sep="\t")

heart_hubTan_BP_simp <- clusterProfiler::simplify(heart_hubTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(heart_hubTan_BP_simp@result[heart_hubTan_BP_simp@result$p.adjust <0.05,], "heart_hubTan_BP_simp.txt", quote=FALSE, sep="\t")

heart_hubTan_kegg <- enrichKEGG(gene          = heart_hubTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = heart_universe_egid,
              pvalueCutoff  = 0.05)
write.table(heart_hubTan_kegg@result[heart_hubTan_kegg@result$p.adjust <0.05,], file="heart_hubTan_kegg.txt", quote=FALSE, sep="\t")

heart_hubTan_react <- enrichPathway(gene          = heart_hubTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = heart_universe_egid,
              readable = TRUE)
write.table(heart_hubTan_react@result[heart_hubTan_react@result$p.adjust <0.05,], file="heart_hubTan_react.txt", quote=FALSE, sep="\t")

#moduleTan
heart_moduleTan_BP <- enrichGO(gene          = heart_moduleTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = heart_universe_egid,
                readable = TRUE)
write.table(heart_moduleTan_BP@result[heart_moduleTan_BP@result$p.adjust <0.05,], "heart_moduleTan_BP.txt", quote=FALSE, sep="\t")

heart_moduleTan_BP_simp <- clusterProfiler::simplify(heart_moduleTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(heart_moduleTan_BP_simp@result[heart_moduleTan_BP_simp@result$p.adjust <0.05,], "heart_moduleTan_BP_simp.txt", quote=FALSE, sep="\t")

heart_moduleTan_kegg <- enrichKEGG(gene          = heart_moduleTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = heart_universe_egid,
              pvalueCutoff  = 0.05)
write.table(heart_moduleTan_kegg@result[heart_moduleTan_kegg@result$p.adjust <0.05,], file="heart_moduleTan_kegg.txt", quote=FALSE, sep="\t")

heart_moduleTan_react <- enrichPathway(gene          = heart_moduleTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = heart_universe_egid,
              readable = TRUE)
write.table(heart_moduleTan_react@result[heart_moduleTan_react@result$p.adjust <0.05,], file="heart_moduleTan_react.txt", quote=FALSE, sep="\t")

### BLUE module
#hubBlue
heart_hubBlue_BP <- enrichGO(gene          = heart_hubBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = heart_universe_egid,
                readable = TRUE)
write.table(heart_hubBlue_BP@result[heart_hubBlue_BP@result$p.adjust <0.05,], "heart_hubBlue_BP.txt", quote=FALSE, sep="\t")

heart_hubBlue_BP_simp <- clusterProfiler::simplify(heart_hubBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(heart_hubBlue_BP_simp@result[heart_hubBlue_BP_simp@result$p.adjust <0.05,], "heart_hubBlue_BP_simp.txt", quote=FALSE, sep="\t")

heart_hubBlue_kegg <- enrichKEGG(gene          = heart_hubBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = heart_universe_egid,
              pvalueCutoff  = 0.05)
write.table(heart_hubBlue_kegg@result[heart_hubBlue_kegg@result$p.adjust <0.05,], file="heart_hubBlue_kegg.txt", quote=FALSE, sep="\t")

heart_hubBlue_react <- enrichPathway(gene          = heart_hubBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = heart_universe_egid,
              readable = TRUE)
write.table(heart_hubBlue_react@result[heart_hubBlue_react@result$p.adjust <0.05,], file="heart_hubBlue_react.txt", quote=FALSE, sep="\t")

#moduleHub
heart_moduleBlue_BP <- enrichGO(gene          = heart_moduleBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = heart_universe_egid,
                readable = TRUE)
write.table(heart_moduleBlue_BP@result[heart_moduleBlue_BP@result$p.adjust <0.05,], "heart_moduleBlue_BP.txt", quote=FALSE, sep="\t")

heart_moduleBlue_BP_simp <- clusterProfiler::simplify(heart_moduleBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(heart_moduleBlue_BP_simp@result[heart_moduleBlue_BP_simp@result$p.adjust <0.05,], "heart_moduleBlue_BP_simp.txt", quote=FALSE, sep="\t")

heart_moduleBlue_kegg <- enrichKEGG(gene          = heart_moduleBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = heart_universe_egid,
              pvalueCutoff  = 0.05)
write.table(heart_moduleBlue_kegg@result[heart_moduleBlue_kegg@result$p.adjust <0.05,], file="heart_moduleBlue_kegg.txt", quote=FALSE, sep="\t")

heart_moduleBlue_react <- enrichPathway(gene          = heart_moduleBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = heart_universe_egid,
              readable = TRUE)
write.table(heart_moduleBlue_react@result[heart_moduleBlue_react@result$p.adjust <0.05,], file="heart_moduleBlue_react.txt", quote=FALSE, sep="\t")
```

#### Muscle

``` r
### MAGENTA module
#hubMagenta
muscle_hubMagenta_BP <- enrichGO(gene          = muscle_hubMagenta_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubMagenta_BP@result[muscle_hubMagenta_BP@result$p.adjust <0.05,], "muscle_hubMagenta_BP.txt", quote=FALSE, sep="\t")

muscle_hubMagenta_BP_simp <- clusterProfiler::simplify(muscle_hubMagenta_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubMagenta_BP_simp@result[muscle_hubMagenta_BP_simp@result$p.adjust <0.05,], "muscle_hubMagenta_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubMagenta_kegg <- enrichKEGG(gene          = muscle_hubMagenta_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubMagenta_kegg@result[muscle_hubMagenta_kegg@result$p.adjust <0.05,], file="muscle_hubMagenta_kegg.txt", quote=FALSE, sep="\t")

muscle_hubMagenta_react <- enrichPathway(gene          = muscle_hubMagenta_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubMagenta_react@result[muscle_hubMagenta_react@result$p.adjust <0.05,], file="muscle_hubMagenta_react.txt", quote=FALSE, sep="\t")

#moduleMagenta
muscle_moduleMagenta_BP <- enrichGO(gene          = muscle_moduleMagenta_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_moduleMagenta_BP@result[muscle_moduleMagenta_BP@result$p.adjust <0.05,], "muscle_moduleMagenta_BP.txt", quote=FALSE, sep="\t")

muscle_moduleMagenta_BP_simp <- clusterProfiler::simplify(muscle_moduleMagenta_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_moduleMagenta_BP_simp@result[muscle_moduleMagenta_BP_simp@result$p.adjust <0.05,], "muscle_moduleMagenta_BP_simp.txt", quote=FALSE, sep="\t")

muscle_moduleMagenta_kegg <- enrichKEGG(gene          = muscle_moduleMagenta_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_moduleMagenta_kegg@result[muscle_moduleMagenta_kegg@result$p.adjust <0.05,], file="muscle_moduleMagenta_kegg.txt", quote=FALSE, sep="\t")

muscle_moduleMagenta_react <- enrichPathway(gene          = muscle_moduleMagenta_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_moduleMagenta_react@result[muscle_moduleMagenta_react@result$p.adjust <0.05,], file="muscle_moduleMagenta_react.txt", quote=FALSE, sep="\t")

### BROWN module
#hubBrown
muscle_hubBrown_BP <- enrichGO(gene          = muscle_hubBrown_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubBrown_BP@result[muscle_hubBrown_BP@result$p.adjust <0.05,], "muscle_hubBrown_BP.txt", quote=FALSE, sep="\t")

muscle_hubBrown_BP_simp <- clusterProfiler::simplify(muscle_hubBrown_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubBrown_BP_simp@result[muscle_hubBrown_BP_simp@result$p.adjust <0.05,], "muscle_hubBrown_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubBrown_kegg <- enrichKEGG(gene          = muscle_hubBrown_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubBrown_kegg@result[muscle_hubBrown_kegg@result$p.adjust <0.05,], file="muscle_hubBrown_kegg.txt", quote=FALSE, sep="\t")

muscle_hubBrown_react <- enrichPathway(gene          = muscle_hubBrown_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubBrown_react@result[muscle_hubBrown_react@result$p.adjust <0.05,], file="muscle_hubBrown_react.txt", quote=FALSE, sep="\t")

#moduleHub
muscle_moduleBrown_BP <- enrichGO(gene          = muscle_moduleBrown_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_moduleBrown_BP@result[muscle_moduleBrown_BP@result$p.adjust <0.05,], "muscle_moduleBrown_BP.txt", quote=FALSE, sep="\t")

muscle_moduleBrown_BP_simp <- clusterProfiler::simplify(muscle_moduleBrown_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_moduleBrown_BP_simp@result[muscle_moduleBrown_BP_simp@result$p.adjust <0.05,], "muscle_moduleBrown_BP_simp.txt", quote=FALSE, sep="\t")

muscle_moduleBrown_kegg <- enrichKEGG(gene          = muscle_moduleBrown_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_moduleBrown_kegg@result[muscle_moduleBrown_kegg@result$p.adjust <0.05,], file="muscle_moduleBrown_kegg.txt", quote=FALSE, sep="\t")

muscle_moduleBrown_react <- enrichPathway(gene          = muscle_moduleBrown_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_moduleBrown_react@result[muscle_moduleBrown_react@result$p.adjust <0.05,], file="muscle_moduleBrown_react.txt", quote=FALSE, sep="\t")
```

#### Liver

``` r
### SALMON module
#hubSalmon
liver_hubSalmon_BP <- enrichGO(gene          = liver_hubSalmon_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubSalmon_BP@result[liver_hubSalmon_BP@result$p.adjust <0.05,], "liver_hubSalmon_BP.txt", quote=FALSE, sep="\t")

liver_hubSalmon_BP_simp <- clusterProfiler::simplify(liver_hubSalmon_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubSalmon_BP_simp@result[liver_hubSalmon_BP_simp@result$p.adjust <0.05,], "liver_hubSalmon_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubSalmon_kegg <- enrichKEGG(gene          = liver_hubSalmon_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubSalmon_kegg@result[liver_hubSalmon_kegg@result$p.adjust <0.05,], file="liver_hubSalmon_kegg.txt", quote=FALSE, sep="\t")

liver_hubSalmon_react <- enrichPathway(gene          = liver_hubSalmon_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubSalmon_react@result[liver_hubSalmon_react@result$p.adjust <0.05,], file="liver_hubSalmon_react.txt", quote=FALSE, sep="\t")

#moduleSalmon
liver_moduleSalmon_BP <- enrichGO(gene          = liver_moduleSalmon_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleSalmon_BP@result[liver_moduleSalmon_BP@result$p.adjust <0.05,], "liver_moduleSalmon_BP.txt", quote=FALSE, sep="\t")

liver_moduleSalmon_BP_simp <- clusterProfiler::simplify(liver_moduleSalmon_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleSalmon_BP_simp@result[liver_moduleSalmon_BP_simp@result$p.adjust <0.05,], "liver_moduleSalmon_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleSalmon_kegg <- enrichKEGG(gene          = liver_moduleSalmon_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleSalmon_kegg@result[liver_moduleSalmon_kegg@result$p.adjust <0.05,], file="liver_moduleSalmon_kegg.txt", quote=FALSE, sep="\t")

liver_moduleSalmon_react <- enrichPathway(gene          = liver_moduleSalmon_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleSalmon_react@result[liver_moduleSalmon_react@result$p.adjust <0.05,], file="liver_moduleSalmon_react.txt", quote=FALSE, sep="\t")

### DARKTURQUOISE module
#hubDarkturquoise
liver_hubDarkturquoise_BP <- enrichGO(gene          = liver_hubDarkturquoise_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubDarkturquoise_BP@result[liver_hubDarkturquoise_BP@result$p.adjust <0.05,], "liver_hubDarkturquoise_BP.txt", quote=FALSE, sep="\t")

liver_hubDarkturquoise_BP_simp <- clusterProfiler::simplify(liver_hubDarkturquoise_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubDarkturquoise_BP_simp@result[liver_hubDarkturquoise_BP_simp@result$p.adjust <0.05,], "liver_hubDarkturquoise_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubDarkturquoise_kegg <- enrichKEGG(gene          = liver_hubDarkturquoise_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubDarkturquoise_kegg@result[liver_hubDarkturquoise_kegg@result$p.adjust <0.05,], file="liver_hubDarkturquoise_kegg.txt", quote=FALSE, sep="\t")

liver_hubDarkturquoise_react <- enrichPathway(gene          = liver_hubDarkturquoise_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubDarkturquoise_react@result[liver_hubDarkturquoise_react@result$p.adjust <0.05,], file="liver_hubDarkturquoise_react.txt", quote=FALSE, sep="\t")

#moduleHub
liver_moduleDarkturquoise_BP <- enrichGO(gene          = liver_moduleDarkturquoise_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleDarkturquoise_BP@result[liver_moduleDarkturquoise_BP@result$p.adjust <0.05,], "liver_moduleDarkturquoise_BP.txt", quote=FALSE, sep="\t")

liver_moduleDarkturquoise_BP_simp <- clusterProfiler::simplify(liver_moduleDarkturquoise_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleDarkturquoise_BP_simp@result[liver_moduleDarkturquoise_BP_simp@result$p.adjust <0.05,], "liver_moduleDarkturquoise_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleDarkturquoise_kegg <- enrichKEGG(gene          = liver_moduleDarkturquoise_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleDarkturquoise_kegg@result[liver_moduleDarkturquoise_kegg@result$p.adjust <0.05,], file="liver_moduleDarkturquoise_kegg.txt", quote=FALSE, sep="\t")

liver_moduleDarkturquoise_react <- enrichPathway(gene          = liver_moduleDarkturquoise_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleDarkturquoise_react@result[liver_moduleDarkturquoise_react@result$p.adjust <0.05,], file="liver_moduleDarkturquoise_react.txt", quote=FALSE, sep="\t")
```

# 3. Sex-associated modules

## 3.1 DEGs-Hub genes' overlap data input

``` r
### BRAIN
#Grey60
upsetList_brain_grey60 = read.table(file = "upsetList_brain_trendyGrey60.txt", header = T, sep="\t",check.names = F)
brain_hubGrey60 = na.omit(upsetList_brain_grey60$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(brain_hubGrey60)
```

    ## [1] 2

``` r
brain_moduleGrey60 = na.omit(unique(c(upsetList_brain_grey60$`Trendy genes:Module genes`,upsetList_brain_grey60$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(brain_moduleGrey60)
```

    ## [1] 11

``` r
### MUSCLE
#red
upsetList_muscle_red = read.table(file = "upsetList_muscle_trendyRed.txt", header = T, sep="\t",check.names = F)
muscle_hubRed = na.omit(upsetList_muscle_red$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubRed)
```

    ## [1] 21

``` r
muscle_moduleRed = na.omit(unique(c(upsetList_muscle_red$`Trendy genes:Module genes`,upsetList_muscle_red$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_moduleRed)
```

    ## [1] 135

``` r
#purple
upsetList_muscle_purple = read.table(file = "upsetList_muscle_trendyPurple.txt", header = T, sep="\t",check.names = F)
muscle_hubPurple = na.omit(upsetList_muscle_purple$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubPurple)
```

    ## [1] 9

``` r
muscle_modulePurple = na.omit(unique(c(upsetList_muscle_purple$`Trendy genes:Module genes`,upsetList_muscle_purple$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_modulePurple)
```

    ## [1] 49

``` r
#greenyellow
upsetList_muscle_greenyellow = read.table(file = "upsetList_muscle_trendyGreenyellow.txt", header = T, sep="\t",check.names = F)
muscle_hubGreenyellow = na.omit(upsetList_muscle_greenyellow$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubGreenyellow)
```

    ## [1] 1

``` r
muscle_moduleGreenyellow = na.omit(unique(c(upsetList_muscle_greenyellow$`Trendy genes:Module genes`,upsetList_muscle_greenyellow$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_moduleGreenyellow)
```

    ## [1] 41

``` r
#blue
upsetList_muscle_blue = read.table(file = "upsetList_muscle_trendyBlue.txt", header = T, sep="\t",check.names = F)
muscle_hubBlue = na.omit(upsetList_muscle_blue$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(muscle_hubBlue)
```

    ## [1] 21

``` r
muscle_moduleBlue = na.omit(unique(c(upsetList_muscle_blue$`Trendy genes:Module genes`,upsetList_muscle_blue$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(muscle_moduleBlue)
```

    ## [1] 58

``` r
### LIVER
#tan
upsetList_liver_tan = read.table(file = "upsetList_liver_trendyTan.txt", header = T, sep="\t",check.names = F)
liver_hubTan = na.omit(upsetList_liver_tan$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubTan)
```

    ## [1] 13

``` r
liver_moduleTan = na.omit(unique(c(upsetList_liver_tan$`Trendy genes:Module genes`,upsetList_liver_tan$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleTan)
```

    ## [1] 146

``` r
#red
upsetList_liver_red = read.table(file = "upsetList_liver_trendyRed.txt", header = T, sep="\t",check.names = F)
liver_hubRed = na.omit(upsetList_liver_red$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubRed)
```

    ## [1] 18

``` r
liver_moduleRed = na.omit(unique(c(upsetList_liver_red$`Trendy genes:Module genes`,upsetList_liver_red$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleRed)
```

    ## [1] 183

``` r
#darkolivegreen
upsetList_liver_darkolivegreen = read.table(file = "upsetList_liver_trendyDarkolivegreen.txt", header = T, sep="\t",check.names = F)
liver_hubDarkolivegreen = na.omit(upsetList_liver_darkolivegreen$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubDarkolivegreen)
```

    ## [1] 4

``` r
liver_moduleDarkolivegreen = na.omit(unique(c(upsetList_liver_darkolivegreen$`Trendy genes:Module genes`,upsetList_liver_darkolivegreen$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleDarkolivegreen)
```

    ## [1] 17

``` r
#darkgrey
upsetList_liver_darkgrey = read.table(file = "upsetList_liver_trendyDarkgrey.txt", header = T, sep="\t",check.names = F)
liver_hubDarkgrey = na.omit(upsetList_liver_darkgrey$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubDarkgrey)
```

    ## [1] 20

``` r
liver_moduleDarkgrey = na.omit(unique(c(upsetList_liver_darkgrey$`Trendy genes:Module genes`,upsetList_liver_darkgrey$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleDarkgrey)
```

    ## [1] 76

``` r
#cyan
upsetList_liver_cyan = read.table(file = "upsetList_liver_trendyCyan.txt", header = T, sep="\t",check.names = F)
liver_hubCyan = na.omit(upsetList_liver_cyan$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubCyan)
```

    ## [1] 12

``` r
liver_moduleCyan = na.omit(unique(c(upsetList_liver_cyan$`Trendy genes:Module genes`,upsetList_liver_cyan$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleCyan)
```

    ## [1] 138

``` r
#brown
upsetList_liver_brown = read.table(file = "upsetList_liver_trendyBrown.txt", header = T, sep="\t",check.names = F)
liver_hubBrown = na.omit(upsetList_liver_brown$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubBrown)
```

    ## [1] 40

``` r
liver_moduleBrown = na.omit(unique(c(upsetList_liver_brown$`Trendy genes:Module genes`,upsetList_liver_brown$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleBrown)
```

    ## [1] 271

``` r
#blue
upsetList_liver_blue = read.table(file = "upsetList_liver_trendyBlue.txt", header = T, sep="\t",check.names = F)
liver_hubBlue = na.omit(upsetList_liver_blue$`Trendy genes:Module genes:Hub genes`) #TRENDY-Hub genes' overlap
length(liver_hubBlue)
```

    ## [1] 92

``` r
liver_moduleBlue = na.omit(unique(c(upsetList_liver_blue$`Trendy genes:Module genes`,upsetList_liver_blue$`Trendy genes:Module genes:Hub genes`))) #TRENDY-Module genes' overlap
length(liver_moduleBlue)
```

    ## [1] 438

## 3.2 Functional analysis

### 3.2.1 Conversion of gene symbols to EGIDs

This package requires as input a vector of Entrez Gene IDs (EGIDs). For that reason, we used the Bioconductor's R package [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) to perform this conversion.

#### Brain

``` r
### GREY60 module
length(brain_hubGrey60) # number of gene symbols
```

    ## [1] 2

``` r
brain_hubGrey60_conv <- bitr(brain_hubGrey60, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(brain_hubGrey60_conv)
```

    ## [1] 2

``` r
brain_hubGrey60_egid <- as.vector(brain_hubGrey60_conv[,2])
length(brain_hubGrey60_egid)
```

    ## [1] 2

``` r
length(brain_moduleGrey60) # number of gene symbols
```

    ## [1] 11

``` r
brain_moduleGrey60_conv <- bitr(brain_moduleGrey60, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(brain_moduleGrey60_conv)
```

    ## [1] 7

``` r
setdiff(brain_moduleGrey60,brain_moduleGrey60_conv$SYMBOL)
```

    ## [1] "Gm28476" "Gm29392" "Gm28673" "Gm29158"

``` r
#[1] "Gm28476" "Gm29392" "Gm28673" "Gm29158" have no entrez gene ids
brain_moduleGrey60_egid <- as.vector(brain_moduleGrey60_conv[,2])
length(brain_moduleGrey60_egid)
```

    ## [1] 7

#### Muscle

``` r
### RED module
length(muscle_hubRed) # number of gene symbols
```

    ## [1] 21

``` r
muscle_hubRed_conv <- bitr(muscle_hubRed, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubRed_conv)
```

    ## [1] 21

``` r
muscle_hubRed_egid <- as.vector(muscle_hubRed_conv[,2])
length(muscle_hubRed_egid)
```

    ## [1] 21

``` r
length(muscle_moduleRed) # number of gene symbols
```

    ## [1] 135

``` r
muscle_moduleRed_conv <- bitr(muscle_moduleRed, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_moduleRed_conv)
```

    ## [1] 132

``` r
setdiff(muscle_moduleRed,muscle_moduleRed_conv$SYMBOL)
```

    ## [1] "Gm38158" "Gm7265"  "Platr31"

``` r
#[1] "Gm38158" "Gm7265"  "Platr31" ("Platr31" has no entrez gene ids)
muscle_moduleRed_conv[nrow(muscle_moduleRed_conv)+1,] <- c("Gm38158","118567646")
muscle_moduleRed_conv[nrow(muscle_moduleRed_conv)+1,] <- c("Nedd4","17999") #previously known as Gm7265; updated recently
muscle_moduleRed_egid <- as.vector(muscle_moduleRed_conv[,2])
length(muscle_moduleRed_egid) 
```

    ## [1] 134

``` r
### PURPLE module
length(muscle_hubPurple) # number of gene symbols
```

    ## [1] 9

``` r
muscle_hubPurple_conv <- bitr(muscle_hubPurple, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubPurple_conv)
```

    ## [1] 9

``` r
muscle_hubPurple_egid <- as.vector(muscle_hubPurple_conv[,2])
length(muscle_hubPurple_egid)
```

    ## [1] 9

``` r
length(muscle_modulePurple) # number of gene symbols
```

    ## [1] 49

``` r
muscle_modulePurple_conv <- bitr(muscle_modulePurple, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_modulePurple_conv)
```

    ## [1] 49

``` r
muscle_modulePurple_egid <- as.vector(muscle_modulePurple_conv[,2])
length(muscle_modulePurple_egid)
```

    ## [1] 49

``` r
### GREENYELLOW module
length(muscle_hubGreenyellow) # number of gene symbols
```

    ## [1] 1

``` r
muscle_hubGreenyellow_conv <- bitr(muscle_hubGreenyellow, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubGreenyellow_conv)
```

    ## [1] 1

``` r
muscle_hubGreenyellow_egid <- as.vector(muscle_hubGreenyellow_conv[,2])
length(muscle_hubGreenyellow_egid)
```

    ## [1] 1

``` r
length(muscle_moduleGreenyellow) # number of gene symbols
```

    ## [1] 41

``` r
muscle_moduleGreenyellow_conv <- bitr(muscle_moduleGreenyellow, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_moduleGreenyellow_conv)
```

    ## [1] 41

``` r
muscle_moduleGreenyellow_egid <- as.vector(muscle_moduleGreenyellow_conv[,2])
length(muscle_moduleGreenyellow_egid)
```

    ## [1] 41

``` r
### BLUE module
length(muscle_hubBlue) # number of gene symbols
```

    ## [1] 21

``` r
muscle_hubBlue_conv <- bitr(muscle_hubBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_hubBlue_conv)
```

    ## [1] 20

``` r
setdiff(muscle_hubBlue,muscle_hubBlue_conv$SYMBOL)
```

    ## [1] "mt-Nd2"

``` r
#[1] "mt-Nd2"
muscle_hubBlue_conv[nrow(muscle_hubBlue_conv)+1,] <- c("mt-Nd2","17717")
muscle_hubBlue_egid <- as.vector(muscle_hubBlue_conv[,2])
length(muscle_hubBlue_egid)
```

    ## [1] 21

``` r
length(muscle_moduleBlue) # number of gene symbols
```

    ## [1] 58

``` r
muscle_moduleBlue_conv <- bitr(muscle_moduleBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(muscle_moduleBlue_conv)
```

    ## [1] 56

``` r
setdiff(muscle_moduleBlue,muscle_moduleBlue_conv$SYMBOL)
```

    ## [1] "CT010467.1" "mt-Nd2"

``` r
#[1] "CT010467.1"  "mt-Nd2" ("CT010467.1" has no entrez gene ids)
muscle_moduleBlue_conv[nrow(muscle_moduleBlue_conv)+1,] <- c("mt-Nd2","17717")
muscle_moduleBlue_egid <- as.vector(muscle_moduleBlue_conv[,2])
length(muscle_moduleBlue_egid)
```

    ## [1] 57

#### Liver

``` r
### TAN module
length(liver_hubTan) # number of gene symbols
```

    ## [1] 13

``` r
liver_hubTan_conv <- bitr(liver_hubTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubTan_conv)
```

    ## [1] 13

``` r
liver_hubTan_egid <- as.vector(liver_hubTan_conv[,2])
length(liver_hubTan_egid)
```

    ## [1] 13

``` r
length(liver_moduleTan) # number of gene symbols
```

    ## [1] 146

``` r
liver_moduleTan_conv <- bitr(liver_moduleTan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleTan_conv)
```

    ## [1] 146

``` r
liver_moduleTan_egid <- as.vector(liver_moduleTan_conv[,2])
length(liver_moduleTan_egid)
```

    ## [1] 146

``` r
### RED module
length(liver_hubRed) # number of gene symbols
```

    ## [1] 18

``` r
liver_hubRed_conv <- bitr(liver_hubRed, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubRed_conv)
```

    ## [1] 18

``` r
liver_hubRed_egid <- as.vector(liver_hubRed_conv[,2])
length(liver_hubRed_egid)
```

    ## [1] 18

``` r
length(liver_moduleRed) # number of gene symbols
```

    ## [1] 183

``` r
liver_moduleRed_conv <- bitr(liver_moduleRed, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleRed_conv)
```

    ## [1] 180

``` r
setdiff(liver_moduleRed,liver_moduleRed_conv$SYMBOL)
```

    ## [1] "Gm48660" "Gm43449" "Gm48488"

``` r
#[1] "Gm48660" "Gm43449" "Gm48488" have no entrez gene ids
liver_moduleRed_egid <- as.vector(liver_moduleRed_conv[,2])
length(liver_moduleRed_egid)
```

    ## [1] 180

``` r
### DARKOLIVEGREEN module
length(liver_hubDarkolivegreen) # number of gene symbols
```

    ## [1] 4

``` r
liver_hubDarkolivegreen_conv <- bitr(liver_hubDarkolivegreen, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubDarkolivegreen_conv)
```

    ## [1] 4

``` r
liver_hubDarkolivegreen_egid <- as.vector(liver_hubDarkolivegreen_conv[,2])
length(liver_hubDarkolivegreen_egid)
```

    ## [1] 4

``` r
length(liver_moduleDarkolivegreen) # number of gene symbols
```

    ## [1] 17

``` r
liver_moduleDarkolivegreen_conv <- bitr(liver_moduleDarkolivegreen, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleDarkolivegreen_conv)
```

    ## [1] 17

``` r
liver_moduleDarkolivegreen_egid <- as.vector(liver_moduleDarkolivegreen_conv[,2])
length(liver_moduleDarkolivegreen_egid)
```

    ## [1] 17

``` r
### DARKGREY module
length(liver_hubDarkgrey) # number of gene symbols
```

    ## [1] 20

``` r
liver_hubDarkgrey_conv <- bitr(liver_hubDarkgrey, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubDarkgrey_conv)
```

    ## [1] 20

``` r
liver_hubDarkgrey_egid <- as.vector(liver_hubDarkgrey_conv[,2])
length(liver_hubDarkgrey_egid)
```

    ## [1] 20

``` r
length(liver_moduleDarkgrey) # number of gene symbols
```

    ## [1] 76

``` r
liver_moduleDarkgrey_conv <- bitr(liver_moduleDarkgrey, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleDarkgrey_conv)
```

    ## [1] 75

``` r
setdiff(liver_moduleDarkgrey,liver_moduleDarkgrey_conv$SYMBOL)
```

    ## [1] "Gm48120"

``` r
#[1] "Gm48120" has no entrez gene ids
liver_moduleDarkgrey_egid <- as.vector(liver_moduleDarkgrey_conv[,2])
length(liver_moduleDarkgrey_egid)
```

    ## [1] 75

``` r
### CYAN module
length(liver_hubCyan) # number of gene symbols
```

    ## [1] 12

``` r
liver_hubCyan_conv <- bitr(liver_hubCyan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubCyan_conv)
```

    ## [1] 12

``` r
liver_hubCyan_egid <- as.vector(liver_hubCyan_conv[,2])
length(liver_hubCyan_egid)
```

    ## [1] 12

``` r
length(liver_moduleCyan) # number of gene symbols
```

    ## [1] 138

``` r
liver_moduleCyan_conv <- bitr(liver_moduleCyan, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleCyan_conv)
```

    ## [1] 138

``` r
liver_moduleCyan_egid <- as.vector(liver_moduleCyan_conv[,2])
length(liver_moduleCyan_egid)
```

    ## [1] 138

``` r
### BROWN module
length(liver_hubBrown) # number of gene symbols
```

    ## [1] 40

``` r
liver_hubBrown_conv <- bitr(liver_hubBrown, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubBrown_conv)
```

    ## [1] 40

``` r
liver_hubBrown_egid <- as.vector(liver_hubBrown_conv[,2])
length(liver_hubBrown_egid)
```

    ## [1] 40

``` r
length(liver_moduleBrown) # number of gene symbols
```

    ## [1] 271

``` r
liver_moduleBrown_conv <- bitr(liver_moduleBrown, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleBrown_conv)
```

    ## [1] 265

``` r
setdiff(liver_moduleBrown,liver_moduleBrown_conv$SYMBOL)
```

    ## [1] "Gm15261"   "Gm37969"   "Gm37332"   "Gm47376"   "Gm22499"   "Pard3bos3"

``` r
#[1] "Gm15261"   "Gm37969"   "Gm37332"   "Gm47376"   "Gm22499"   "Pard3bos3" ("Gm15261"   "Gm37969"   "Gm37332"   "Gm47376"  "Pard3bos3" have no entrez gene ids)
liver_moduleBrown_conv[nrow(liver_moduleBrown_conv)+1,] <- c("Gm22499","115485942") #previously known as Gm15261; updated recently
liver_moduleBrown_egid <- as.vector(liver_moduleBrown_conv[,2])
length(liver_moduleBrown_egid)
```

    ## [1] 266

``` r
### BLUE module
length(liver_hubBlue) # number of gene symbols
```

    ## [1] 92

``` r
liver_hubBlue_conv <- bitr(liver_hubBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_hubBlue_conv)
```

    ## [1] 91

``` r
setdiff(liver_hubBlue,liver_hubBlue_conv$SYMBOL)
```

    ## [1] "Gm49338"

``` r
#[1] "Gm49338" has no entry gene id
liver_hubBlue_egid <- as.vector(liver_hubBlue_conv[,2])
length(liver_hubBlue_egid)
```

    ## [1] 91

``` r
length(liver_moduleBlue) # number of gene symbols
```

    ## [1] 438

``` r
liver_moduleBlue_conv <- bitr(liver_moduleBlue, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Mm.eg.db)
nrow(liver_moduleBlue_conv)
```

    ## [1] 429

``` r
setdiff(liver_moduleBlue,liver_moduleBlue_conv$SYMBOL)
```

    ## [1] "Gm26967"  "Mup-ps20" "Gm48855"  "Gm49405"  "Gm47465"  "Gm48199"  "Gm45724" 
    ## [8] "Gm38266"  "Gm49338"

``` r
#[1] "Gm26967"  "Mup-ps20" "Gm48855"  "Gm49405"  "Gm47465"  "Gm48199" 
#[7] "Gm45724"  "Gm38266"  "Gm49338"  "Gm38266"  "Gm49338" ("Gm26967"  "Mup-ps20" "Gm48855"  "Gm49405"  "Gm45724"  "Gm38266"  "Gm49338"  "Gm38266"  "Gm49338" have no entrez gene id)
liver_moduleBlue_conv[nrow(liver_moduleBlue_conv)+1,] <- c("Gm47465","118567729") #updated recently
liver_moduleBlue_conv[nrow(liver_moduleBlue_conv)+1,] <- c("Gm48199","118568049")
liver_moduleBlue_egid <- as.vector(liver_moduleBlue_conv[,2])
length(liver_moduleBlue_egid)
```

    ## [1] 431

### 3.2.2 Gene Ontology over-representation analysis

#### Brain

``` r
### grey60 module
#hubGrey60
brain_hubGrey60_BP <- enrichGO(gene          = brain_hubGrey60_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = brain_universe_egid,
                readable = TRUE)
# write.table(brain_hubGrey60_BP@result[brain_hubGrey60_BP@result$p.adjust <0.05,], "brain_hubGrey60_BP.txt", quote=FALSE, sep="\t")

# brain_hubGrey60_BP_simp <- clusterProfiler::simplify(brain_hubGrey60_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
# write.table(brain_hubGrey60_BP_simp@result[brain_hubGrey60_BP_simp@result$p.adjust <0.05,], "brain_hubGrey60_BP_simp.txt", quote=FALSE, sep="\t")

brain_hubGrey60_kegg <- enrichKEGG(gene          = brain_hubGrey60_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = brain_universe_egid,
              pvalueCutoff  = 0.05)
# write.table(brain_hubGrey60_kegg@result[brain_hubGrey60_kegg@result$p.adjust <0.05,], file="brain_hubGrey60_kegg.txt", quote=FALSE, sep="\t")

brain_hubGrey60_react <- enrichPathway(gene          = brain_hubGrey60_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = brain_universe_egid,
              readable = TRUE)
# write.table(brain_hubGrey60_react@result[brain_hubGrey60_react@result$p.adjust <0.05,], file="brain_hubGrey60_react.txt", quote=FALSE, sep="\t")

#moduleGrey60
brain_moduleGrey60_BP <- enrichGO(gene          = brain_moduleGrey60_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = brain_universe_egid,
                readable = TRUE)
# write.table(brain_moduleGrey60_BP@result[brain_moduleGrey60_BP@result$p.adjust <0.05,], "brain_moduleGrey60_BP.txt", quote=FALSE, sep="\t")

# brain_moduleGrey60_BP_simp <- clusterProfiler::simplify(brain_moduleGrey60_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
# write.table(brain_moduleGrey60_BP_simp@result[brain_moduleGrey60_BP_simp@result$p.adjust <0.05,], "brain_moduleGrey60_BP_simp.txt", quote=FALSE, sep="\t")

brain_moduleGrey60_kegg <- enrichKEGG(gene          = brain_moduleGrey60_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = brain_universe_egid,
              pvalueCutoff  = 0.05)
# write.table(brain_moduleGrey60_kegg@result[brain_moduleGrey60_kegg@result$p.adjust <0.05,], file="brain_moduleGrey60_kegg.txt", quote=FALSE, sep="\t")

brain_moduleGrey60_react <- enrichPathway(gene          = brain_moduleGrey60_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = brain_universe_egid,
              readable = TRUE)
# write.table(brain_moduleGrey60_react@result[brain_moduleGrey60_react@result$p.adjust <0.05,], file="brain_moduleGrey60_react.txt", quote=FALSE, sep="\t")
```

#### Muscle

``` r
### RED module
#hubRed
muscle_hubRed_BP <- enrichGO(gene          = muscle_hubRed_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubRed_BP@result[muscle_hubRed_BP@result$p.adjust <0.05,], "muscle_hubRed_BP.txt", quote=FALSE, sep="\t")

muscle_hubRed_BP_simp <- clusterProfiler::simplify(muscle_hubRed_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubRed_BP_simp@result[muscle_hubRed_BP_simp@result$p.adjust <0.05,], "muscle_hubRed_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubRed_kegg <- enrichKEGG(gene          = muscle_hubRed_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubRed_kegg@result[muscle_hubRed_kegg@result$p.adjust <0.05,], file="muscle_hubRed_kegg.txt", quote=FALSE, sep="\t")

muscle_hubRed_react <- enrichPathway(gene          = muscle_hubRed_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubRed_react@result[muscle_hubRed_react@result$p.adjust <0.05,], file="muscle_hubRed_react.txt", quote=FALSE, sep="\t")

#moduleRed
muscle_moduleRed_BP <- enrichGO(gene          = muscle_moduleRed_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_moduleRed_BP@result[muscle_moduleRed_BP@result$p.adjust <0.05,], "muscle_moduleRed_BP.txt", quote=FALSE, sep="\t")

muscle_moduleRed_BP_simp <- clusterProfiler::simplify(muscle_moduleRed_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_moduleRed_BP_simp@result[muscle_moduleRed_BP_simp@result$p.adjust <0.05,], "muscle_moduleRed_BP_simp.txt", quote=FALSE, sep="\t")

muscle_moduleRed_kegg <- enrichKEGG(gene          = muscle_moduleRed_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_moduleRed_kegg@result[muscle_moduleRed_kegg@result$p.adjust <0.05,], file="muscle_moduleRed_kegg.txt", quote=FALSE, sep="\t")

muscle_moduleRed_react <- enrichPathway(gene          = muscle_moduleRed_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_moduleRed_react@result[muscle_moduleRed_react@result$p.adjust <0.05,], file="muscle_moduleRed_react.txt", quote=FALSE, sep="\t")

### PURPLE module
#hubPurple
muscle_hubPurple_BP <- enrichGO(gene          = muscle_hubPurple_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubPurple_BP@result[muscle_hubPurple_BP@result$p.adjust <0.05,], "muscle_hubPurple_BP.txt", quote=FALSE, sep="\t")

muscle_hubPurple_BP_simp <- clusterProfiler::simplify(muscle_hubPurple_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubPurple_BP_simp@result[muscle_hubPurple_BP_simp@result$p.adjust <0.05,], "muscle_hubPurple_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubPurple_kegg <- enrichKEGG(gene          = muscle_hubPurple_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubPurple_kegg@result[muscle_hubPurple_kegg@result$p.adjust <0.05,], file="muscle_hubPurple_kegg.txt", quote=FALSE, sep="\t")

muscle_hubPurple_react <- enrichPathway(gene          = muscle_hubPurple_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubPurple_react@result[muscle_hubPurple_react@result$p.adjust <0.05,], file="muscle_hubPurple_react.txt", quote=FALSE, sep="\t")

#moduleHub
muscle_modulePurple_BP <- enrichGO(gene          = muscle_modulePurple_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_modulePurple_BP@result[muscle_modulePurple_BP@result$p.adjust <0.05,], "muscle_modulePurple_BP.txt", quote=FALSE, sep="\t")

muscle_modulePurple_BP_simp <- clusterProfiler::simplify(muscle_modulePurple_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_modulePurple_BP_simp@result[muscle_modulePurple_BP_simp@result$p.adjust <0.05,], "muscle_modulePurple_BP_simp.txt", quote=FALSE, sep="\t")

muscle_modulePurple_kegg <- enrichKEGG(gene          = muscle_modulePurple_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_modulePurple_kegg@result[muscle_modulePurple_kegg@result$p.adjust <0.05,], file="muscle_modulePurple_kegg.txt", quote=FALSE, sep="\t")

muscle_modulePurple_react <- enrichPathway(gene          = muscle_modulePurple_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_modulePurple_react@result[muscle_modulePurple_react@result$p.adjust <0.05,], file="muscle_modulePurple_react.txt", quote=FALSE, sep="\t")

### GREENYELLOW module
#hubGreenyellow
muscle_hubGreenyellow_BP <- enrichGO(gene          = muscle_hubGreenyellow_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubGreenyellow_BP@result[muscle_hubGreenyellow_BP@result$p.adjust <0.05,], "muscle_hubGreenyellow_BP.txt", quote=FALSE, sep="\t")

muscle_hubGreenyellow_BP_simp <- clusterProfiler::simplify(muscle_hubGreenyellow_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubGreenyellow_BP_simp@result[muscle_hubGreenyellow_BP_simp@result$p.adjust <0.05,], "muscle_hubGreenyellow_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubGreenyellow_kegg <- enrichKEGG(gene          = muscle_hubGreenyellow_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubGreenyellow_kegg@result[muscle_hubGreenyellow_kegg@result$p.adjust <0.05,], file="muscle_hubGreenyellow_kegg.txt", quote=FALSE, sep="\t")

muscle_hubGreenyellow_react <- enrichPathway(gene          = muscle_hubGreenyellow_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubGreenyellow_react@result[muscle_hubGreenyellow_react@result$p.adjust <0.05,], file="muscle_hubGreenyellow_react.txt", quote=FALSE, sep="\t")

#moduleHub
muscle_moduleGreenyellow_BP <- enrichGO(gene          = muscle_moduleGreenyellow_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_moduleGreenyellow_BP@result[muscle_moduleGreenyellow_BP@result$p.adjust <0.05,], "muscle_moduleGreenyellow_BP.txt", quote=FALSE, sep="\t")

muscle_moduleGreenyellow_BP_simp <- clusterProfiler::simplify(muscle_moduleGreenyellow_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_moduleGreenyellow_BP_simp@result[muscle_moduleGreenyellow_BP_simp@result$p.adjust <0.05,], "muscle_moduleGreenyellow_BP_simp.txt", quote=FALSE, sep="\t")

muscle_moduleGreenyellow_kegg <- enrichKEGG(gene          = muscle_moduleGreenyellow_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_moduleGreenyellow_kegg@result[muscle_moduleGreenyellow_kegg@result$p.adjust <0.05,], file="muscle_moduleGreenyellow_kegg.txt", quote=FALSE, sep="\t")

muscle_moduleGreenyellow_react <- enrichPathway(gene          = muscle_moduleGreenyellow_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_moduleGreenyellow_react@result[muscle_moduleGreenyellow_react@result$p.adjust <0.05,], file="muscle_moduleGreenyellow_react.txt", quote=FALSE, sep="\t")

### BLUE module
#hubBlue
muscle_hubBlue_BP <- enrichGO(gene          = muscle_hubBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_hubBlue_BP@result[muscle_hubBlue_BP@result$p.adjust <0.05,], "muscle_hubBlue_BP.txt", quote=FALSE, sep="\t")

muscle_hubBlue_BP_simp <- clusterProfiler::simplify(muscle_hubBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_hubBlue_BP_simp@result[muscle_hubBlue_BP_simp@result$p.adjust <0.05,], "muscle_hubBlue_BP_simp.txt", quote=FALSE, sep="\t")

muscle_hubBlue_kegg <- enrichKEGG(gene          = muscle_hubBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_hubBlue_kegg@result[muscle_hubBlue_kegg@result$p.adjust <0.05,], file="muscle_hubBlue_kegg.txt", quote=FALSE, sep="\t")

muscle_hubBlue_react <- enrichPathway(gene          = muscle_hubBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_hubBlue_react@result[muscle_hubBlue_react@result$p.adjust <0.05,], file="muscle_hubBlue_react.txt", quote=FALSE, sep="\t")

#moduleHub
muscle_moduleBlue_BP <- enrichGO(gene          = muscle_moduleBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = muscle_universe_egid,
                readable = TRUE)
write.table(muscle_moduleBlue_BP@result[muscle_moduleBlue_BP@result$p.adjust <0.05,], "muscle_moduleBlue_BP.txt", quote=FALSE, sep="\t")

muscle_moduleBlue_BP_simp <- clusterProfiler::simplify(muscle_moduleBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(muscle_moduleBlue_BP_simp@result[muscle_moduleBlue_BP_simp@result$p.adjust <0.05,], "muscle_moduleBlue_BP_simp.txt", quote=FALSE, sep="\t")

muscle_moduleBlue_kegg <- enrichKEGG(gene          = muscle_moduleBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = muscle_universe_egid,
              pvalueCutoff  = 0.05)
write.table(muscle_moduleBlue_kegg@result[muscle_moduleBlue_kegg@result$p.adjust <0.05,], file="muscle_moduleBlue_kegg.txt", quote=FALSE, sep="\t")

muscle_moduleBlue_react <- enrichPathway(gene          = muscle_moduleBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = muscle_universe_egid,
              readable = TRUE)
write.table(muscle_moduleBlue_react@result[muscle_moduleBlue_react@result$p.adjust <0.05,], file="muscle_moduleBlue_react.txt", quote=FALSE, sep="\t")
```

#### Liver

``` r
### TAN module
#hubTan
liver_hubTan_BP <- enrichGO(gene          = liver_hubTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubTan_BP@result[liver_hubTan_BP@result$p.adjust <0.05,], "liver_hubTan_BP.txt", quote=FALSE, sep="\t")

liver_hubTan_BP_simp <- clusterProfiler::simplify(liver_hubTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubTan_BP_simp@result[liver_hubTan_BP_simp@result$p.adjust <0.05,], "liver_hubTan_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubTan_kegg <- enrichKEGG(gene          = liver_hubTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubTan_kegg@result[liver_hubTan_kegg@result$p.adjust <0.05,], file="liver_hubTan_kegg.txt", quote=FALSE, sep="\t")

liver_hubTan_react <- enrichPathway(gene          = liver_hubTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubTan_react@result[liver_hubTan_react@result$p.adjust <0.05,], file="liver_hubTan_react.txt", quote=FALSE, sep="\t")

#moduleTan
liver_moduleTan_BP <- enrichGO(gene          = liver_moduleTan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleTan_BP@result[liver_moduleTan_BP@result$p.adjust <0.05,], "liver_moduleTan_BP.txt", quote=FALSE, sep="\t")

liver_moduleTan_BP_simp <- clusterProfiler::simplify(liver_moduleTan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleTan_BP_simp@result[liver_moduleTan_BP_simp@result$p.adjust <0.05,], "liver_moduleTan_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleTan_kegg <- enrichKEGG(gene          = liver_moduleTan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleTan_kegg@result[liver_moduleTan_kegg@result$p.adjust <0.05,], file="liver_moduleTan_kegg.txt", quote=FALSE, sep="\t")

liver_moduleTan_react <- enrichPathway(gene          = liver_moduleTan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleTan_react@result[liver_moduleTan_react@result$p.adjust <0.05,], file="liver_moduleTan_react.txt", quote=FALSE, sep="\t")

### RED module
#hubRed
liver_hubRed_BP <- enrichGO(gene          = liver_hubRed_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubRed_BP@result[liver_hubRed_BP@result$p.adjust <0.05,], "liver_hubRed_BP.txt", quote=FALSE, sep="\t")

liver_hubRed_BP_simp <- clusterProfiler::simplify(liver_hubRed_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubRed_BP_simp@result[liver_hubRed_BP_simp@result$p.adjust <0.05,], "liver_hubRed_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubRed_kegg <- enrichKEGG(gene          = liver_hubRed_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubRed_kegg@result[liver_hubRed_kegg@result$p.adjust <0.05,], file="liver_hubRed_kegg.txt", quote=FALSE, sep="\t")

liver_hubRed_react <- enrichPathway(gene          = liver_hubRed_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubRed_react@result[liver_hubRed_react@result$p.adjust <0.05,], file="liver_hubRed_react.txt", quote=FALSE, sep="\t")

#moduleRed
liver_moduleRed_BP <- enrichGO(gene          = liver_moduleRed_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleRed_BP@result[liver_moduleRed_BP@result$p.adjust <0.05,], "liver_moduleRed_BP.txt", quote=FALSE, sep="\t")

liver_moduleRed_BP_simp <- clusterProfiler::simplify(liver_moduleRed_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleRed_BP_simp@result[liver_moduleRed_BP_simp@result$p.adjust <0.05,], "liver_moduleRed_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleRed_kegg <- enrichKEGG(gene          = liver_moduleRed_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleRed_kegg@result[liver_moduleRed_kegg@result$p.adjust <0.05,], file="liver_moduleRed_kegg.txt", quote=FALSE, sep="\t")

liver_moduleRed_react <- enrichPathway(gene          = liver_moduleRed_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleRed_react@result[liver_moduleRed_react@result$p.adjust <0.05,], file="liver_moduleRed_react.txt", quote=FALSE, sep="\t")

### Darkolivegreen module
#hubDarkolivegreen
liver_hubDarkolivegreen_BP <- enrichGO(gene          = liver_hubDarkolivegreen_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubDarkolivegreen_BP@result[liver_hubDarkolivegreen_BP@result$p.adjust <0.05,], "liver_hubDarkolivegreen_BP.txt", quote=FALSE, sep="\t")

liver_hubDarkolivegreen_BP_simp <- clusterProfiler::simplify(liver_hubDarkolivegreen_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubDarkolivegreen_BP_simp@result[liver_hubDarkolivegreen_BP_simp@result$p.adjust <0.05,], "liver_hubDarkolivegreen_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubDarkolivegreen_kegg <- enrichKEGG(gene          = liver_hubDarkolivegreen_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubDarkolivegreen_kegg@result[liver_hubDarkolivegreen_kegg@result$p.adjust <0.05,], file="liver_hubDarkolivegreen_kegg.txt", quote=FALSE, sep="\t")

liver_hubDarkolivegreen_react <- enrichPathway(gene          = liver_hubDarkolivegreen_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubDarkolivegreen_react@result[liver_hubDarkolivegreen_react@result$p.adjust <0.05,], file="liver_hubDarkolivegreen_react.txt", quote=FALSE, sep="\t")

#moduleDarkolivegreen
liver_moduleDarkolivegreen_BP <- enrichGO(gene          = liver_moduleDarkolivegreen_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleDarkolivegreen_BP@result[liver_moduleDarkolivegreen_BP@result$p.adjust <0.05,], "liver_moduleDarkolivegreen_BP.txt", quote=FALSE, sep="\t")

liver_moduleDarkolivegreen_BP_simp <- clusterProfiler::simplify(liver_moduleDarkolivegreen_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleDarkolivegreen_BP_simp@result[liver_moduleDarkolivegreen_BP_simp@result$p.adjust <0.05,], "liver_moduleDarkolivegreen_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleDarkolivegreen_kegg <- enrichKEGG(gene          = liver_moduleDarkolivegreen_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleDarkolivegreen_kegg@result[liver_moduleDarkolivegreen_kegg@result$p.adjust <0.05,], file="liver_moduleDarkolivegreen_kegg.txt", quote=FALSE, sep="\t")

liver_moduleDarkolivegreen_react <- enrichPathway(gene          = liver_moduleDarkolivegreen_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleDarkolivegreen_react@result[liver_moduleDarkolivegreen_react@result$p.adjust <0.05,], file="liver_moduleDarkolivegreen_react.txt", quote=FALSE, sep="\t")

### DARKGREY module
#hubDarkgrey
liver_hubDarkgrey_BP <- enrichGO(gene          = liver_hubDarkgrey_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubDarkgrey_BP@result[liver_hubDarkgrey_BP@result$p.adjust <0.05,], "liver_hubDarkgrey_BP.txt", quote=FALSE, sep="\t")

liver_hubDarkgrey_BP_simp <- clusterProfiler::simplify(liver_hubDarkgrey_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubDarkgrey_BP_simp@result[liver_hubDarkgrey_BP_simp@result$p.adjust <0.05,], "liver_hubDarkgrey_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubDarkgrey_kegg <- enrichKEGG(gene          = liver_hubDarkgrey_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubDarkgrey_kegg@result[liver_hubDarkgrey_kegg@result$p.adjust <0.05,], file="liver_hubDarkgrey_kegg.txt", quote=FALSE, sep="\t")

liver_hubDarkgrey_react <- enrichPathway(gene          = liver_hubDarkgrey_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubDarkgrey_react@result[liver_hubDarkgrey_react@result$p.adjust <0.05,], file="liver_hubDarkgrey_react.txt", quote=FALSE, sep="\t")

#moduleDarkgrey
liver_moduleDarkgrey_BP <- enrichGO(gene          = liver_moduleDarkgrey_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleDarkgrey_BP@result[liver_moduleDarkgrey_BP@result$p.adjust <0.05,], "liver_moduleDarkgrey_BP.txt", quote=FALSE, sep="\t")

liver_moduleDarkgrey_BP_simp <- clusterProfiler::simplify(liver_moduleDarkgrey_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleDarkgrey_BP_simp@result[liver_moduleDarkgrey_BP_simp@result$p.adjust <0.05,], "liver_moduleDarkgrey_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleDarkgrey_kegg <- enrichKEGG(gene          = liver_moduleDarkgrey_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleDarkgrey_kegg@result[liver_moduleDarkgrey_kegg@result$p.adjust <0.05,], file="liver_moduleDarkgrey_kegg.txt", quote=FALSE, sep="\t")

liver_moduleDarkgrey_react <- enrichPathway(gene          = liver_moduleDarkgrey_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleDarkgrey_react@result[liver_moduleDarkgrey_react@result$p.adjust <0.05,], file="liver_moduleDarkgrey_react.txt", quote=FALSE, sep="\t")

### CYAN module
#hubCyan
liver_hubCyan_BP <- enrichGO(gene          = liver_hubCyan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubCyan_BP@result[liver_hubCyan_BP@result$p.adjust <0.05,], "liver_hubCyan_BP.txt", quote=FALSE, sep="\t")

liver_hubCyan_BP_simp <- clusterProfiler::simplify(liver_hubCyan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubCyan_BP_simp@result[liver_hubCyan_BP_simp@result$p.adjust <0.05,], "liver_hubCyan_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubCyan_kegg <- enrichKEGG(gene          = liver_hubCyan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubCyan_kegg@result[liver_hubCyan_kegg@result$p.adjust <0.05,], file="liver_hubCyan_kegg.txt", quote=FALSE, sep="\t")

liver_hubCyan_react <- enrichPathway(gene          = liver_hubCyan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubCyan_react@result[liver_hubCyan_react@result$p.adjust <0.05,], file="liver_hubCyan_react.txt", quote=FALSE, sep="\t")

#moduleCyan
liver_moduleCyan_BP <- enrichGO(gene          = liver_moduleCyan_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleCyan_BP@result[liver_moduleCyan_BP@result$p.adjust <0.05,], "liver_moduleCyan_BP.txt", quote=FALSE, sep="\t")

liver_moduleCyan_BP_simp <- clusterProfiler::simplify(liver_moduleCyan_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleCyan_BP_simp@result[liver_moduleCyan_BP_simp@result$p.adjust <0.05,], "liver_moduleCyan_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleCyan_kegg <- enrichKEGG(gene          = liver_moduleCyan_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleCyan_kegg@result[liver_moduleCyan_kegg@result$p.adjust <0.05,], file="liver_moduleCyan_kegg.txt", quote=FALSE, sep="\t")

liver_moduleCyan_react <- enrichPathway(gene          = liver_moduleCyan_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleCyan_react@result[liver_moduleCyan_react@result$p.adjust <0.05,], file="liver_moduleCyan_react.txt", quote=FALSE, sep="\t")

### BROWN module
#hubBrown
liver_hubBrown_BP <- enrichGO(gene          = liver_hubBrown_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubBrown_BP@result[liver_hubBrown_BP@result$p.adjust <0.05,], "liver_hubBrown_BP.txt", quote=FALSE, sep="\t")

liver_hubBrown_BP_simp <- clusterProfiler::simplify(liver_hubBrown_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubBrown_BP_simp@result[liver_hubBrown_BP_simp@result$p.adjust <0.05,], "liver_hubBrown_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubBrown_kegg <- enrichKEGG(gene          = liver_hubBrown_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubBrown_kegg@result[liver_hubBrown_kegg@result$p.adjust <0.05,], file="liver_hubBrown_kegg.txt", quote=FALSE, sep="\t")

liver_hubBrown_react <- enrichPathway(gene          = liver_hubBrown_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubBrown_react@result[liver_hubBrown_react@result$p.adjust <0.05,], file="liver_hubBrown_react.txt", quote=FALSE, sep="\t")

#moduleBrown
liver_moduleBrown_BP <- enrichGO(gene          = liver_moduleBrown_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleBrown_BP@result[liver_moduleBrown_BP@result$p.adjust <0.05,], "liver_moduleBrown_BP.txt", quote=FALSE, sep="\t")

liver_moduleBrown_BP_simp <- clusterProfiler::simplify(liver_moduleBrown_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleBrown_BP_simp@result[liver_moduleBrown_BP_simp@result$p.adjust <0.05,], "liver_moduleBrown_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleBrown_kegg <- enrichKEGG(gene          = liver_moduleBrown_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleBrown_kegg@result[liver_moduleBrown_kegg@result$p.adjust <0.05,], file="liver_moduleBrown_kegg.txt", quote=FALSE, sep="\t")

liver_moduleBrown_react <- enrichPathway(gene          = liver_moduleBrown_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleBrown_react@result[liver_moduleBrown_react@result$p.adjust <0.05,], file="liver_moduleBrown_react.txt", quote=FALSE, sep="\t")

### BLUE module
#hubBlue
liver_hubBlue_BP <- enrichGO(gene          = liver_hubBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_hubBlue_BP@result[liver_hubBlue_BP@result$p.adjust <0.05,], "liver_hubBlue_BP.txt", quote=FALSE, sep="\t")

liver_hubBlue_BP_simp <- clusterProfiler::simplify(liver_hubBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_hubBlue_BP_simp@result[liver_hubBlue_BP_simp@result$p.adjust <0.05,], "liver_hubBlue_BP_simp.txt", quote=FALSE, sep="\t")

liver_hubBlue_kegg <- enrichKEGG(gene          = liver_hubBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_hubBlue_kegg@result[liver_hubBlue_kegg@result$p.adjust <0.05,], file="liver_hubBlue_kegg.txt", quote=FALSE, sep="\t")

liver_hubBlue_react <- enrichPathway(gene          = liver_hubBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_hubBlue_react@result[liver_hubBlue_react@result$p.adjust <0.05,], file="liver_hubBlue_react.txt", quote=FALSE, sep="\t")

#moduleBlue
liver_moduleBlue_BP <- enrichGO(gene          = liver_moduleBlue_egid,
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                universe = liver_universe_egid,
                readable = TRUE)
write.table(liver_moduleBlue_BP@result[liver_moduleBlue_BP@result$p.adjust <0.05,], "liver_moduleBlue_BP.txt", quote=FALSE, sep="\t")

liver_moduleBlue_BP_simp <- clusterProfiler::simplify(liver_moduleBlue_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
write.table(liver_moduleBlue_BP_simp@result[liver_moduleBlue_BP_simp@result$p.adjust <0.05,], "liver_moduleBlue_BP_simp.txt", quote=FALSE, sep="\t")

liver_moduleBlue_kegg <- enrichKEGG(gene          = liver_moduleBlue_egid,
              keyType       = "ncbi-geneid",
              organism         = "mmu",
              pAdjustMethod = "BH",
              universe = liver_universe_egid,
              pvalueCutoff  = 0.05)
write.table(liver_moduleBlue_kegg@result[liver_moduleBlue_kegg@result$p.adjust <0.05,], file="liver_moduleBlue_kegg.txt", quote=FALSE, sep="\t")

liver_moduleBlue_react <- enrichPathway(gene          = liver_moduleBlue_egid,
              organism         = "mouse",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              universe = liver_universe_egid,
              readable = TRUE)
write.table(liver_moduleBlue_react@result[liver_moduleBlue_react@result$p.adjust <0.05,], file="liver_moduleBlue_react.txt", quote=FALSE, sep="\t")
```

# 4. Upset plots

## Fig 5 - right

``` r
# create list comprising each tissue's DEG-HUB overlap; in tissues with more than one module, we merged the vectors of each module and extracted the unique elements 

#TRENDY-HUB  INTERSECTION
allTissues_hub <- list(brain_hubTan_BP@result$ID[brain_hubTan_BP@result$p.adjust<0.05],
                       
                      c(heart_hubTan_BP@result$ID[heart_hubTan_BP@result$p.adjust<0.05],
                        heart_hubBlue_BP@result$ID[heart_hubBlue_BP@result$p.adjust<0.05]),
                      
                      c(liver_hubSalmon_BP@result$ID[liver_hubSalmon_BP@result$p.adjust<0.05],
                        liver_hubDarkturquoise_BP@result$ID[liver_hubDarkturquoise_BP@result$p.adjust<0.05],
                        liver_hubBlue_BP@result$ID[liver_hubBlue_BP@result$p.adjust<0.05],
                        liver_hubCyan_BP@result$ID[liver_hubCyan_BP@result$p.adjust<0.05],
                        liver_hubDarkgrey_BP@result$ID[liver_hubDarkgrey_BP@result$p.adjust<0.05],
                        liver_hubTan_BP@result$ID[liver_hubTan_BP@result$p.adjust<0.05]),
                      
                      c(muscle_hubMagenta_BP@result$ID[muscle_hubMagenta_BP@result$p.adjust<0.05],
                        muscle_hubBrown_BP@result$ID[muscle_hubBrown_BP@result$p.adjust<0.05],
                        muscle_hubBlue_BP@result$ID[muscle_hubBlue_BP@result$p.adjust<0.05]))

names(allTissues_hub) <- c('Brain','Heart','Liver','Muscle')
```

The code below was adapted from [here](http://research.libd.org/rstatsclub/post/hacking-our-way-through-upsetr/#.YG2YGuhKiUk).

``` r
nsets = 5; nintersects = 40; sets = NULL;
  keep.order = F; set.metadata = NULL; intersections = NULL;
  matrix.color = "gray23"; main.bar.color = "gray23";
  mainbar.y.label = "Intersection Size"; mainbar.y.max = NULL;
  sets.bar.color = "gray23"; sets.x.label = "Set Size";
  point.size = 2.2; line.size = 0.7; mb.ratio = c(0.7, 0.3);
  expression = NULL; att.pos = NULL; att.color = main.bar.color;
  order.by = c("freq", "degree"); decreasing = c(T, F);
  show.numbers = "yes"; number.angles = 0; group.by = "degree";
  cutoff = NULL; queries = NULL; query.legend = "none";
  shade.color = "gray88"; shade.alpha = 0.25; matrix.dot.alpha = 0.5;
  empty.intersections = NULL; color.pal = 1; boxplot.summary = NULL;
  attribute.plots = NULL; scale.intersections = "identity";
  scale.sets = "identity"; text.scale = 1; set_size.angles = 0;
  set_size.show = FALSE; set_size.numbers_size = NULL;
  set_size.scale_max = NULL

data=fromList(allTissues_hub);
sets = c("Brain","Heart","Liver","Muscle");
order.by = "freq";
text.scale = c(2.5,2,2,2,2.5,3);
point.size = 8; 
line.size = 1; 
sets.bar.color= c("#EEDD88","#77AADD","#FFAABB","#EE8866")#Liver-Brain-Muscle-Heart


startend <- UpSetR:::FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]

  if(color.pal == 1){
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF")
  } else{
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                 "#CC79A7")
  }

  if(is.null(intersections) == F){
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::specific_intersections(data, first.col, last.col, intersections, order.by, group.by, decreasing,
                                        cutoff, main.bar.color, Set_names)
  } else if(is.null(intersections) == T){
    Set_names <- sets
    if(is.null(Set_names) == T || length(Set_names) == 0 ){
      Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, nsets)
    }
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
    Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, Set_names, nintersects, main.bar.color,
                         order.by, group.by, cutoff, empty.intersections, decreasing)
  }
  Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
  labels <- UpSetR:::Make_labels(Matrix_setup)
  
  att.x <- c(); att.y <- c();
  if(is.null(attribute.plots) == F){
    for(i in seq_along(attribute.plots$plots)){
      if(length(attribute.plots$plots[[i]]$x) != 0){
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if(length(attribute.plots$plots[[i]]$x) == 0){
        att.x[i] <- NA
      }
      if(length(attribute.plots$plots[[i]]$y) != 0){
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if(length(attribute.plots$plots[[i]]$y) == 0){
        att.y[i] <- NA
      }
    }
  }

  BoxPlots <- NULL
  if(is.null(boxplot.summary) == F){
    BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, Set_names)
    BoxPlots <- list()
    for(i in seq_along(boxplot.summary)){
      BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], att.color)
    }
  }

  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if(is.null(queries) == F){
    custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
    customDat <- UpSetR:::customQueries(New_data, custom.queries, Set_names)
    legend <- UpSetR:::GuideGenerator(queries, palette)
    legend <- UpSetR:::Make_legend(legend)
    if(is.null(att.x) == F && is.null(customDat) == F){
      customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
    }
    customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, All_Freqs, custom.queries)
  }
  if(is.null(queries) == F){
    Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
    Matrix_col <- UpSetR:::intersects(QuerieInterData, Intersection, New_data, first.col, Num_of_set,
                             All_Freqs, expression, Set_names, palette)
    Element <- UpSetR:::SeperateQueries(queries, 1, palette)
    EBar_data <-UpSetR:::ElemBarDat(Element, New_data, first.col, expression, Set_names,palette, All_Freqs)
  } else{
    Matrix_col <- NULL
  }
  
  Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, Matrix_col, matrix.dot.alpha)
  
  # #Liver-Brain-Heart-Muscle
   for(i in 1:4) {
       j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
       if(length(j) > 0) Matrix_layout$color[j] <- c("#EEDD88","#77AADD","#FFAABB","#EE8866")[i]
   }

  
  
Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #EEDD88   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  4 1     0  gray83   0.5          4No
    ## 5  1 2     0  gray83   0.5          5No
    ## 6  2 2     1 #77AADD   1.0         2yes
    ## 7  3 2     0  gray83   0.5          7No
    ## 8  4 2     0  gray83   0.5          8No
    ## 9  1 3     0  gray83   0.5          9No
    ## 10 2 3     0  gray83   0.5         10No
    ## 11 3 3     1 #FFAABB   1.0         3yes
    ## 12 4 3     0  gray83   0.5         12No
    ## 13 1 4     0  gray83   0.5         13No
    ## 14 2 4     0  gray83   0.5         14No
    ## 15 3 4     0  gray83   0.5         15No
    ## 16 4 4     1 #EE8866   1.0         4yes
    ## 17 1 5     1 #EEDD88   1.0         5yes
    ## 18 2 5     1 #77AADD   1.0         5yes
    ## 19 3 5     0  gray83   0.5         19No
    ## 20 4 5     0  gray83   0.5         20No
    ## 21 1 6     0  gray83   0.5         21No
    ## 22 2 6     0  gray83   0.5         22No
    ## 23 3 6     1 #FFAABB   1.0         6yes
    ## 24 4 6     1 #EE8866   1.0         6yes
    ## 25 1 7     0  gray83   0.5         25No
    ## 26 2 7     1 #77AADD   1.0         7yes
    ## 27 3 7     1 #FFAABB   1.0         7yes
    ## 28 4 7     0  gray83   0.5         28No
    ## 29 1 8     1 #EEDD88   1.0         8yes
    ## 30 2 8     0  gray83   0.5         30No
    ## 31 3 8     0  gray83   0.5         31No
    ## 32 4 8     1 #EE8866   1.0         8yes
    ## 33 1 9     1 #EEDD88   1.0         9yes
    ## 34 2 9     0  gray83   0.5         34No
    ## 35 3 9     1 #FFAABB   1.0         9yes
    ## 36 4 9     0  gray83   0.5         36No

``` r
 Set_sizes <- UpSetR:::FindSetFreqs(New_data, first.col, Num_of_set, Set_names, keep.order)
  Bar_Q <- NULL
  if(is.null(queries) == F){
    Bar_Q <- UpSetR:::intersects(QuerieInterBar, Intersection, New_data, first.col, Num_of_set, All_Freqs, expression, Set_names, palette)
  }
  QInter_att_data <- NULL
  QElem_att_data <- NULL
  if((is.null(queries) == F) & (is.null(att.x) == F)){
    QInter_att_data <- UpSetR:::intersects(QuerieInterAtt, Intersection, New_data, first.col, Num_of_set, att.x, att.y,
                                  expression, Set_names, palette)
    QElem_att_data <- UpSetR:::elements(QuerieElemAtt, Element, New_data, first.col, expression, Set_names, att.x, att.y,
                               palette)
  }
  AllQueryData <- UpSetR:::combineQueriesData(QInter_att_data, QElem_att_data, customAttDat, att.x, att.y)

  ShadingData <- NULL

  if(is.null(set.metadata) == F){
    ShadingData <- get_shade_groups(set.metadata, Set_names, Matrix_layout, shade.alpha)
    output <- Make_set_metadata_plot(set.metadata, Set_names)
    set.metadata.plots <- output[[1]]
    set.metadata <- output[[2]]

    if(is.null(ShadingData) == FALSE){
    shade.alpha <- unique(ShadingData$alpha)
    }
  } else {
    set.metadata.plots <- NULL
  }
  if(is.null(ShadingData) == TRUE){
  ShadingData <- UpSetR:::MakeShading(Matrix_layout, shade.color)
  }
  Main_bar <- suppressMessages(UpSetR:::Make_main_bar(All_Freqs, Bar_Q, show.numbers, mb.ratio, customQBar, number.angles, EBar_data, mainbar.y.label,
                            mainbar.y.max, scale.intersections, text.scale, attribute.plots))
  Matrix <- UpSetR:::Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, point.size, line.size,
                             text.scale, labels, ShadingData, shade.alpha)
  Sizes <- UpSetR:::Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, sets.x.label, scale.sets, text.scale, set_size.angles,set_size.show,set_size.numbers_size, set_size.scale_max)

  structure(class = "upset",
    .Data=list(
      Main_bar = Main_bar,
      Matrix = Matrix,
      Sizes = Sizes,
      labels = labels,
      mb.ratio = mb.ratio,
      att.x = att.x,
      att.y = att.y,
      New_data = New_data,
      expression = expression,
      att.pos = att.pos,
      first.col = first.col,
      att.color = att.color,
      AllQueryData = AllQueryData,
      attribute.plots = attribute.plots,
      legend = legend,
      query.legend = query.legend,
      BoxPlots = BoxPlots,
      Set_names = Set_names,
      set.metadata = set.metadata,
      set.metadata.plots = set.metadata.plots)
  )
```

![](functional_analysis_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
  Make_matrix_plot <- function(Mat_data,Set_size_data, Main_bar_data, point_size, line_size, text_scale, labels,
                             shading_data, shade_alpha){

  if(length(text_scale) == 1){
    name_size_scale <- text_scale
  }
  if(length(text_scale) > 1 && length(text_scale) <= 6){
    name_size_scale <- text_scale[5]
  }
  
  Mat_data$line_col <- 'black'

  Matrix_plot <- (ggplot()
                  + theme(panel.background = element_rect(fill = "white"),
                          plot.margin=unit(c(-0.2,0.5,0.5,0.5), "lines"),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_text(colour = "gray0",
                                                     size = 7*name_size_scale, hjust = 0.4),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                  + xlab(NULL) + ylab("   ")
                  + scale_y_continuous(breaks = c(1:nrow(Set_size_data)),
                                       limits = c(0.5,(nrow(Set_size_data) +0.5)),
                                       labels = labels, expand = c(0,0))
                  + scale_x_continuous(limits = c(0,(nrow(Main_bar_data)+1 )), expand = c(0,0))
                  + geom_rect(data = shading_data, aes_string(xmin = "min", xmax = "max",
                                                              ymin = "y_min", ymax = "y_max"),
                              fill = shading_data$shade_color, alpha = shade_alpha)
                  + geom_line(data= Mat_data, aes_string(group = "Intersection", x="x", y="y",
                                                         colour = "line_col"), size = line_size)
                 + geom_point(data= Mat_data, aes_string(x= "x", y= "y"), colour = Mat_data$color,
                     size= point_size, alpha = Mat_data$alpha, shape=16)
                  + scale_color_identity())
  Matrix_plot <- ggplot_gtable(ggplot_build(Matrix_plot))
  return(Matrix_plot)
  }
  
  Matrix <- Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, point.size, line.size,
                             text.scale, labels, ShadingData, shade.alpha)
  Sizes <- UpSetR:::Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, sets.x.label, scale.sets, text.scale, set_size.angles,set_size.show,set_size.numbers_size, set_size.scale_max)

  structure(class = "upset",
    .Data=list(
      Main_bar = Main_bar,
      Matrix = Matrix,
      Sizes = Sizes,
      labels = labels,
      mb.ratio = mb.ratio,
      att.x = att.x,
      att.y = att.y,
      New_data = New_data,
      expression = expression,
      att.pos = att.pos,
      first.col = first.col,
      att.color = att.color,
      AllQueryData = AllQueryData,
      attribute.plots = attribute.plots,
      legend = legend,
      query.legend = query.legend,
      BoxPlots = BoxPlots,
      Set_names = Set_names,
      set.metadata = set.metadata,
      set.metadata.plots = set.metadata.plots)
  )
```

![](functional_analysis_files/figure-markdown_github/unnamed-chunk-20-2.png)

``` r
tiff('Fig-5-right-hub.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
structure(class = "upset",
    .Data=list(
      Main_bar = Main_bar,
      Matrix = Matrix,
      Sizes = Sizes,
      labels = labels,
      mb.ratio = mb.ratio,
      att.x = att.x,
      att.y = att.y,
      New_data = New_data,
      expression = expression,
      att.pos = att.pos,
      first.col = first.col,
      att.color = att.color,
      AllQueryData = AllQueryData,
      attribute.plots = attribute.plots,
      legend = legend,
      query.legend = query.legend,
      BoxPlots = BoxPlots,
      Set_names = Set_names,
      set.metadata = set.metadata,
      set.metadata.plots = set.metadata.plots))
grid.text("GO term overlap", gp=gpar(fontsize=25), x = 0.65, y=0.93)
grid.text("*", gp=gpar(fontsize=30), x = 0.648, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.715, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.782, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.855, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.92, y=0.05)
dev.off()
```

    ## png 
    ##   2

The code below allows to fetch the genes present in each overlap and was retrieved from [here](https://github.com/hms-dbmi/UpSetR/issues/85).

``` r
###function
overlapGroups <- function (listInput, sort = TRUE) {
  listInputmat    <- fromList(listInput) == 1
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}
```

``` r
#TRENDY-HUB
list_allTissues_hub <- overlapGroups(allTissues_hub)
names(list_allTissues_hub)
```

    ## [1] "Liver"        "Brain"        "Muscle"       "Heart"        "Brain:Liver" 
    ## [6] "Brain:Muscle" "Heart:Muscle" "Heart:Liver"  "Liver:Muscle"

``` r
list1_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[1]]]
list2_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[2]]]
list3_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[3]]]
list4_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[4]]]
list5_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[5]]]
list6_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[6]]]
list7_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[7]]]
list8_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[8]]]
list9_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[9]]]

upsetList_allTissues_hub <- list(list1_allTissues_hub,list2_allTissues_hub,list3_allTissues_hub,list4_allTissues_hub,list5_allTissues_hub,list6_allTissues_hub,list7_allTissues_hub,list8_allTissues_hub,list9_allTissues_hub)

upsetList_allTissues_hub <- as.data.frame(sapply(upsetList_allTissues_hub, '[', seq(max(sapply(upsetList_allTissues_hub, length)))))
colnames(upsetList_allTissues_hub) <- c("Liver","Brain","Muscle","Heart","Brain:Liver","Brain:Muscle","Heart:Muscle","Heart:Liver","Liver:Muscle")

write.table(upsetList_allTissues_hub, file="upsetList_BP_allTissues_hub.txt", row.names = F, quote=FALSE, sep="\t")
```

Yu, Guangchuang, Li Gen Wang, Yanyan Han, and Qing Yu He. 2012. “ClusterProfiler: An R package for comparing biological themes among gene clusters.” *OMICS A Journal of Integrative Biology* 16 (5): 284–87. doi:[10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118).
