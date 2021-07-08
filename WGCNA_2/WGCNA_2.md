WGCNA - Module-trait associations
================
Margarida Ferreira
Last updated July 2021

-   [1. Initial setup](#initial-setup)
-   [2. Expression and network data input](#expression-and-network-data-input)
    -   [Brain](#brain)
    -   [Heart](#heart)
    -   [Muscle](#muscle)
    -   [Liver](#liver)
    -   [Pancreas](#pancreas)
-   [3. Relating modules to external information and identifying important genes](#relating-modules-to-external-information-and-identifying-important-genes)
    -   [3.1 Calculation of module-trait associations](#calculation-of-module-trait-associations)
        -   [Brain](#brain-1)
        -   [Heart](#heart-1)
        -   [Muscle](#muscle-1)
        -   [Liver](#liver-1)
        -   [Pancreas](#pancreas-1)
        -   [Fig 2A](#fig-2a)
    -   [3.2 Significant module-age (and module-sex) associations](#significant-module-age-and-module-sex-associations)
        -   [Fig 3A/Fig S2](#fig-3afig-s2)
        -   [Brain](#brain-2)
        -   [Heart](#heart-2)
        -   [Muscle](#muscle-2)
        -   [Liver](#liver-2)
        -   [Pancreas](#pancreas-2)
    -   [3.3 Gene relationship to trait and important modules: Gene Significance and Module Membership](#gene-relationship-to-trait-and-important-modules-gene-significance-and-module-membership)
        -   [Brain](#brain-3)
        -   [Heart](#heart-3)
        -   [Muscle](#muscle-3)
        -   [Liver](#liver-3)
        -   [Pancreas](#pancreas-3)
        -   [Fig 2B](#fig-2b)
        -   [Fig 3B /Fig S3](#fig-3b-fig-s3)
    -   [3.4 Intramodular analysis: identifying genes with high GS and MM](#intramodular-analysis-identifying-genes-with-high-gs-and-mm)
        -   [Brain](#brain-4)
        -   [Heart](#heart-4)
        -   [Muscle](#muscle-4)
        -   [Liver](#liver-4)
        -   [Pancreas](#pancreas-4)
-   [4. Upset plots](#upset-plots)
    -   [Fig 4/Fig S4](#fig-4fig-s4)
        -   [Brain](#brain-5)
        -   [Heart](#heart-5)
        -   [Muscle](#muscle-5)
        -   [Liver](#liver-5)
    -   [Fig 5 - left](#fig-5---left)
-   [Session Info](#session-info)

# 1. Initial setup

Loading the required packages:

``` r
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
library(ggplot2)
library(grid)
library(UpSetR)
```

# 2. Expression and network data input

## Brain

``` r
# Expression and trait data 
load(file = "Brain-dataInput.RData");

# Load network data 
load(file = "Brain-networkConstruction-auto.RData"); # default parameters: deepsplit=2; mergeCutHeight=0.15

# Check the sample order
all(rownames(brain_traitData) %in% rownames(brainDataExpr))
```

    ## [1] TRUE

``` r
all(rownames(brain_traitData) == rownames(brainDataExpr))
```

    ## [1] TRUE

``` r
# Convert Age and Sex variables into factors
brain_traitData$Age <- factor(brain_traitData$Age) # 3,6,9,12,15,18,21,24,27
brain_traitData$Sex <- factor(brain_traitData$Sex) # Female == 0; Male == 1

# Order the expression data frame by increasing time point
brain_traitData <- brain_traitData[order(brain_traitData$Age, decreasing=FALSE),]
brainDataExpr <- brainDataExpr[rownames(brain_traitData),]

# Check the number and size of the identified network modules
length(unique(brain_moduleColors))
```

    ## [1] 25

``` r
table(brain_moduleColors)
```

    ## brain_moduleColors
    ##         black          blue         brown          cyan     darkgreen 
    ##           411          1818          1006           120            69 
    ##      darkgrey       darkred darkturquoise         green   greenyellow 
    ##            59            73            64           778           283 
    ##          grey        grey60     lightcyan    lightgreen   lightyellow 
    ##         17870            97           108            97            85 
    ##       magenta  midnightblue          pink        purple           red 
    ##           369           108           372           368           613 
    ##     royalblue        salmon           tan     turquoise        yellow 
    ##            82           129           220          8126           839

``` r
utils::capture.output(length(unique(brain_moduleColors)),table(brain_moduleColors), file = "brain_modules.txt")
```

## Heart

``` r
# Expression and trait data 
load(file = "Heart-dataInput.RData");

# Load network data 
load(file = "Heart-networkConstruction-auto.RData"); # default parameters: deepsplit=2; mergeCutHeight=0.15

# Check the sample order
all(rownames(heart_traitData) %in% rownames(heartDataExpr))
```

    ## [1] TRUE

``` r
all(rownames(heart_traitData) == rownames(heartDataExpr))
```

    ## [1] TRUE

``` r
# Convert Age and Sex variables into factors
heart_traitData$Age <- factor(heart_traitData$Age) # 3,6,9,12,15,18,21,24,27
heart_traitData$Sex <- factor(heart_traitData$Sex) # Female == 0; Male == 1

# Order the expression data frame by increasing time point
heart_traitData <- heart_traitData[order(heart_traitData$Age, decreasing=FALSE),]
heartDataExpr <- heartDataExpr[rownames(heart_traitData),]

# Check the number and size of the identified network modules
length(unique(heart_moduleColors))
```

    ## [1] 20

``` r
table(heart_moduleColors)
```

    ## heart_moduleColors
    ##        black         blue        brown         cyan        green  greenyellow 
    ##          245         1209          904           99          415          135 
    ##         grey       grey60    lightcyan   lightgreen  lightyellow      magenta 
    ##        14055           64           71           63           61          154 
    ## midnightblue         pink       purple          red       salmon          tan 
    ##           90          181          143          257          107          125 
    ##    turquoise       yellow 
    ##         9039          656

``` r
utils::capture.output(length(unique(heart_moduleColors)),table(heart_moduleColors), file = "heart_modules.txt")
```

## Muscle

``` r
# Expression and trait data 
load(file = "Muscle-dataInput.RData");

# Load network data 
load(file = "Muscle-networkConstruction-auto.RData"); # default parameters: deepsplit=2; mergeCutHeight=0.15

# Check the sample order
all(rownames(muscle_traitData) %in% rownames(muscleDataExpr))
```

    ## [1] TRUE

``` r
all(rownames(muscle_traitData) == rownames(muscleDataExpr))
```

    ## [1] TRUE

``` r
# Convert Age and Sex variables into factors
muscle_traitData$Age <- factor(muscle_traitData$Age) # 3,6,9,12,15,18,21,24,27
muscle_traitData$Sex <- factor(muscle_traitData$Sex) # Female == 0; Male == 1

# Order the expression data frame by increasing time point
muscle_traitData <- muscle_traitData[order(muscle_traitData$Age, decreasing=FALSE),]
muscleDataExpr <- muscleDataExpr[rownames(muscle_traitData),]

# Check the number and size of the identified network modules
length(unique(muscle_moduleColors))
```

    ## [1] 22

``` r
table(muscle_moduleColors)
```

    ## muscle_moduleColors
    ##        black         blue        brown         cyan      darkred        green 
    ##          241          949          806           88           59          518 
    ##  greenyellow         grey       grey60    lightcyan   lightgreen  lightyellow 
    ##          140         9336           81           84           73           64 
    ##      magenta midnightblue         pink       purple          red    royalblue 
    ##          203           88          214          174          484           62 
    ##       salmon          tan    turquoise       yellow 
    ##           92          135         4338          748

``` r
utils::capture.output(length(unique(muscle_moduleColors)),table(muscle_moduleColors), file = "muscle_modules.txt")
```

## Liver

``` r
# Expression and trait data 
load(file = "Liver-dataInput.RData");

# Load network data 
load(file = "Liver-networkConstruction-auto.RData"); # default parameters: deepsplit=2; mergeCutHeight=0.15

# Check the sample order
all(rownames(liver_traitData) %in% rownames(liverDataExpr))
```

    ## [1] TRUE

``` r
all(rownames(liver_traitData) == rownames(liverDataExpr))
```

    ## [1] TRUE

``` r
# Convert Age and Sex variables into factors
liver_traitData$Age <- factor(liver_traitData$Age) # 3,6,9,12,15,18,21,24,27
liver_traitData$Sex <- factor(liver_traitData$Sex) # Female == 0; Male == 1

# Order the expression data frame by increasing time point
liver_traitData <- liver_traitData[order(liver_traitData$Age, decreasing=FALSE),]
liverDataExpr <- liverDataExpr[rownames(liver_traitData),]

# Check the number and size of the identified network modules
length(unique(liver_moduleColors))
```

    ## [1] 37

``` r
table(liver_moduleColors)
```

    ## liver_moduleColors
    ##          black           blue          brown           cyan      darkgreen 
    ##            482            954            791            223            122 
    ##       darkgrey    darkmagenta darkolivegreen     darkorange        darkred 
    ##            116             66             70             97            123 
    ##  darkturquoise          green    greenyellow           grey         grey60 
    ##            118            595            293           6603            182 
    ##      lightcyan     lightgreen    lightyellow        magenta   midnightblue 
    ##            192            178            158            416            217 
    ##         orange  paleturquoise           pink         purple            red 
    ##            110             71            449            346            536 
    ##      royalblue    saddlebrown         salmon        sienna3        skyblue 
    ##            146             87            232             64             87 
    ##      steelblue            tan      turquoise         violet          white 
    ##             82            267           4772             70             93 
    ##         yellow    yellowgreen 
    ##            687             62

``` r
utils::capture.output(length(unique(liver_moduleColors)),table(liver_moduleColors), file = "liver_modules.txt")
```

## Pancreas

``` r
# Expression and trait data 
load(file = "Pancreas-dataInput.RData");

# Load network data 
load(file = "Pancreas-networkConstruction-auto.RData"); # default parameters: deepsplit=2; mergeCutHeight=0.15

# Check the sample order
all(rownames(pancreas_traitData) %in% rownames(pancreasDataExpr))
```

    ## [1] TRUE

``` r
all(rownames(pancreas_traitData) == rownames(pancreasDataExpr))
```

    ## [1] TRUE

``` r
# Convert Age and Sex variables into factors
pancreas_traitData$Age <- factor(pancreas_traitData$Age) # 3,6,9,12,15,18,21,24,27
pancreas_traitData$Sex <- factor(pancreas_traitData$Sex) # Female == 0; Male == 1

# Order the expression data frame by increasing time point
pancreas_traitData <- pancreas_traitData[order(pancreas_traitData$Age, decreasing=FALSE),]
pancreasDataExpr <- pancreasDataExpr[rownames(pancreas_traitData),]

# Check the number and size of the identified network modules
length(unique(pancreas_moduleColors))
```

    ## [1] 11

``` r
table(pancreas_moduleColors)
```

    ## pancreas_moduleColors
    ##     black      blue     brown     green      grey   magenta      pink    purple 
    ##       141      1596       981       496      9162        94       129        55 
    ##       red turquoise    yellow 
    ##       475      4686       596

``` r
utils::capture.output(length(unique(pancreas_moduleColors)),table(pancreas_moduleColors), file = "pancreas_modules.txt")
```

# 3. Relating modules to external information and identifying important genes

## 3.1 Calculation of module-trait associations

When constructing the network, for each identified module a summary expression profile (eigengene) was calculated, corresponding to the *"first principal component of the expression matrix of the corresponding module."*

As a first step to identify modules of interest, we correlate each module's eigengene with both traits.

### Brain

``` r
# Recalculate MEs with color labels
brain_MEs = moduleEigengenes(brainDataExpr, brain_moduleColors, excludeGrey = TRUE)$eigengenes
brain_MEs = orderMEs(brain_MEs)

brain_MEcors <- bicorAndPvalue(brain_MEs, brain_traitData[3:4], maxPOutliers=0.1, robustY=F) #maxPOutliers and robustY parameters were set to 0.1 and F, respectively, as suggested in the FAQs
brain_moduleTraitCor = as.data.frame(brain_MEcors$bicor)
names(brain_moduleTraitCor) <- c("Age.bicor","Sex.bicor")
brain_moduleTraitPvalue = as.data.frame(brain_MEcors$p)
names(brain_moduleTraitPvalue) <- c("Age.pval","Sex.pval")

brain_moduleTrait = cbind(brain_moduleTraitCor,brain_moduleTraitPvalue)

brain_moduleTrait$Age.padj = p.adjust(brain_moduleTrait$Age.pval, method="fdr")
brain_moduleTrait$Sex.padj = p.adjust(brain_moduleTrait$Sex.pval, method="fdr")

write.table(brain_moduleTrait, file="brain_moduleTrait.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Heart

``` r
# Recalculate MEs with color labels
heart_MEs = moduleEigengenes(heartDataExpr, heart_moduleColors, excludeGrey = TRUE)$eigengenes
heart_MEs = orderMEs(heart_MEs)

heart_MEcors <- bicorAndPvalue(heart_MEs, heart_traitData[3:4], maxPOutliers=0.1, robustY=F) #maxPOutliers and robustY parameters were set to 0.1 and F, respectively, as suggested in the FAQs
heart_moduleTraitCor = as.data.frame(heart_MEcors$bicor)
names(heart_moduleTraitCor) <- c("Age.bicor","Sex.bicor")
heart_moduleTraitPvalue = as.data.frame(heart_MEcors$p)
names(heart_moduleTraitPvalue) <- c("Age.pval","Sex.pval")

heart_moduleTrait = cbind(heart_moduleTraitCor,heart_moduleTraitPvalue)

heart_moduleTrait$Age.padj = p.adjust(heart_moduleTrait$Age.pval, method="fdr")
heart_moduleTrait$Sex.padj = p.adjust(heart_moduleTrait$Sex.pval, method="fdr")

write.table(heart_moduleTrait, file="heart_moduleTrait.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Muscle

``` r
# Recalculate MEs with color labels
muscle_MEs = moduleEigengenes(muscleDataExpr, muscle_moduleColors, excludeGrey = TRUE)$eigengenes
muscle_MEs = orderMEs(muscle_MEs)

muscle_MEcors <- bicorAndPvalue(muscle_MEs, muscle_traitData[3:4], maxPOutliers=0.1, robustY=F) #maxPOutliers and robustY parameters were set to 0.1 and F, respectively, as suggested in the FAQs
muscle_moduleTraitCor = as.data.frame(muscle_MEcors$bicor)
names(muscle_moduleTraitCor) <- c("Age.bicor","Sex.bicor")
muscle_moduleTraitPvalue = as.data.frame(muscle_MEcors$p)
names(muscle_moduleTraitPvalue) <- c("Age.pval","Sex.pval")

muscle_moduleTrait = cbind(muscle_moduleTraitCor,muscle_moduleTraitPvalue)

muscle_moduleTrait$Age.padj = p.adjust(muscle_moduleTrait$Age.pval, method="fdr")
muscle_moduleTrait$Sex.padj = p.adjust(muscle_moduleTrait$Sex.pval, method="fdr")

write.table(muscle_moduleTrait, file="muscle_moduleTrait.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Liver

``` r
# Recalculate MEs with color labels
liver_MEs = moduleEigengenes(liverDataExpr, liver_moduleColors, excludeGrey = TRUE)$eigengenes
liver_MEs = orderMEs(liver_MEs)

liver_MEcors <- bicorAndPvalue(liver_MEs, liver_traitData[3:4], maxPOutliers=0.1, robustY=F) #maxPOutliers and robustY parameters were set to 0.1 and F, respectively, as suggested in the FAQs
liver_moduleTraitCor = as.data.frame(liver_MEcors$bicor)
names(liver_moduleTraitCor) <- c("Age.bicor","Sex.bicor")
liver_moduleTraitPvalue = as.data.frame(liver_MEcors$p)
names(liver_moduleTraitPvalue) <- c("Age.pval","Sex.pval")

liver_moduleTrait = cbind(liver_moduleTraitCor,liver_moduleTraitPvalue)

liver_moduleTrait$Age.padj = p.adjust(liver_moduleTrait$Age.pval, method="fdr")
liver_moduleTrait$Sex.padj = p.adjust(liver_moduleTrait$Sex.pval, method="fdr")

write.table(liver_moduleTrait, file="liver_moduleTrait.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Pancreas

``` r
# Recalculate MEs with color labels
pancreas_MEs = moduleEigengenes(pancreasDataExpr, pancreas_moduleColors, excludeGrey = TRUE)$eigengenes
pancreas_MEs = orderMEs(pancreas_MEs)

pancreas_MEcors <- bicorAndPvalue(pancreas_MEs, pancreas_traitData[3:4], maxPOutliers=0.1, robustY=F) #maxPOutliers and robustY parameters were set to 0.1 and F, respectively, as suggested in the FAQs
pancreas_moduleTraitCor = as.data.frame(pancreas_MEcors$bicor)
names(pancreas_moduleTraitCor) <- c("Age.bicor","Sex.bicor")
pancreas_moduleTraitPvalue = as.data.frame(pancreas_MEcors$p)
names(pancreas_moduleTraitPvalue) <- c("Age.pval","Sex.pval")

pancreas_moduleTrait = cbind(pancreas_moduleTraitCor,pancreas_moduleTraitPvalue)

pancreas_moduleTrait$Age.padj = p.adjust(pancreas_moduleTrait$Age.pval, method="fdr")
pancreas_moduleTrait$Sex.padj = p.adjust(pancreas_moduleTrait$Sex.pval, method="fdr")

write.table(pancreas_moduleTrait, file="pancreas_moduleTrait.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Fig 2A

Next we plot the calculated correlations in an heatmap where positive and negative correlations are depicted in red and blue, respectively.

``` r
# prepare columns of data frame to input to ggplot 
brain_Modules <- rownames(brain_moduleTrait)
brain_Tissue <- rep("Brain",nrow(brain_moduleTrait))
brain_Age.bicor <- brain_moduleTrait$Age.bicor
brain_Age.padj <- brain_moduleTrait$Age.padj
brain_Sex.bicor <- brain_moduleTrait$Sex.bicor
brain_Sex.padj <- brain_moduleTrait$Sex.padj

heart_Modules <- rownames(heart_moduleTrait)
heart_Tissue <- rep("Heart",nrow(heart_moduleTrait))
heart_Age.bicor <- heart_moduleTrait$Age.bicor
heart_Age.padj <- heart_moduleTrait$Age.padj
heart_Sex.bicor <- heart_moduleTrait$Sex.bicor
heart_Sex.padj <- heart_moduleTrait$Sex.padj

muscle_Modules <- rownames(muscle_moduleTrait)
muscle_Tissue <- rep("Muscle",nrow(muscle_moduleTrait))
muscle_Age.bicor <- muscle_moduleTrait$Age.bicor
muscle_Age.padj <- muscle_moduleTrait$Age.padj
muscle_Sex.bicor <- muscle_moduleTrait$Sex.bicor
muscle_Sex.padj <- muscle_moduleTrait$Sex.padj

liver_Modules <- rownames(liver_moduleTrait)
liver_Tissue <- rep("Liver",nrow(liver_moduleTrait))
liver_Age.bicor <- liver_moduleTrait$Age.bicor
liver_Age.padj <- liver_moduleTrait$Age.padj
liver_Sex.bicor <- liver_moduleTrait$Sex.bicor
liver_Sex.padj <- liver_moduleTrait$Sex.padj

pancreas_Modules <- rownames(pancreas_moduleTrait)
pancreas_Tissue <- rep("Pancreas",nrow(pancreas_moduleTrait))
pancreas_Age.bicor <- pancreas_moduleTrait$Age.bicor
pancreas_Age.padj <- pancreas_moduleTrait$Age.padj
pancreas_Sex.bicor <- pancreas_moduleTrait$Sex.bicor
pancreas_Sex.padj <- pancreas_moduleTrait$Sex.padj

Modules <- c(brain_Modules,heart_Modules,muscle_Modules,liver_Modules,pancreas_Modules)
Modules <- sapply(strsplit(Modules, "\\ME"), paste, collapse = "ME-") # to improve readability of the modules' names
Tissue <- c(brain_Tissue,heart_Tissue,muscle_Tissue,liver_Tissue,pancreas_Tissue)

# create data frame for the age variable
Bicor <- c(brain_Age.bicor,heart_Age.bicor,muscle_Age.bicor,liver_Age.bicor,pancreas_Age.bicor)
pAdj <- c(brain_Age.padj,heart_Age.padj,muscle_Age.padj,liver_Age.padj,pancreas_Age.padj)
Label <- paste(signif(Bicor, 2), " (", signif(pAdj, 1),")", sep = "")
Trait <- rep("Age",nrow(brain_moduleTrait)+nrow(heart_moduleTrait)+nrow(muscle_moduleTrait)+nrow(liver_moduleTrait)+nrow(pancreas_moduleTrait))

data1 <- data.frame(Trait,Modules,Bicor,pAdj,Tissue,Label)

# create data frame for the sex variable
Bicor <- c(brain_Sex.bicor,heart_Sex.bicor,muscle_Sex.bicor,liver_Sex.bicor,pancreas_Sex.bicor)
pAdj <- c(brain_Sex.padj,heart_Sex.padj,muscle_Sex.padj,liver_Sex.padj,pancreas_Sex.padj)
Label <- paste(signif(Bicor, 2), " (", signif(pAdj, 1),")", sep = "")
Trait <- rep("Sex",nrow(brain_moduleTrait)+nrow(heart_moduleTrait)+nrow(muscle_moduleTrait)+nrow(liver_moduleTrait)+nrow(pancreas_moduleTrait))

data2 <- data.frame(Trait,Modules,Bicor,pAdj,Tissue,Label)

# merge variable data frames to input to ggplot 
data <- rbind(data1, data2)

# add * or º to highlight modules significantly associated with age or sex, respectively (bicor > 0.5 and FDR < 0.05)
data$Label <- ifelse(data$Trait=="Age" & data$pAdj<=0.05 & abs(signif(data$Bicor,2))>=0.5, paste(data$Label, " *", sep=""), ifelse(data$Trait=="Sex" & data$pAdj<=0.05 & abs(signif(data$Bicor,2))>=0.5, paste(data$Label, " º", sep=""), data$Label))

# Plot heatmap with all modules
p1 <- ggplot(data, aes(Trait, Modules, fill= Bicor)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white",  high = "red", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", breaks = c(-1,-0.5,0,0.5,1), limits=c(-1,1)) +
  labs(fill="Correlation", x="\n \nTrait", y="Module Eigengene") +
  geom_text(aes(label=Label), size=5.5) +
  theme_light() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=25,face="bold"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title =element_text(size=20,face="bold"), legend.text = element_text(size=20), axis.text.x = element_blank(), strip.text = element_text(size=25, face="bold"), legend.position="right", rect = element_rect(fill = "transparent")) +
  facet_grid(cols=vars(Tissue))

p1
```

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
g1 <- ggplot_gtable(ggplot_build(p1))
stript <- which(grepl('strip-t', g1$layout$name))
fills <- c("#77AADD","#EE8866","#EEDD88","#FFAABB","#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff('Fig-2A.tiff', units="cm", width=60, height=35, res=300, compression = 'lzw')
grid.draw(g1)
grid.text(label = "Age", rot=0, x = unit(0.18, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.25, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.34, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.41, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.49, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.56, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.65, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.72, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.81, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.88, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
dev.off()
```

    ## png 
    ##   2

## 3.2 Significant module-age (and module-sex) associations

``` r
# change colors to a more color-blindness friendly palette (blue-white-red instead of green-black-red)
redblue<-colorRampPalette(c("blue","white","red"))

#rebuild the original function to support new color palette                     
plotMat_blue<-function(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="", ...)
{
#  X <-x
  n<-nrow(x)
  p<-ncol(x)      

  image(1:p,1:n,t(x[n:1,]),col=redblue(nrgcols),axes=FALSE, xlab="", ylab="", ... ) 

  if(length(ccols)==1){
    axis(3,at=1:p,labels=clabels,las=2,cex.axis=0.6,col.axis=ccols)
      }

  if(length(ccols)==p){
    cols<-unique(ccols)
    for(i in 1:length(cols)){
      which<-(1:p)[ccols==cols[i]]
      axis(3,at=which,labels=clabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  if(length(rcols)==1){
    axis(2,at=n:1,labels=rlabels,las=2,cex.axis=0.6,col.axis=rcols)
      }

  if(length(rcols)==n){
    cols<-unique(rcols)
    for(i in 1:length(cols)){
      which<-(1:n)[rcols==cols[i]]
      axis(2,at=(n:1)[which],labels=rlabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  mtext(title,side=3,line=3)
  box()
}
```

### Fig 3A/Fig S2

### Brain

``` r
# attribute fixed colors to the different time points
colors <- rep(c('#e6194B','#f58231','#ffe119','#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6'))

colors_sex <- rep(c('lightcoral','cornflowerblue'))

# TAN module - AGE

sizeGrWindow(8,9);
brain_which.module="tan"
brain_ME=brain_MEs[, paste("ME",brain_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(brainDataExpr[,brain_moduleColors==brain_which.module ])),
        nrgcols=30,rlabels=F,rcols=brain_which.module,
        main=expression(paste("Brain Tan  ", "bicor=0.86; ", italic("p-value"), "=2e-13")), cex.main=1.5)

legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(brain_ME, col=colors[brain_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(brain_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)


# GREY60 module - GENDER
sizeGrWindow(8,9);
brain_which.module="grey60"
brain_ME=brain_MEs[, paste("ME",brain_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(brainDataExpr[,brain_moduleColors==brain_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=brain_which.module,
        main=expression(paste("Brain Grey60  ", "bicor=0.84; ", italic("p-value"), "=9e-12")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(brain_ME, col=colors_sex[brain_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)


# save plots
tiff('Fig-3A-heat-brain-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
brain_which.module="tan"
brain_ME=brain_MEs[, paste("ME",brain_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(brainDataExpr[,brain_moduleColors==brain_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=brain_which.module,
        main=expression(paste("Brain Tan  ", "bicor=0.86; ", italic("p-value"), "=2e-13")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-brain-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(brain_ME, col=colors[brain_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(brain_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-brain-grey60.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
brain_which.module="grey60"
brain_ME=brain_MEs[, paste("ME",brain_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(brainDataExpr[,brain_moduleColors==brain_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=brain_which.module,
        main=expression(paste("Brain Grey60  ", "bicor=0.84; ", italic("p-value"), "=9e-12")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-brain-grey60.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(brain_ME, col=colors_sex[brain_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

### Heart

``` r
# attribute fixed colors to the different time points
colors <- rep(c('#e6194B','#f58231','#ffe119','#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6'))

# TAN module - AGE
sizeGrWindow(8,9);
heart_which.module="tan"
heart_ME=heart_MEs[, paste("ME",heart_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(heartDataExpr[,heart_moduleColors==heart_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=heart_which.module,
        main=expression(paste("Heart Tan  ", "bicor=0.75; ", italic("p-value"), "=6e-08")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(heart_ME, col=colors[heart_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(heart_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# BLUE module - AGE
sizeGrWindow(8,9);
heart_which.module="blue"
heart_ME=heart_MEs[, paste("ME",heart_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(heartDataExpr[,heart_moduleColors==heart_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=heart_which.module,
        main=expression(paste("Heart Blue  ", "bicor=-0.57; ", italic("p-value"), "=5e-04")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(heart_ME, col=colors[heart_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(heart_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)


# save plots
tiff('Fig-3A-heat-heart-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
heart_which.module="tan"
heart_ME=heart_MEs[, paste("ME",heart_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(heartDataExpr[,heart_moduleColors==heart_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=heart_which.module,
        main=expression(paste("Heart Tan  ", "bicor=0.75; ", italic("p-value"), "=6e-08")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-heart-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(heart_ME, col=colors[heart_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(heart_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-heat-heart-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
heart_which.module="blue"
heart_ME=heart_MEs[, paste("ME",heart_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(heartDataExpr[,heart_moduleColors==heart_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=heart_which.module,
        main=expression(paste("Heart Blue  ", "bicor=-0.57; ", italic("p-value"), "=5e-04")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-heart-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(heart_ME, col=colors[heart_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(heart_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

### Muscle

``` r
# attribute fixed colors to the different time points
colors <- rep(c('#e6194B','#f58231','#ffe119','#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6'))

colors_sex <- c("lightcoral","cornflowerblue")

# MAGENTA module - AGE
sizeGrWindow(8,9);
muscle_which.module="magenta"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Magenta  ", "bicor=0.73; ", italic("p-value"), "=3e-07")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors[muscle_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(muscle_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# BROWN module - AGE
sizeGrWindow(8,9);
muscle_which.module="brown"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Brown  ", "bicor=-0.76; ", italic("p-value"), "=8e-08")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors[muscle_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(muscle_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# RED module - GENDER
sizeGrWindow(8,9);
muscle_which.module="red"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Red  ", "bicor=0.86; ", italic("p-value"), "=3e-12")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# PURPLE module - GENDER
sizeGrWindow(8,9);
muscle_which.module="purple"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Purple  ", "bicor=-0.76; ", italic("p-value"), "=2e-08")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# GREENYELLOW module - GENDER
sizeGrWindow(8,9);
muscle_which.module="greenyellow"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Greenyellow  ", "bicor=-0.85; ", italic("p-value"), "=4e-12")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# BLUE module - GENDER
sizeGrWindow(8,9);
muscle_which.module="blue"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Blue  ", "bicor=0.5; ", italic("p-value"), "=0.004")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# save plots
tiff('Fig-3A-heat-muscle-magenta.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="magenta"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Magenta  ", "bicor=0.73; ", italic("p-value"), "=3e-07")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-muscle-magenta.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors[muscle_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(muscle_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-heat-muscle-brown.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="brown"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Brown  ", "bicor=-0.76; ", italic("p-value"), "=8e-08")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-muscle-brown.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors[muscle_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(muscle_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-muscle-red.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="red"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ), nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Red  ", "bicor=0.86; ", italic("p-value"), "=3e-12")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-muscle-red.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-muscle-purple.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="purple"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Purple  ", "bicor=-0.76; ", italic("p-value"), "=2e-08")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-muscle-purple.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-muscle-greenyellow.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="greenyellow"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Greenyellow ", "bicor=-0.85; ", italic("p-value"), "=4e-12")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-muscle-greenyellow.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-muscle-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
muscle_which.module="blue"
muscle_ME=muscle_MEs[, paste("ME",muscle_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(muscleDataExpr[,muscle_moduleColors==muscle_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=muscle_which.module,
        main=expression(paste("Muscle Blue  ", "bicor=0.5; ", italic("p-value"), "=0.004")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplott-muscle-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(muscle_ME, col=colors_sex[muscle_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

### Liver

``` r
# attribute fixed colors to the different time points
colors <- rep(c('#e6194B','#f58231','#ffe119','#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6'))

colors_sex <- c("lightcoral","cornflowerblue")

# SALMON module - AGE
sizeGrWindow(8,9);
liver_which.module="salmon"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Salmon  ", "bicor=0.71; ", italic("p-value"), "=1e-06")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors[liver_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(liver_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# DARKTURQUOISE module - AGE
sizeGrWindow(8,9);
liver_which.module="darkturquoise"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkturquoise  ", "bicor=0.54; ", italic("p-value"), "=0.001")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors[liver_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(liver_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# DARKGREY module - AGE | GENDER
sizeGrWindow(8,9);
liver_which.module="darkgrey"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkgrey  ", "bicor=0.79; ", italic("p-value"), "=1e-09")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(liver_traitData$Sex), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)

# TAN module - GENDER
sizeGrWindow(8,9);
liver_which.module="tan"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Tan  ", "bicor=0.55; ", italic("p-value"), "=6e-04")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# RED module - GENDER
sizeGrWindow(8,9);
liver_which.module="red"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Red  ", "bicor=-0.68; ", italic("p-value"), "=2e-06")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# DARKOLIVEGREEN module - GENDER
sizeGrWindow(8,9);
liver_which.module="darkolivegreen"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkolivegreen  ", "bicor=-0.53; ", italic("p-value"), "=7e-04")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# CYAN module - GENDER
sizeGrWindow(8,9);
liver_which.module="cyan"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Cyan  ", "bicor=0.55; ", italic("p-value"), "=6e-04")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# BROWN module - GENDER
sizeGrWindow(8,9);
liver_which.module="brown"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Brown  ", "bicor=-0.89; ", italic("p-value"), "=1e-14")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# BLUE module - GENDER
sizeGrWindow(8,9);
liver_which.module="blue"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]

layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Blue  ", "bicor=0.84; ", italic("p-value"), "=1e-11")), cex.main=1.5)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)

# save plots
tiff('Fig-3A-heat-liver-salmon.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="salmon"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Salmon  ", "bicor=0.71; ", italic("p-value"), "=1e-06")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-liver-salmon.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors[liver_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(liver_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-heat-liver-darkturquoise.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="darkturquoise"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkturquoise  ", "bicor=0.54; ", italic("p-value"), "=0.001")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-3A-barplot-liver-darkturquoise.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors[liver_traitData$Age], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(0,-0.2),legend=levels(liver_traitData$Age), pch=15, col=colors, bty="n", ncol=9, cex = 2, xpd=NA, pt.cex=2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-darkgrey.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="darkgrey"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkgrey  ", "bicor=0.79; ", italic("p-value"), "=1e-09")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-darkgrey.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="tan"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Tan  ", "bicor=0.55; ", italic("p-value"), "=6e-04")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-tan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-red.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="red"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Red  ", "bicor=-0.68; ", italic("p-value"), "=2e-06")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-red.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-darkolivegreen.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="darkolivegreen"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Darkolivegreen ", "bicor=-0.53; ", italic("p-value"), "=7e-4")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-darkolivegreen.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-cyan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="cyan"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Cyan  ", "bicor=0.55; ", italic("p-value"), "=6e-04")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-cyan.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-brown.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="brown"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Brown  ", "bicor=-0.89; ", italic("p-value"), "=1e-14")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-brown.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-heat-liver-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
liver_which.module="blue"
liver_ME=liver_MEs[, paste("ME",liver_which.module, sep="")]
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plotMat_blue(t(scale(liverDataExpr[,liver_moduleColors==liver_which.module ]) ),
        nrgcols=30,rlabels=F,rcols=liver_which.module,
        main=expression(paste("Liver Blue  ", "bicor=0.84; ", italic("p-value"), "=1e-11")), cex.main=2)
legend_image <- as.raster(matrix(rev(redblue(30)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-Score', cex.main=1.5)
text(x=1.5, y = seq(0,1,l=5), labels = c(-4,-2,0,2,4), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()
```

    ## pdf 
    ##   2

``` r
tiff('Fig-S2-barplot-liver-blue.tiff', units="cm", width=25, height=15, res=300, compression = 'lzw')
barplot(liver_ME, col=colors_sex[liver_traitData$Sex], main="", cex.lab=2,
        ylab="eigengene expression", names.arg = F)
legend("bottom",inset=c(-0.05,-0.15),legend=c("Females","Males"), pch=15, col=colors_sex, bty="n", ncol=2, cex = 1.2, xpd=NA, pt.cex=1.2)
dev.off()
```

    ## pdf 
    ##   2

### Pancreas

There are no modules significantly associated with age nor with sex.

## 3.3 Gene relationship to trait and important modules: Gene Significance and Module Membership

So far we have evaluated the association of each module's eigengene with age and sex and than selected all modules that exhibited a significant (FDR &lt; 0.05) correlation with age equal or higher than 0.5.

Next, we proceed our analysis at the level of individual genes by calculating the correlation of each gene's expression with age (Gene Significance - GS), as well as the correlation of each gene's expression with the eigengene's expression (Module Membership - MM).

### Brain

``` r
# Define variable Age containing the age column of the trait data frame
brain_age = as.data.frame(brain_traitData$Age)
rownames(brain_age) <- rownames(brain_traitData)
names(brain_age) = "Age"

brain_sex = as.data.frame(brain_traitData$Sex)
rownames(brain_sex) <- rownames(brain_traitData)
names(brain_sex) = "Sex"

# Extract the names of the modules
brain_modNames = substring(names(brain_MEs), 3)

# Calculate module membership
brain_MMs <- bicorAndPvalue(brainDataExpr,brain_MEs,maxPOutliers=0.1, robustY=F)

brain_geneModuleMembership = as.data.frame(brain_MMs$bicor)
names(brain_geneModuleMembership) = paste("MM.", brain_modNames, sep="")
brain_MMPvalue = as.data.frame(brain_MMs$p)
names(brain_MMPvalue) = paste("p.MM.", brain_modNames, sep="")

brain_MMPadj = sapply(brain_MMPvalue, p.adjust, method="fdr")
rownames(brain_MMPadj) <- rownames(brain_MMPvalue)
colnames(brain_MMPadj) = paste("pAdj.MM.", brain_modNames, sep="")

brain_moduleMembership = cbind(brain_geneModuleMembership,brain_MMPvalue,brain_MMPadj)

write.table(brain_moduleMembership, file="brain_moduleMembership.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

# Calculate gene significance
#age
brain_GSs_age <- bicorAndPvalue(brainDataExpr,brain_age,maxPOutliers=0.1, robustY=F)

brain_geneTraitSignificance_age = as.data.frame(brain_GSs_age$bicor)
names(brain_geneTraitSignificance_age) = paste("GS.", names(brain_age), sep="")
brain_GSPvalue_age = as.data.frame(brain_GSs_age$p)
names(brain_GSPvalue_age) = paste("p.GS.", names(brain_age), sep="")

brain_GSPadj_age <- as.data.frame(p.adjust(brain_GSPvalue_age$p.GS.Age, method="fdr"))
rownames(brain_GSPadj_age) <- rownames(brain_GSPvalue_age)
colnames(brain_GSPadj_age) = paste("pAdj.GS.", names(brain_age), sep="")

#sex
brain_GSs_sex <- bicorAndPvalue(brainDataExpr,brain_sex,maxPOutliers=0.1, robustY=F)

brain_geneTraitSignificance_sex = as.data.frame(brain_GSs_sex$bicor)
names(brain_geneTraitSignificance_sex) = paste("GS.", names(brain_sex), sep="")
brain_GSPvalue_sex = as.data.frame(brain_GSs_sex$p)
names(brain_GSPvalue_sex) = paste("p.GS.", names(brain_sex), sep="")

brain_GSPadj_sex <- as.data.frame(p.adjust(brain_GSPvalue_sex$p.GS.Sex, method="fdr"))
rownames(brain_GSPadj_sex) <- rownames(brain_GSPvalue_sex)
colnames(brain_GSPadj_sex) = paste("pAdj.GS.", names(brain_sex), sep="")

brain_geneSignificance = cbind(brain_geneTraitSignificance_age,brain_GSPvalue_age,brain_GSPadj_age,brain_geneTraitSignificance_sex,brain_GSPvalue_sex,brain_GSPadj_sex)

write.table(brain_geneSignificance, file="brain_geneSignificance.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Heart

``` r
# Define variable Age containing the age column of the trait data frame
heart_age = as.data.frame(heart_traitData$Age)
rownames(heart_age) <- rownames(heart_traitData)
names(heart_age) = "Age"

heart_sex = as.data.frame(heart_traitData$Sex)
rownames(heart_sex) <- rownames(heart_traitData)
names(heart_sex) = "Sex"

# Extract the names of the modules
heart_modNames = substring(names(heart_MEs), 3)

# Calculate module membership
heart_MMs <- bicorAndPvalue(heartDataExpr,heart_MEs,maxPOutliers=0.1, robustY=F)

heart_geneModuleMembership = as.data.frame(heart_MMs$bicor)
names(heart_geneModuleMembership) = paste("MM.", heart_modNames, sep="")
heart_MMPvalue = as.data.frame(heart_MMs$p)
names(heart_MMPvalue) = paste("p.MM.", heart_modNames, sep="")

heart_MMPadj = sapply(heart_MMPvalue, p.adjust, method="fdr")
rownames(heart_MMPadj) <- rownames(heart_MMPvalue)
colnames(heart_MMPadj) = paste("pAdj.MM.", heart_modNames, sep="")

heart_moduleMembership = cbind(heart_geneModuleMembership,heart_MMPvalue,heart_MMPadj)

write.table(heart_moduleMembership, file="heart_moduleMembership.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

# Calculate gene significance
#age
heart_GSs_age <- bicorAndPvalue(heartDataExpr,heart_age,maxPOutliers=0.1, robustY=F)

heart_geneTraitSignificance_age = as.data.frame(heart_GSs_age$bicor)
names(heart_geneTraitSignificance_age) = paste("GS.", names(heart_age), sep="")
heart_GSPvalue_age = as.data.frame(heart_GSs_age$p)
names(heart_GSPvalue_age) = paste("p.GS.", names(heart_age), sep="")

heart_GSPadj_age <- as.data.frame(p.adjust(heart_GSPvalue_age$p.GS.Age, method="fdr"))
rownames(heart_GSPadj_age) <- rownames(heart_GSPvalue_age)
colnames(heart_GSPadj_age) = paste("pAdj.GS.", names(heart_age), sep="")

#sex
heart_GSs_sex <- bicorAndPvalue(heartDataExpr,heart_sex,maxPOutliers=0.1, robustY=F)

heart_geneTraitSignificance_sex = as.data.frame(heart_GSs_sex$bicor)
names(heart_geneTraitSignificance_sex) = paste("GS.", names(heart_sex), sep="")
heart_GSPvalue_sex = as.data.frame(heart_GSs_sex$p)
names(heart_GSPvalue_sex) = paste("p.GS.", names(heart_sex), sep="")

heart_GSPadj_sex <- as.data.frame(p.adjust(heart_GSPvalue_sex$p.GS.Sex, method="fdr"))
rownames(heart_GSPadj_sex) <- rownames(heart_GSPvalue_sex)
colnames(heart_GSPadj_sex) = paste("pAdj.GS.", names(heart_sex), sep="")

heart_geneSignificance = cbind(heart_geneTraitSignificance_age,heart_GSPvalue_age,heart_GSPadj_age,heart_geneTraitSignificance_sex,heart_GSPvalue_sex,heart_GSPadj_sex)

write.table(heart_geneSignificance, file="heart_geneSignificance.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Muscle

``` r
# Define variable Age containing the age column of the trait data frame
muscle_age = as.data.frame(muscle_traitData$Age)
rownames(muscle_age) <- rownames(muscle_traitData)
names(muscle_age) = "Age"

muscle_sex = as.data.frame(muscle_traitData$Sex)
rownames(muscle_sex) <- rownames(muscle_traitData)
names(muscle_sex) = "Sex"

# Extract the names of the modules
muscle_modNames = substring(names(muscle_MEs), 3)

# Calculate module membership
muscle_MMs <- bicorAndPvalue(muscleDataExpr,muscle_MEs,maxPOutliers=0.1, robustY=F)

muscle_geneModuleMembership = as.data.frame(muscle_MMs$bicor)
names(muscle_geneModuleMembership) = paste("MM.", muscle_modNames, sep="")
muscle_MMPvalue = as.data.frame(muscle_MMs$p)
names(muscle_MMPvalue) = paste("p.MM.", muscle_modNames, sep="")

muscle_MMPadj = sapply(muscle_MMPvalue, p.adjust, method="fdr")
rownames(muscle_MMPadj) <- rownames(muscle_MMPvalue)
colnames(muscle_MMPadj) = paste("pAdj.MM.", muscle_modNames, sep="")

muscle_moduleMembership = cbind(muscle_geneModuleMembership,muscle_MMPvalue,muscle_MMPadj)

write.table(muscle_moduleMembership, file="muscle_moduleMembership.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

# Calculate gene significance
#age
muscle_GSs_age <- bicorAndPvalue(muscleDataExpr,muscle_age,maxPOutliers=0.1, robustY=F)

muscle_geneTraitSignificance_age = as.data.frame(muscle_GSs_age$bicor)
names(muscle_geneTraitSignificance_age) = paste("GS.", names(muscle_age), sep="")
muscle_GSPvalue_age = as.data.frame(muscle_GSs_age$p)
names(muscle_GSPvalue_age) = paste("p.GS.", names(muscle_age), sep="")

muscle_GSPadj_age <- as.data.frame(p.adjust(muscle_GSPvalue_age$p.GS.Age, method="fdr"))
rownames(muscle_GSPadj_age) <- rownames(muscle_GSPvalue_age)
colnames(muscle_GSPadj_age) = paste("pAdj.GS.", names(muscle_age), sep="")

#sex
muscle_GSs_sex <- bicorAndPvalue(muscleDataExpr,muscle_sex,maxPOutliers=0.1, robustY=F)

muscle_geneTraitSignificance_sex = as.data.frame(muscle_GSs_sex$bicor)
names(muscle_geneTraitSignificance_sex) = paste("GS.", names(muscle_sex), sep="")
muscle_GSPvalue_sex = as.data.frame(muscle_GSs_sex$p)
names(muscle_GSPvalue_sex) = paste("p.GS.", names(muscle_sex), sep="")

muscle_GSPadj_sex <- as.data.frame(p.adjust(muscle_GSPvalue_sex$p.GS.Sex, method="fdr"))
rownames(muscle_GSPadj_sex) <- rownames(muscle_GSPvalue_sex)
colnames(muscle_GSPadj_sex) = paste("pAdj.GS.", names(muscle_sex), sep="")

muscle_geneSignificance = cbind(muscle_geneTraitSignificance_age,muscle_GSPvalue_age,muscle_GSPadj_age,muscle_geneTraitSignificance_sex,muscle_GSPvalue_sex,muscle_GSPadj_sex)

write.table(muscle_geneSignificance, file="muscle_geneSignificance.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Liver

``` r
# Define variable Age containing the age column of the trait data frame
liver_age = as.data.frame(liver_traitData$Age)
rownames(liver_age) <- rownames(liver_traitData)
names(liver_age) = "Age"

liver_sex = as.data.frame(liver_traitData$Sex)
rownames(liver_sex) <- rownames(liver_traitData)
names(liver_sex) = "Sex"

# Extract the names of the modules
liver_modNames = substring(names(liver_MEs), 3)

# Calculate module membership
liver_MMs <- bicorAndPvalue(liverDataExpr,liver_MEs,maxPOutliers=0.1, robustY=F)

liver_geneModuleMembership = as.data.frame(liver_MMs$bicor)
names(liver_geneModuleMembership) = paste("MM.", liver_modNames, sep="")
liver_MMPvalue = as.data.frame(liver_MMs$p)
names(liver_MMPvalue) = paste("p.MM.", liver_modNames, sep="")

liver_MMPadj = sapply(liver_MMPvalue, p.adjust, method="fdr")
rownames(liver_MMPadj) <- rownames(liver_MMPvalue)
colnames(liver_MMPadj) = paste("pAdj.MM.", liver_modNames, sep="")

liver_moduleMembership = cbind(liver_geneModuleMembership,liver_MMPvalue,liver_MMPadj)

write.table(liver_moduleMembership, file="liver_moduleMembership.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

# Calculate gene significance
#age
liver_GSs_age <- bicorAndPvalue(liverDataExpr,liver_age,maxPOutliers=0.1, robustY=F)

liver_geneTraitSignificance_age = as.data.frame(liver_GSs_age$bicor)
names(liver_geneTraitSignificance_age) = paste("GS.", names(liver_age), sep="")
liver_GSPvalue_age = as.data.frame(liver_GSs_age$p)
names(liver_GSPvalue_age) = paste("p.GS.", names(liver_age), sep="")

liver_GSPadj_age <- as.data.frame(p.adjust(liver_GSPvalue_age$p.GS.Age, method="fdr"))
rownames(liver_GSPadj_age) <- rownames(liver_GSPvalue_age)
colnames(liver_GSPadj_age) = paste("pAdj.GS.", names(liver_age), sep="")

#sex
liver_GSs_sex <- bicorAndPvalue(liverDataExpr,liver_sex,maxPOutliers=0.1, robustY=F)

liver_geneTraitSignificance_sex = as.data.frame(liver_GSs_sex$bicor)
names(liver_geneTraitSignificance_sex) = paste("GS.", names(liver_sex), sep="")
liver_GSPvalue_sex = as.data.frame(liver_GSs_sex$p)
names(liver_GSPvalue_sex) = paste("p.GS.", names(liver_sex), sep="")

liver_GSPadj_sex <- as.data.frame(p.adjust(liver_GSPvalue_sex$p.GS.Sex, method="fdr"))
rownames(liver_GSPadj_sex) <- rownames(liver_GSPvalue_sex)
colnames(liver_GSPadj_sex) = paste("pAdj.GS.", names(liver_sex), sep="")

liver_geneSignificance = cbind(liver_geneTraitSignificance_age,liver_GSPvalue_age,liver_GSPadj_age,liver_geneTraitSignificance_sex,liver_GSPvalue_sex,liver_GSPadj_sex)

write.table(liver_geneSignificance, file="liver_geneSignificance.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")
```

### Pancreas

There are no modules significantly associated with age nor with sex.

### Fig 2B

Next we plot the calculated gene significances where the darker the red, the higher the significance of the correlation of a given gene with age.

``` r
### BRAIN 

# tan module
brain_module = "tan"
brain_column = match(brain_module, brain_modNames)
brain_moduleGenes = brain_moduleColors==brain_module

# Calculate the correlation between gene significance for each trait and module membership
brain_age_tan <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_age[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

brain_sex_tan <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_sex[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)


# lightyellow module
brain_module = "lightyellow"
brain_column = match(brain_module, brain_modNames)
brain_moduleGenes = brain_moduleColors==brain_module

# Calculate the correlation between gene significance for each trait and module membership
brain_age_lightyellow <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_age[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

brain_sex_lightyellow <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_sex[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# cyan module
brain_module = "cyan"
brain_column = match(brain_module, brain_modNames)
brain_moduleGenes = brain_moduleColors==brain_module

# Calculate the correlation between gene significance for each trait and module membership
brain_age_cyan <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_age[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

brain_sex_cyan <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_sex[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# grey60 module
brain_module = "grey60"
brain_column = match(brain_module, brain_modNames)
brain_moduleGenes = brain_moduleColors==brain_module

# Calculate the correlation between gene significance for each trait and module membership
brain_age_grey60 <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_age[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

brain_sex_grey60 <- bicorAndPvalue(abs(brain_geneModuleMembership[brain_moduleGenes, brain_column]),abs(brain_geneTraitSignificance_sex[brain_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)


# Prepare columns of data frame to input to ggplot
brain_Modules <- c("MM-tan","MM-lightyellow","MM-cyan","MM-grey60")
brain_Tissue <- rep("Brain",4)

brain_age_Bicor <- c(brain_age_tan$bicor,brain_age_lightyellow$bicor,brain_age_cyan$bicor,brain_age_grey60$bicor)
brain_age_pVal <- c(brain_age_tan$p,brain_age_lightyellow$p,brain_age_cyan$p,brain_age_grey60$p)

brain_sex_Bicor <- c(brain_sex_tan$bicor,brain_sex_lightyellow$bicor,brain_sex_cyan$bicor,brain_sex_grey60$bicor)
brain_sex_pVal <- c(brain_sex_tan$p,brain_sex_lightyellow$p,brain_sex_cyan$p,brain_sex_grey60$p)

### Heart 

# tan module
heart_module = "tan"
heart_column = match(heart_module, heart_modNames)
heart_moduleGenes = heart_moduleColors==heart_module

# Calculate the correlation between gene significance for each trait and module membership
heart_age_tan <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_age[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

heart_sex_tan <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_sex[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# cyan module
heart_module = "cyan"
heart_column = match(heart_module, heart_modNames)
heart_moduleGenes = heart_moduleColors==heart_module

# Calculate the correlation between gene significance for each trait and module membership
heart_age_cyan <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_age[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

heart_sex_cyan <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_sex[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# blue module
heart_module = "blue"
heart_column = match(heart_module, heart_modNames)
heart_moduleGenes = heart_moduleColors==heart_module

# Calculate the correlation between gene significance for each trait and module membership
heart_age_blue <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_age[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

heart_sex_blue <- bicorAndPvalue(abs(heart_geneModuleMembership[heart_moduleGenes, heart_column]),abs(heart_geneTraitSignificance_sex[heart_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# Prepare columns of data frame to input to ggplot
heart_Modules <- c("MM-tan","MM-cyan","MM-blue")
heart_Tissue <- rep("Heart",3)

heart_age_Bicor <- c(heart_age_tan$bicor,heart_age_cyan$bicor,heart_age_blue$bicor)
heart_age_pVal <- c(heart_age_tan$p,heart_age_cyan$p,heart_age_blue$p)

heart_sex_Bicor <- c(heart_sex_tan$bicor,heart_sex_cyan$bicor,heart_sex_blue$bicor)
heart_sex_pVal <- c(heart_sex_tan$p,heart_sex_cyan$p,heart_sex_blue$p)

### MUSCLE 

# tan module
muscle_module = "tan"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_tan <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_tan <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# magenta module
muscle_module = "magenta"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_magenta <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_magenta <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# lightyellow module
muscle_module = "lightyellow"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_lightyellow <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_lightyellow <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# darkred module
muscle_module = "darkred"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_darkred <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_darkred <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# brown module
muscle_module = "brown"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_brown <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_brown <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# red module
muscle_module = "red"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_red <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_red <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# purple module
muscle_module = "purple"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_purple <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_purple <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# greenyellow module
muscle_module = "greenyellow"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_greenyellow <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_greenyellow <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# blue module
muscle_module = "blue"
muscle_column = match(muscle_module, muscle_modNames)
muscle_moduleGenes = muscle_moduleColors==muscle_module

# Calculate the correlation between gene significance for each trait and module membership
muscle_age_blue <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_age[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

muscle_sex_blue <- bicorAndPvalue(abs(muscle_geneModuleMembership[muscle_moduleGenes, muscle_column]),abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# Prepare columns of data frame to input to ggplot
muscle_Modules <- c("MM-tan","MM-magenta","MM-lightyellow","MM-darkred","MM-brown","MM-red","MM-purple","MM-greenyellow","MM-blue")
muscle_Tissue <- rep("Muscle",9)

muscle_age_Bicor <- c(muscle_age_tan$bicor,muscle_age_magenta$bicor,muscle_age_lightyellow$bicor,muscle_age_darkred$bicor,muscle_age_brown$bicor,muscle_age_red$bicor,muscle_age_purple$bicor,muscle_age_greenyellow$bicor,muscle_age_blue$bicor)
muscle_age_pVal <- c(muscle_age_tan$p,muscle_age_magenta$p,muscle_age_lightyellow$p,muscle_age_darkred$p,muscle_age_brown$p,muscle_age_red$p,muscle_age_purple$p,muscle_age_greenyellow$p,muscle_age_blue$p)

muscle_sex_Bicor <- c(muscle_sex_tan$bicor,muscle_sex_magenta$bicor,muscle_sex_lightyellow$bicor,muscle_sex_darkred$bicor,muscle_sex_brown$bicor,muscle_sex_red$bicor,muscle_sex_purple$bicor,muscle_sex_greenyellow$bicor,muscle_sex_blue$bicor)
muscle_sex_pVal <- c(muscle_sex_tan$p,muscle_sex_magenta$p,muscle_sex_lightyellow$p,muscle_sex_darkred$p,muscle_sex_brown$p,muscle_sex_red$p,muscle_sex_purple$p,muscle_sex_greenyellow$p,muscle_sex_blue$p)

### LIVER 

# salmon module
liver_module = "salmon"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_salmon <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_salmon <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# lightcyan module
liver_module = "lightcyan"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_lightcyan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_lightcyan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# darkturquoise module
liver_module = "darkturquoise"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_darkturquoise <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_darkturquoise <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# darkgrey module
liver_module = "darkgrey"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_darkgrey <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_darkgrey <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# tan module
liver_module = "tan"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_tan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_tan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# red module
liver_module = "red"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_red <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_red <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# orange module
liver_module = "orange"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_orange <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_orange <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# darkolivegreen module
liver_module = "darkolivegreen"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_darkolivegreen <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_darkolivegreen <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# cyan module
liver_module = "cyan"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_cyan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_cyan <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# brown module
liver_module = "brown"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_brown <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_brown <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# blue module
liver_module = "blue"
liver_column = match(liver_module, liver_modNames)
liver_moduleGenes = liver_moduleColors==liver_module

# Calculate the correlation between gene significance for each trait and module membership
liver_age_blue <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_age[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

liver_sex_blue <- bicorAndPvalue(abs(liver_geneModuleMembership[liver_moduleGenes, liver_column]),abs(liver_geneTraitSignificance_sex[liver_moduleGenes, 1]),maxPOutliers=0.1, robustY=F)

# Prepare columns of data frame to input to ggplot
liver_Modules <- c("MM-salmon","MM-lightcyan","MM-darkturquoise","MM-darkgrey","MM-tan","MM-red","MM-orange","MM-darkolivegreen","MM-cyan","MM-brown","MM-blue")
liver_Tissue <- rep("Liver",11)

liver_age_Bicor <- c(liver_age_salmon$bicor,liver_age_lightcyan$bicor,liver_age_darkturquoise$bicor,liver_age_darkgrey$bicor,liver_age_tan$bicor,liver_age_red$bicor,liver_age_orange$bicor,liver_age_darkolivegreen$bicor,liver_age_cyan$bicor,liver_age_brown$bicor,liver_age_blue$bicor)
liver_age_pVal <- c(liver_age_salmon$p,liver_age_lightcyan$p,liver_age_darkturquoise$p,liver_age_darkgrey$p,liver_age_tan$p,liver_age_red$p,liver_age_orange$p,liver_age_darkolivegreen$p,liver_age_cyan$p,liver_age_brown$p,liver_age_blue$p)

liver_sex_Bicor <- c(liver_sex_salmon$bicor,liver_sex_lightcyan$bicor,liver_sex_darkturquoise$bicor,liver_sex_darkgrey$bicor,liver_sex_tan$bicor,liver_sex_red$bicor,liver_sex_orange$bicor,liver_sex_darkolivegreen$bicor,liver_sex_cyan$bicor,liver_sex_brown$bicor,liver_sex_blue$bicor)
liver_sex_pVal <- c(liver_sex_salmon$p,liver_sex_lightcyan$p,liver_sex_darkturquoise$p,liver_sex_darkgrey$p,liver_sex_tan$p,liver_sex_red$p,liver_sex_orange$p,liver_sex_darkolivegreen$p,liver_sex_cyan$p,liver_sex_brown$p,liver_sex_blue$p)


# create data frame for the age variable
Trait <- rep("Age",27)
Modules <- c(brain_Modules,heart_Modules,muscle_Modules,liver_Modules)
Tissue <- c(brain_Tissue,heart_Tissue,muscle_Tissue,liver_Tissue)
Bicor <- c(brain_age_Bicor,heart_age_Bicor,muscle_age_Bicor,liver_age_Bicor)
pVal <- c(brain_age_pVal,heart_age_pVal,muscle_age_pVal,liver_age_pVal)
Label <- paste(signif(Bicor, 2), " (", signif(pVal, 1),")", sep = "")

data1 <- data.frame(Trait,Modules,Bicor,pVal,Tissue,Label)

# create data frame for the sex variable
Trait <- rep("Sex",27)
Modules <- c(brain_Modules,heart_Modules,muscle_Modules,liver_Modules)
Tissue <- c(brain_Tissue,heart_Tissue,muscle_Tissue,liver_Tissue)
Bicor <- c(brain_sex_Bicor,heart_sex_Bicor,muscle_sex_Bicor,liver_sex_Bicor)
pVal <- c(brain_sex_pVal,heart_sex_pVal,muscle_sex_pVal,liver_sex_pVal)
Label <- paste(signif(Bicor, 2), " (", signif(pVal, 1),")", sep = "")

data2 <- data.frame(Trait,Modules,Bicor,pVal,Tissue,Label)

# merge variable data frames to input to ggplot 
data <- rbind(data1, data2)

# add * or º to highlight modules significantly associated with age or sex, respectively (bicor > 0.5 and FDR < 0.05)
data$Label <- ifelse(data$Trait=="Age" & data$pVal<=0.05 & abs(signif(data$Bicor,2))>=0.5, paste(data$Label, " *", sep=""), ifelse(data$Trait=="Sex" & data$pVal<=0.05 & abs(signif(data$Bicor,2))>=0.5, paste(data$Label, " º", sep=""), data$Label))

# Plot heatmap with all modules
p2 <- ggplot(data, aes(Trait, Modules, fill= Bicor)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white",  high = "red", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", breaks = c(-1,-0.5,0,0.5,1), limits=c(-1,1)) +
  labs(fill="Correlation", x="\n \nTrait", y="Module Eigengene") +
  geom_text(aes(label=Label), size=5.5) +
  theme_light() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=25,face="bold"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title =element_text(size=20,face="bold"), legend.text = element_text(size=20), axis.text.x = element_blank(), strip.text = element_text(size=25, face="bold"), legend.position="right", rect = element_rect(fill = "transparent")) +
  facet_grid(cols=vars(Tissue))

p2
```

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
g2 <- ggplot_gtable(ggplot_build(p2))
stript <- which(grepl('strip-t', g2$layout$name))
fills <- c("#77AADD","#EE8866","#EEDD88","#FFAABB","#44BB99")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff('Fig-2B.tiff', units="cm", width=60, height=35, res=300, compression = 'lzw')
grid.draw(g2)
grid.text(label = "Age", rot=0, x = unit(0.19, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.28, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.39, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.48, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.58, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.67, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Age", rot=0, x = unit(0.78, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
grid.text(label = "Sex", rot=0, x = unit(0.87, "npc"), y = unit(0.095, "npc"),just = "centre", gp=gpar(fontsize=20,fontface="bold"))
dev.off()
```

    ## png 
    ##   2

### Fig 3B /Fig S3

``` r
brain_module1 = "tan"
brain_column1 = match(brain_module1, brain_modNames);
brain_moduleGenes1 = brain_moduleColors==brain_module1;

brain_module2 = "grey60"
brain_column2 = match(brain_module2, brain_modNames);
brain_moduleGenes2 = brain_moduleColors==brain_module2;
 
heart_module1 = "tan"
heart_column1 = match(heart_module1, heart_modNames);
heart_moduleGenes1 = heart_moduleColors==heart_module1;

heart_module2 = "blue"
heart_column2 = match(heart_module2, heart_modNames);
heart_moduleGenes2 = heart_moduleColors==heart_module2;

liver_module1 = "salmon"
liver_column1 = match(liver_module1, liver_modNames);
liver_moduleGenes1 = liver_moduleColors==liver_module1;

liver_module2 = "darkturquoise"
liver_column2 = match(liver_module2, liver_modNames);
liver_moduleGenes2 = liver_moduleColors==liver_module2;

liver_module3 = "tan"
liver_column3 = match(liver_module3, liver_modNames);
liver_moduleGenes3 = liver_moduleColors==liver_module3;

liver_module4 = "red"
liver_column4 = match(liver_module4, liver_modNames);
liver_moduleGenes4 = liver_moduleColors==liver_module4;

liver_module5 = "darkolivegreen"
liver_column5 = match(liver_module5, liver_modNames);
liver_moduleGenes5 = liver_moduleColors==liver_module5;

liver_module6 = "darkgrey"
liver_column6 = match(liver_module6, liver_modNames);
liver_moduleGenes6 = liver_moduleColors==liver_module6;

liver_module7 = "cyan"
liver_column7 = match(liver_module7, liver_modNames);
liver_moduleGenes7 = liver_moduleColors==liver_module7;

liver_module8 = "brown"
liver_column8 = match(liver_module8, liver_modNames);
liver_moduleGenes8 = liver_moduleColors==liver_module8;

liver_module9 = "blue"
liver_column9 = match(liver_module9, liver_modNames);
liver_moduleGenes9 = liver_moduleColors==liver_module9;

muscle_module1 = "magenta"
muscle_column1 = match(muscle_module1, muscle_modNames);
muscle_moduleGenes1 = muscle_moduleColors==muscle_module1;

muscle_module2 = "brown"
muscle_column2 = match(muscle_module2, muscle_modNames);
muscle_moduleGenes2 = muscle_moduleColors==muscle_module2;

muscle_module3 = "red"
muscle_column3 = match(muscle_module3, muscle_modNames);
muscle_moduleGenes3 = muscle_moduleColors==muscle_module3;

muscle_module4 = "purple"
muscle_column4 = match(muscle_module4, muscle_modNames);
muscle_moduleGenes4 = muscle_moduleColors==muscle_module4;

muscle_module5 = "greenyellow"
muscle_column5 = match(muscle_module5, muscle_modNames);
muscle_moduleGenes5 = muscle_moduleColors==muscle_module5;

muscle_module6 = "blue"
muscle_column6 = match(muscle_module6, muscle_modNames);
muscle_moduleGenes6 = muscle_moduleColors==muscle_module6;


###Main figures - AGE
tiff('Fig-3B-brain-heart.tiff', units="cm", width=35, height=30, res=300, compression = 'lzw')
par(mfrow = c(2,2), mar=c(5,5,4,4));
plot(abs(brain_geneModuleMembership[brain_moduleGenes1, brain_column1]),
                   abs(brain_geneTraitSignificance_age[brain_moduleGenes1, 1]),
                   xlab = paste("Module membership", brain_module1, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Brain Tan  ", "bicor=0.82; ", italic("p-value"), "=6e-55")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = brain_module1, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(heart_geneModuleMembership[heart_moduleGenes1, heart_column1]),
                   abs(heart_geneTraitSignificance_age[heart_moduleGenes1, 1]),
                   xlab = paste("Module membership", heart_module1, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Heart Tan  ", "bicor=0.63; ", italic("p-value"), "=7e-15")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = heart_module1, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(heart_geneModuleMembership[heart_moduleGenes2, heart_column2]),
                   abs(heart_geneTraitSignificance_age[heart_moduleGenes2, 1]),
                   xlab = paste("Module membership", heart_module2, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Heart Blue  ", "bicor=0.51; ", italic("p-value"), "=1e-82")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = heart_module2, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
dev.off()
```

    ## png 
    ##   2

``` r
tiff('Fig-3B-liver-muscle.tiff', units="cm", width=35, height=30, res=300, compression = 'lzw')
par(mfrow = c(2,2), mar=c(5,5,4,4));
plot(abs(liver_geneModuleMembership[liver_moduleGenes1, liver_column1]),
                   abs(liver_geneTraitSignificance_age[liver_moduleGenes1, 1]),
                   xlab = paste("Module membership", liver_module1, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Liver Salmon  ", "bicor=0.65; ", italic("p-value"), "=5e-29")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module1, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes2, liver_column2]),
                   abs(liver_geneTraitSignificance_age[liver_moduleGenes2, 1]),
                   xlab = paste("Module membership", liver_module2, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Liver Darkturquoise  ", "bicor=0.64; ", italic("p-value"), "=5e-15")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8,col = liver_module2, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes1, muscle_column1]),
                   abs(muscle_geneTraitSignificance_age[muscle_moduleGenes1, 1]),
                   xlab = paste("Module membership", muscle_module1, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Muscle Magenta  ", "bicor=0.66; ", italic("p-value"), "=2e-26")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = muscle_module1, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes2, muscle_column2]),
                   abs(muscle_geneTraitSignificance_age[muscle_moduleGenes2, 1]),
                   xlab = paste("Module membership", muscle_module2, "module"),
                   ylab = "Gene significance for Age",
                   main = expression(paste("Muscle Brown  ", "bicor=0.71; ", italic("p-value"), "=2e-122")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = muscle_module2, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
dev.off()
```

    ## png 
    ##   2

``` r
###Supplemental figures - SEX
tiff('Fig-S3-brain-liver.tiff', units="cm", width=35, height=30, res=300, compression = 'lzw')
par(mfrow = c(2,2), mar=c(5,5,4,4));
plot(abs(brain_geneModuleMembership[brain_moduleGenes2, brain_column2]),
                   abs(brain_geneTraitSignificance_sex[brain_moduleGenes2, 1]),
                   xlab = paste("Module membership", brain_module2, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Brain Grey60  ", "bicor=0.89; ", italic("p-value"), "=8e-34")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = brain_module2, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes3, liver_column3]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes3, 1]),
                   xlab = paste("Module membership", liver_module3, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Tan  ", "bicor=0.54; ", italic("p-value"), "=2e-21")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8,col = liver_module3, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes4, liver_column4]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes4, 1]),
                   xlab = paste("Module membership", liver_module4, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Red  ", "bicor=0.73; ", italic("p-value"), "=5e-91")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module4, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes5, liver_column5]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes5, 1]),
                   xlab = paste("Module membership", liver_module5, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Darkolivegreen  ", "bicor=0.61; ", italic("p-value"), "=3e-08")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module5, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
dev.off()
```

    ## png 
    ##   2

``` r
tiff('Fig-S3-liver.tiff.tiff', units="cm", width=35, height=30, res=300, compression = 'lzw')
par(mfrow = c(2,2), mar=c(5,5,4,4));
plot(abs(liver_geneModuleMembership[liver_moduleGenes6, liver_column6]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes6, 1]),
                   xlab = paste("Module membership", liver_module6, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Darkgrey  ", "bicor=0.89; ", italic("p-value"), "=1e-40")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module6, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes7, liver_column7]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes7, 1]),
                   xlab = paste("Module membership", liver_module7, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Cyan  ", "bicor=0.56; ", italic("p-value"), "=7e-20")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8,col = liver_module7, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes8, liver_column8]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes8, 1]),
                   xlab = paste("Module membership", liver_module8, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Brown  ", "bicor=0.94; ", italic("p-value"), "=0")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module8, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(liver_geneModuleMembership[liver_moduleGenes9, liver_column9]),
                   abs(liver_geneTraitSignificance_sex[liver_moduleGenes9, 1]),
                   xlab = paste("Module membership", liver_module9, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Liver Blue  ", "bicor=0.87; ", italic("p-value"), "=4e-294")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = liver_module9, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
dev.off()
```

    ## png 
    ##   2

``` r
tiff('Fig-S3-muscle.tiff.tiff', units="cm", width=35, height=30, res=300, compression = 'lzw')
par(mfrow = c(2,2), mar=c(5,5,4,4));
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes3, muscle_column3]),
                   abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes3, 1]),
                   xlab = paste("Module membership", muscle_module3, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Muscle Red  ", "bicor=0.93; ", italic("p-value"), "=6e-210")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = muscle_module3, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes4, muscle_column4]),
                   abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes4, 1]),
                   xlab = paste("Module membership", muscle_module4, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Muscle Purple  ", "bicor=0.66; ", italic("p-value"), "=2e-23")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8,col = muscle_module4, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes5, muscle_column5]),
                   abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes5, 1]),
                   xlab = paste("Module membership", muscle_module5, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Muscle Greenyellow  ", "bicor=0.91; ", italic("p-value"), "=4e-53")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = muscle_module5, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
plot(abs(muscle_geneModuleMembership[muscle_moduleGenes6, muscle_column6]),
                   abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes6, 1]),
                   xlab = paste("Module membership", muscle_module6, "module"),
                   ylab = "Gene significance for Sex",
                   main = expression(paste("Muscle Blue  ", "bicor=0.56; ", italic("p-value"), "=5e-79")),
                   cex.main = 2.2, cex.lab = 1.8, cex.axis = 1.8, col = muscle_module6, pch=19, xlim=c(0,1),ylim=c(0,1))
abline(h=0.2, v=0.8, col="gray80", lty=1)
rect(0.8,0.2,2,2, col= rgb(0.8,0.8,0.8,alpha=0.5), border="gray80", lty=1)
dev.off()
```

    ## png 
    ##   2

## 3.4 Intramodular analysis: identifying genes with high GS and MM

Now we extract the list of each module's hub genes.

### Brain

``` r
###TAN module
# Retrieve each module's genes
brain_allTan <- as.data.frame(names(brainDataExpr)[brain_moduleColors=="tan"])
colnames(brain_allTan) <- "SYMBOL" 
rownames(brain_allTan) <- brain_allTan$SYMBOL 
nrow(brain_allTan) #number of genes included in the module
```

    ## [1] 220

``` r
write.table(brain_allTan, file="brain_allTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
brain_tan_filterGenes = abs(brain_geneTraitSignificance_age[brain_moduleGenes1, 1])>.2 & abs(brain_geneModuleMembership[brain_moduleGenes1, brain_column1])>.8
# Check how many genes are hub genes 
table(brain_tan_filterGenes) #TRUE indicates the number of hub genes
```

    ## brain_tan_filterGenes
    ## FALSE  TRUE 
    ##   186    34

``` r
# Check which genes are hub genes
brain_hubTan <- rownames(brain_allTan)[brain_tan_filterGenes]
brain_hubTan
```

    ##  [1] "B2m"       "C1qa"      "C1qb"      "C1qc"      "C3"        "C4b"      
    ##  [7] "Csf1"      "Ctsd"      "Ctsh"      "Ctss"      "Ctsz"      "Cx3cr1"   
    ## [13] "Gbp3"      "Gfap"      "H2-D1"     "H2-K1"     "H2-T23"    "Hexb"     
    ## [19] "Ifi27"     "Ifit3"     "Il33"      "Irf7"      "Itgb2"     "Lag3"     
    ## [25] "Laptm5"    "Lgals3"    "Lgals3bp"  "Lyz2"      "Neat1"     "Psmb8"    
    ## [31] "Serpina3n" "Slc11a1"   "Tap2"      "Tapbp"

``` r
all(brain_hubTan %in% brain_allTan$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(brain_hubTan, file="brain_hubTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###GREY60 module
# Retrieve each module's genes
brain_allGrey60 <- as.data.frame(names(brainDataExpr)[brain_moduleColors=="grey60"])
colnames(brain_allGrey60) <- "SYMBOL" 
rownames(brain_allGrey60) <- brain_allGrey60$SYMBOL 
nrow(brain_allGrey60) #number of genes included in the module
```

    ## [1] 97

``` r
write.table(brain_allGrey60, file="brain_allGrey60.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
brain_grey60_filterGenes = abs(brain_geneTraitSignificance_sex[brain_moduleGenes2, 1])>.2 & abs(brain_geneModuleMembership[brain_moduleGenes2, brain_column2])>.8
# Check how many genes are hub genes 
table(brain_grey60_filterGenes) #TRUE indicates the number of hub genes
```

    ## brain_grey60_filterGenes
    ## FALSE  TRUE 
    ##    60    37

``` r
# Check which genes are hub genes
brain_hubGrey60 <- rownames(brain_allGrey60)[brain_grey60_filterGenes]
brain_hubGrey60
```

    ##  [1] "Gm18796"     "Gm18797"     "Gm20788"     "Gm20830"     "Gm21064"    
    ##  [6] "Gm21292"     "Gm21719"     "Gm21721"     "Gm21854"     "Gm21865"    
    ## [11] "Gm21874"     "Gm28278"     "Gm28348"     "Gm28356"     "Gm28444"    
    ## [16] "Gm28445"     "Gm28507"     "Gm28510"     "Gm28587"     "Gm28597"    
    ## [21] "Gm28649"     "Gm28674"     "Gm28919"     "Gm29049"     "Gm29274"    
    ## [26] "Gm37222"     "Gm37236"     "Gm8446"      "Kdm5d"       "Tspy-ps"    
    ## [31] "Uba1y"       "Uba1y-ps2"   "Usp9y"       "Uty"         "Vmn2r-ps139"
    ## [36] "Zfy1"        "Zfy2"

``` r
all(brain_hubGrey60 %in% brain_allGrey60$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(brain_hubGrey60, file="brain_hubGrey60.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

#create table with module membership and gene significance values comprising only hub genes
brain_hubGrey60_geneSignificance <- brain_geneSignificance[rownames(brain_geneSignificance) %in% brain_hubGrey60,]
brain_hubGrey60_geneSignificance$Genes <- rownames(brain_hubGrey60_geneSignificance)
brain_hubGrey60_geneSignificance$Module <- "Grey60"

brain_hubTan_geneSignificance <- brain_geneSignificance[rownames(brain_geneSignificance) %in% brain_hubTan,]
brain_hubTan_geneSignificance$Genes <- rownames(brain_hubTan_geneSignificance)
brain_hubTan_geneSignificance$Module <- "Tan"

brain_hub_geneSignificance <- rbind(brain_hubTan_geneSignificance,brain_hubGrey60_geneSignificance)
write.table(brain_hub_geneSignificance, file="brain_hub_geneSignificance.txt", row.names = F, quote=FALSE, sep="\t")

brain_hubGrey60_moduleMembership <- brain_moduleMembership[rownames(brain_moduleMembership) %in% brain_hubGrey60, c("MM.grey60","MM.tan")]
brain_hubGrey60_moduleMembership$Genes <- rownames(brain_hubGrey60_moduleMembership)
brain_hubGrey60_moduleMembership$Module <- "Grey60"

brain_hubTan_moduleMembership <- brain_moduleMembership[rownames(brain_moduleMembership) %in% brain_hubTan, c("MM.grey60","MM.tan")]
brain_hubTan_moduleMembership$Genes <- rownames(brain_hubTan_moduleMembership)
brain_hubTan_moduleMembership$Module <- "Tan"

brain_hub_moduleMembership <- rbind(brain_hubTan_moduleMembership,brain_hubGrey60_moduleMembership)
write.table(brain_hub_moduleMembership, file="brain_hub_moduleMembership.txt", row.names = F, quote=FALSE, sep="\t")
```

### Heart

``` r
###TAN module
# Retrieve each module's genes
heart_allTan <- as.data.frame(names(heartDataExpr)[heart_moduleColors=="tan"])
colnames(heart_allTan) <- "SYMBOL" 
rownames(heart_allTan) <- heart_allTan$SYMBOL 
nrow(heart_allTan) #number of genes included in the module
```

    ## [1] 125

``` r
write.table(heart_allTan, file="heart_allTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
heart_tan_filterGenes = abs(heart_geneTraitSignificance_age[heart_moduleGenes1, 1])>.2 & abs(heart_geneModuleMembership[heart_moduleGenes1, heart_column1])>.8
# Check how many genes are hub genes 
table(heart_tan_filterGenes) #TRUE indicates the number of hub genes
```

    ## heart_tan_filterGenes
    ## FALSE  TRUE 
    ##   114    11

``` r
# Check which genes are hub genes
heart_hubTan <- rownames(heart_allTan)[heart_tan_filterGenes]
heart_hubTan
```

    ##  [1] "Acsm5"   "Amy1"    "Cd209f"  "Cds1"    "Ighg2c"  "Kcnk1"   "Pcdhb20"
    ##  [8] "Prkcq"   "Scn4b"   "Skap2"   "Vgll2"

``` r
all(heart_hubTan %in% heart_allTan$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(heart_hubTan, file="heart_hubTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###BLUE module
# Retrieve each module's genes
heart_allBlue <- as.data.frame(names(heartDataExpr)[heart_moduleColors=="blue"])
colnames(heart_allBlue) <- "SYMBOL" 
rownames(heart_allBlue) <- heart_allBlue$SYMBOL 
nrow(heart_allBlue) #number of genes included in the module
```

    ## [1] 1209

``` r
write.table(heart_allBlue, file="heart_allBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
heart_blue_filterGenes = abs(heart_geneTraitSignificance_age[heart_moduleGenes2, 1])>.2 & abs(heart_geneModuleMembership[heart_moduleGenes2, heart_column2])>.8
# Check how many genes are hub genes 
table(heart_blue_filterGenes) #TRUE indicates the number of hub genes
```

    ## heart_blue_filterGenes
    ## FALSE  TRUE 
    ##  1081   128

``` r
# Check which genes are hub genes
heart_hubBlue <- rownames(heart_allBlue)[heart_blue_filterGenes]
heart_hubBlue
```

    ##   [1] "Acaa1a"        "Acaa2"         "Actr1a"        "Adipor1"      
    ##   [5] "Adk"           "Adsl"          "Arf1"          "Atp5b"        
    ##   [9] "Atpaf1"        "Auh"           "Bcat2"         "Bsg"          
    ##  [13] "C030006K11Rik" "Capzb"         "Casq2"         "Cct5"         
    ##  [17] "Cd81"          "Cdc34"         "Cdc37"         "Chchd3"       
    ##  [21] "Ckm"           "Ckmt2"         "Cops7a"        "Cops8"        
    ##  [25] "Coq2"          "Coq5"          "Coq6"          "Coq9"         
    ##  [29] "Cox10"         "Ctbp1"         "Cyc1"          "Dele1"        
    ##  [33] "Des"           "Dnajb2"        "Dnpep"         "Echs1"        
    ##  [37] "Eci2"          "Ecsit"         "Egln2"         "Eif5a"        
    ##  [41] "Eno3"          "Etfb"          "Fars2"         "Fastk"        
    ##  [45] "Fbxw5"         "Fh1"           "Fkbp4"         "Gnpat"        
    ##  [49] "Gpi1"          "Hadh"          "Hmgcl"         "Idh3g"        
    ##  [53] "Immt"          "Isca1"         "Klhdc2"        "Ldha"         
    ##  [57] "Ldhb"          "Maf1"          "Map1lc3b"      "Mccc2"        
    ##  [61] "Mdh1"          "Mdh2"          "Mrpl37"        "Mrpl38"       
    ##  [65] "Mrpl4"         "Mrpl45"        "Mtfr1l"        "Myzap"        
    ##  [69] "Napa"          "Ndufa10"       "Ndufa9"        "Ndufs2"       
    ##  [73] "Ndufv1"        "Nfs1"          "Obscn"         "Pdhb"         
    ##  [77] "Pdrg1"         "Pgm2"          "Phb"           "Phyh"         
    ##  [81] "Pkm"           "Pmpcb"         "Poldip2"       "Popdc2"       
    ##  [85] "Ppp1ca"        "Ppp2r5d"       "Ppp5c"         "Prdx3"        
    ##  [89] "Prkaca"        "Psma1"         "Psmc4"         "Psmc5"        
    ##  [93] "Psmd2"         "Psmd3"         "Psmd7"         "Ptcd2"        
    ##  [97] "Ptges2"        "Pygm"          "Rnf187"        "Rpl3l"        
    ## [101] "Rrp1"          "Rxrg"          "Samm50"        "Sdhd"         
    ## [105] "Serf2"         "Sgca"          "Slc25a11"      "Slc25a12"     
    ## [109] "Slc25a39"      "Slc25a5"       "Slc2a4"        "Smpd1"        
    ## [113] "Snta1"         "Sod2"          "Stoml2"        "Suclg1"       
    ## [117] "Tcp1"          "Tmem70"        "Tprgl"         "Tufm"         
    ## [121] "Txn2"          "Ube2b"         "Ube2g2"        "Ubl7"         
    ## [125] "Uqcrc1"        "Uqcrc2"        "Wdr18"         "Yipf3"

``` r
all(heart_hubBlue %in% heart_allBlue$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(heart_hubBlue, file="heart_hubBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

#create table with module membership and gene significance values comprising only hub genes
heart_hubBlue_geneSignificance <- heart_geneSignificance[rownames(heart_geneSignificance) %in% heart_hubBlue,]
heart_hubBlue_geneSignificance$Genes <- rownames(heart_hubBlue_geneSignificance)
heart_hubBlue_geneSignificance$Module <- "Blue"

heart_hubTan_geneSignificance <- heart_geneSignificance[rownames(heart_geneSignificance) %in% heart_hubTan,]
heart_hubTan_geneSignificance$Genes <- rownames(heart_hubTan_geneSignificance)
heart_hubTan_geneSignificance$Module <- "Tan"

heart_hub_geneSignificance <- rbind(heart_hubTan_geneSignificance,heart_hubBlue_geneSignificance)
write.table(heart_hub_geneSignificance, file="heart_hub_geneSignificance.txt", row.names = F, quote=FALSE, sep="\t")

heart_hubBlue_moduleMembership <- heart_moduleMembership[rownames(heart_moduleMembership) %in% heart_hubBlue, c("MM.blue","MM.tan")]
heart_hubBlue_moduleMembership$Genes <- rownames(heart_hubBlue_moduleMembership)
heart_hubBlue_moduleMembership$Module <- "Blue"

heart_hubTan_moduleMembership <- heart_moduleMembership[rownames(heart_moduleMembership) %in% heart_hubTan, c("MM.blue","MM.tan")]
heart_hubTan_moduleMembership$Genes <- rownames(heart_hubTan_moduleMembership)
heart_hubTan_moduleMembership$Module <- "Tan"

heart_hub_moduleMembership <- rbind(heart_hubTan_moduleMembership,heart_hubBlue_moduleMembership)
write.table(heart_hub_moduleMembership, file="heart_hub_moduleMembership.txt", row.names = F, quote=FALSE, sep="\t")
```

### Muscle

``` r
###MAGENTA module
# Retrieve each module's genes
muscle_allMagenta <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="magenta"])
colnames(muscle_allMagenta) <- "SYMBOL" 
rownames(muscle_allMagenta) <- muscle_allMagenta$SYMBOL 
nrow(muscle_allMagenta) #number of genes included in the module
```

    ## [1] 203

``` r
write.table(muscle_allMagenta, file="muscle_allMagenta.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_magenta_filterGenes = abs(muscle_geneTraitSignificance_age[muscle_moduleGenes1, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes1, muscle_column1])>.8
# Check how many genes are hub genes 
table(muscle_magenta_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_magenta_filterGenes
    ## FALSE  TRUE 
    ##   194     9

``` r
# Check which genes are hub genes
muscle_hubMagenta <- rownames(muscle_allMagenta)[muscle_magenta_filterGenes]
muscle_hubMagenta
```

    ## [1] "Cpe"     "Eif3e"   "Itgb5"   "Kcmf1"   "Mib1"    "Plekhb1" "Rab2a"  
    ## [8] "Rasd2"   "Setd3"

``` r
all(muscle_hubMagenta %in% muscle_allMagenta$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubMagenta, file="muscle_hubMagenta.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###BROWN module
# Retrieve each module's genes
muscle_allBrown <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="brown"])
colnames(muscle_allBrown) <- "SYMBOL" 
rownames(muscle_allBrown) <- muscle_allBrown$SYMBOL 
nrow(muscle_allBrown) #number of genes included in the module
```

    ## [1] 806

``` r
write.table(muscle_allBrown, file="muscle_allBrown.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_brown_filterGenes = abs(muscle_geneTraitSignificance_age[muscle_moduleGenes2, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes2, muscle_column2])>.8
# Check how many genes are hub genes 
table(muscle_brown_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_brown_filterGenes
    ## FALSE  TRUE 
    ##   753    53

``` r
# Check which genes are hub genes
muscle_hubBrown <- rownames(muscle_allBrown)[muscle_brown_filterGenes]
muscle_hubBrown
```

    ##  [1] "Adamts2"  "Angptl1"  "Antxr2"   "Anxa2"    "Axl"      "C3"      
    ##  [7] "Cd34"     "Clec3b"   "Col1a1"   "Col1a2"   "Col3a1"   "Col5a1"  
    ## [13] "Col5a2"   "Col6a1"   "Col6a2"   "Col6a3"   "Ddr2"     "Dok2"    
    ## [19] "Dstn"     "Fbn1"     "Fn1"      "Fndc1"    "Fstl1"    "Igfbp6"  
    ## [25] "Islr"     "Itgbl1"   "Lpar1"    "Lrp1"     "Metrnl"   "Mfap5"   
    ## [31] "Mmp2"     "Mrc1"     "Mrc2"     "Ndn"      "Nid1"     "Olfml2b" 
    ## [37] "Olfml3"   "P3h3"     "Pcolce"   "Pcolce2"  "Pi16"     "Pltp"    
    ## [43] "Plxdc2"   "Rcn1"     "Rcn3"     "Rnase4"   "Serpinf1" "Serping1"
    ## [49] "Sod3"     "Sparc"    "Ssc5d"    "Tgfbi"    "Timp2"

``` r
all(muscle_hubBrown %in% muscle_allBrown$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubBrown, file="muscle_hubBrown.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###RED module
# Retrieve each module's genes
muscle_allRed <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="red"])
colnames(muscle_allRed) <- "SYMBOL" 
rownames(muscle_allRed) <- muscle_allRed$SYMBOL 
nrow(muscle_allRed) #number of genes included in the module
```

    ## [1] 484

``` r
write.table(muscle_allRed, file="muscle_allRed.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_red_filterGenes = abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes3, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes3, muscle_column3])>.8
# Check how many genes are hub genes 
table(muscle_red_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_red_filterGenes
    ## FALSE  TRUE 
    ##   444    40

``` r
# Check which genes are hub genes
muscle_hubRed <- rownames(muscle_allRed)[muscle_red_filterGenes]
muscle_hubRed
```

    ##  [1] "1700001O22Rik" "Amd1"          "Amd2"          "Anxa7"        
    ##  [5] "C7"            "Cbr2"          "Cdk19"         "Ddx3y"        
    ##  [9] "Eif2s3y"       "Eif4ebp1"      "Fam131a"       "Gm10032"      
    ## [13] "Gm12240"       "Gm8734"        "Grina"         "Hipk2"        
    ## [17] "Htra4"         "Irx3os"        "Kdm5d"         "Ldlr"         
    ## [21] "Mgst1"         "Musk"          "Nek6"          "Ppp1r14b"     
    ## [25] "Psmd8"         "Samd10"        "Sbk2"          "Serpinb6a"    
    ## [29] "Slc15a5"       "Slc2a3"        "Slc30a2"       "Spns2"        
    ## [33] "Stab2"         "Sub1"          "Tfcp2l1"       "Tmem37"       
    ## [37] "Ubl7"          "Uty"           "Vldlr"         "Wsb2"

``` r
all(muscle_hubRed %in% muscle_allRed$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubRed, file="muscle_hubRed.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###PURPLE module
# Retrieve each module's genes
muscle_allPurple <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="purple"])
colnames(muscle_allPurple) <- "SYMBOL" 
rownames(muscle_allPurple) <- muscle_allPurple$SYMBOL 
nrow(muscle_allPurple) #number of genes included in the module
```

    ## [1] 174

``` r
write.table(muscle_allPurple, file="muscle_allPurple.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_purple_filterGenes = abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes4, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes4, muscle_column4])>.8
# Check how many genes are hub genes 
table(muscle_purple_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_purple_filterGenes
    ## FALSE  TRUE 
    ##   164    10

``` r
# Check which genes are hub genes
muscle_hubPurple <- rownames(muscle_allPurple)[muscle_purple_filterGenes]
muscle_hubPurple
```

    ##  [1] "Acox1"    "Aldh2"    "Ces1d"    "Egf"      "Gcdh"     "Gstm2"   
    ##  [7] "Mybpc1"   "Ppp1r1a"  "Selenbp1" "Suclg2"

``` r
all(muscle_hubPurple %in% muscle_allPurple$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubPurple, file="muscle_hubPurple.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###GREENYELLOW module
# Retrieve each module's genes
muscle_allGreenyellow <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="greenyellow"])
colnames(muscle_allGreenyellow) <- "SYMBOL" 
rownames(muscle_allGreenyellow) <- muscle_allGreenyellow$SYMBOL 
nrow(muscle_allGreenyellow) #number of genes included in the module
```

    ## [1] 140

``` r
write.table(muscle_allGreenyellow, file="muscle_allGreenyellow.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_greenyellow_filterGenes = abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes5, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes5, muscle_column5])>.8
# Check how many genes are hub genes 
table(muscle_greenyellow_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_greenyellow_filterGenes
    ## FALSE  TRUE 
    ##   129    11

``` r
# Check which genes are hub genes
muscle_hubGreenyellow <- rownames(muscle_allGreenyellow)[muscle_greenyellow_filterGenes]
muscle_hubGreenyellow
```

    ##  [1] "Cyfip2"  "Ddx3x"   "Eif2s3x" "Gm15337" "Gm47708" "Homer2"  "Lamb2"  
    ##  [8] "Mybph"   "Neu2"    "Padi2"   "Xist"

``` r
all(muscle_hubGreenyellow %in% muscle_allGreenyellow$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubGreenyellow, file="muscle_hubGreenyellow.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###BLUE module
# Retrieve each module's genes
muscle_allBlue <- as.data.frame(names(muscleDataExpr)[muscle_moduleColors=="blue"])
colnames(muscle_allBlue) <- "SYMBOL" 
rownames(muscle_allBlue) <- muscle_allBlue$SYMBOL 
nrow(muscle_allBlue) #number of genes included in the module
```

    ## [1] 949

``` r
write.table(muscle_allBlue, file="muscle_allBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
muscle_blue_filterGenes = abs(muscle_geneTraitSignificance_sex[muscle_moduleGenes6, 1])>.2 & abs(muscle_geneModuleMembership[muscle_moduleGenes6, muscle_column6])>.8
# Check how many genes are hub genes 
table(muscle_blue_filterGenes) #TRUE indicates the number of hub genes
```

    ## muscle_blue_filterGenes
    ## FALSE  TRUE 
    ##   796   153

``` r
# Check which genes are hub genes
muscle_hubBlue <- rownames(muscle_allBlue)[muscle_blue_filterGenes]
muscle_hubBlue
```

    ##   [1] "0610012G03Rik" "Acot13"        "Atp5c1"        "Atp5d"        
    ##   [5] "Atp5e"         "Atp5g1"        "Atp5g2"        "Atp5g3"       
    ##   [9] "Atp5h"         "Atp5j"         "Atp5j2"        "Atp5l"        
    ##  [13] "Atp5mpl"       "Atp5o"         "Chchd1"        "Chchd10"      
    ##  [17] "Cisd1"         "Cox14"         "Cox4i1"        "Cox5a"        
    ##  [21] "Cox5b"         "Cox6a2"        "Cox6b1"        "Cox6c"        
    ##  [25] "Cox7a1"        "Cox7a2"        "Cox7c"         "Cox8a"        
    ##  [29] "Cox8b"         "D8Ertd738e"    "Elob"          "Fau"          
    ##  [33] "Fmc1"          "Fxyd1"         "Gng5"          "Higd1a"       
    ##  [37] "Hint2"         "Krtcap2"       "Lamtor2"       "Lars2"        
    ##  [41] "Mpc2"          "Mrpl14"        "Mrpl23"        "Mrpl33"       
    ##  [45] "Mrpl51"        "Mrps18a"       "Mrps21"        "Mrps24"       
    ##  [49] "Mylpf"         "Ndufa11"       "Ndufa12"       "Ndufa13"      
    ##  [53] "Ndufa2"        "Ndufa4"        "Ndufa5"        "Ndufa6"       
    ##  [57] "Ndufa7"        "Ndufa8"        "Ndufb10"       "Ndufb11"      
    ##  [61] "Ndufb2"        "Ndufb4"        "Ndufb5"        "Ndufb6"       
    ##  [65] "Ndufb7"        "Ndufb8"        "Ndufb9"        "Ndufc1"       
    ##  [69] "Ndufc2"        "Ndufs5"        "Ndufs6"        "Ndufs7"       
    ##  [73] "Ndufv3"        "Nedd8"         "Nenf"          "Nop10"        
    ##  [77] "Pam16"         "Pdcd5"         "Pfdn5"         "Polr2l"       
    ##  [81] "Psmb1"         "Psmb3"         "Psmb4"         "Psmb5"        
    ##  [85] "Pvalb"         "Romo1"         "Rpl10a"        "Rpl11"        
    ##  [89] "Rpl12"         "Rpl13"         "Rpl14"         "Rpl18"        
    ##  [93] "Rpl18a"        "Rpl19"         "Rpl21"         "Rpl23"        
    ##  [97] "Rpl26"         "Rpl27"         "Rpl28"         "Rpl30"        
    ## [101] "Rpl31"         "Rpl32"         "Rpl34"         "Rpl35a"       
    ## [105] "Rpl36"         "Rpl36a"        "Rpl37"         "Rpl37a"       
    ## [109] "Rpl38"         "Rpl41"         "Rpl7a"         "Rpl8"         
    ## [113] "Rplp1"         "Rplp2"         "Rps11"         "Rps13"        
    ## [117] "Rps14"         "Rps15"         "Rps16"         "Rps17"        
    ## [121] "Rps18"         "Rps19"         "Rps20"         "Rps21"        
    ## [125] "Rps23"         "Rps24"         "Rps25"         "Rps26"        
    ## [129] "Rps27"         "Rps27a"        "Rps27l"        "Rps28"        
    ## [133] "Rps3"          "Rps5"          "Rps8"          "Rps9"         
    ## [137] "S100a1"        "Sem1"          "Smim26"        "Timm13"       
    ## [141] "Timm8b"        "Tmem147"       "Tmem256"       "Tnni2"        
    ## [145] "Tomm7"         "Tpt1"          "Uqcc2"         "Uqcr10"       
    ## [149] "Uqcr11"        "Uqcrb"         "Uqcrh"         "Uqcrq"        
    ## [153] "mt-Nd2"

``` r
all(muscle_hubBlue %in% muscle_allBlue$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(muscle_hubBlue, file="muscle_hubBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

#create table with module membership and gene significance values comprising only hub genes
muscle_hubMagenta_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubMagenta,]
muscle_hubMagenta_geneSignificance$Genes <- rownames(muscle_hubMagenta_geneSignificance)
muscle_hubMagenta_geneSignificance$Module <- "Magenta"

muscle_hubBrown_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubBrown,]
muscle_hubBrown_geneSignificance$Genes <- rownames(muscle_hubBrown_geneSignificance)
muscle_hubBrown_geneSignificance$Module <- "Brown"

muscle_hubRed_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubRed,]
muscle_hubRed_geneSignificance$Genes <- rownames(muscle_hubRed_geneSignificance)
muscle_hubRed_geneSignificance$Module <- "Red"

muscle_hubPurple_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubPurple,]
muscle_hubPurple_geneSignificance$Genes <- rownames(muscle_hubPurple_geneSignificance)
muscle_hubPurple_geneSignificance$Module <- "Purple"

muscle_hubGreenyellow_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubGreenyellow,]
muscle_hubGreenyellow_geneSignificance$Genes <- rownames(muscle_hubGreenyellow_geneSignificance)
muscle_hubGreenyellow_geneSignificance$Module <- "Greenyellow"

muscle_hubBlue_geneSignificance <- muscle_geneSignificance[rownames(muscle_geneSignificance) %in% muscle_hubBlue,]
muscle_hubBlue_geneSignificance$Genes <- rownames(muscle_hubBlue_geneSignificance)
muscle_hubBlue_geneSignificance$Module <- "Blue"

muscle_hub_geneSignificance <- rbind(muscle_hubMagenta_geneSignificance,muscle_hubBrown_geneSignificance,muscle_hubRed_geneSignificance,muscle_hubPurple_geneSignificance,muscle_hubGreenyellow_geneSignificance,muscle_hubBlue_geneSignificance)
write.table(muscle_hub_geneSignificance, file="muscle_hub_geneSignificance.txt", row.names = F, quote=FALSE, sep="\t")

muscle_hubMagenta_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubMagenta, c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubMagenta_moduleMembership$Genes <- rownames(muscle_hubMagenta_moduleMembership)
muscle_hubMagenta_moduleMembership$Module <- "Magenta"

muscle_hubBrown_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubBrown,c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubBrown_moduleMembership$Genes <- rownames(muscle_hubBrown_moduleMembership)
muscle_hubBrown_moduleMembership$Module <- "Brown"

muscle_hubRed_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubRed,c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubRed_moduleMembership$Genes <- rownames(muscle_hubRed_moduleMembership)
muscle_hubRed_moduleMembership$Module <- "Red"

muscle_hubPurple_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubPurple,c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubPurple_moduleMembership$Genes <- rownames(muscle_hubPurple_moduleMembership)
muscle_hubPurple_moduleMembership$Module <- "Purple"

muscle_hubGreenyellow_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubGreenyellow,c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubGreenyellow_moduleMembership$Genes <- rownames(muscle_hubGreenyellow_moduleMembership)
muscle_hubGreenyellow_moduleMembership$Module <- "Greenyellow"

muscle_hubBlue_moduleMembership <- muscle_moduleMembership[rownames(muscle_moduleMembership) %in% muscle_hubBlue,c("MM.magenta","MM.brown","MM.red","MM.purple","MM.greenyellow","MM.blue")]
muscle_hubBlue_moduleMembership$Genes <- rownames(muscle_hubBlue_moduleMembership)
muscle_hubBlue_moduleMembership$Module <- "Blue"

muscle_hub_moduleMembership <- rbind(muscle_hubMagenta_moduleMembership,muscle_hubBrown_moduleMembership,muscle_hubRed_moduleMembership,muscle_hubPurple_moduleMembership,muscle_hubGreenyellow_moduleMembership,muscle_hubBlue_moduleMembership)
write.table(muscle_hub_moduleMembership, file="muscle_hub_moduleMembership.txt", row.names = F, quote=FALSE, sep="\t")
```

### Liver

``` r
###SALMON module
# Retrieve each module's genes
liver_allSalmon <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="salmon"])
colnames(liver_allSalmon) <- "SYMBOL" 
rownames(liver_allSalmon) <- liver_allSalmon$SYMBOL 
nrow(liver_allSalmon) #number of genes included in the module
```

    ## [1] 232

``` r
write.table(liver_allSalmon, file="liver_allSalmon.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_salmon_filterGenes = abs(liver_geneTraitSignificance_age[liver_moduleGenes1, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes1, liver_column1])>.8
# Check how many genes are hub genes 
table(liver_salmon_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_salmon_filterGenes
    ## FALSE  TRUE 
    ##   227     5

``` r
# Check which genes are hub genes
liver_hubSalmon <- rownames(liver_allSalmon)[liver_salmon_filterGenes]
liver_hubSalmon
```

    ## [1] "Ccl5"   "H2-Aa"  "H2-Eb1" "Ntrk2"  "Slamf7"

``` r
all(liver_hubSalmon %in% liver_allSalmon$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubSalmon, file="liver_hubSalmon.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###DARKTURQUOISE module
# Retrieve each module's genes
liver_allDarkturquoise <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="darkturquoise"])
colnames(liver_allDarkturquoise) <- "SYMBOL" 
rownames(liver_allDarkturquoise) <- liver_allDarkturquoise$SYMBOL 
nrow(liver_allDarkturquoise) #number of genes included in the module
```

    ## [1] 118

``` r
write.table(liver_allDarkturquoise, file="liver_allDarkturquoise.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_darkturquoise_filterGenes = abs(liver_geneTraitSignificance_age[liver_moduleGenes2, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes2, liver_column2])>.8
# Check how many genes are hub genes 
table(liver_darkturquoise_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_darkturquoise_filterGenes
    ## FALSE  TRUE 
    ##   103    15

``` r
# Check which genes are hub genes
liver_hubDarkturquoise <- rownames(liver_allDarkturquoise)[liver_darkturquoise_filterGenes]
liver_hubDarkturquoise
```

    ##  [1] "Cd19"     "Cd79a"    "Cd79b"    "Ighg2b"   "Ighg2c"   "Ighm"    
    ##  [7] "Ighv1-53" "Igkc"     "Igkv3-2"  "Igkv3-5"  "Iglc1"    "Iglc2"   
    ## [13] "Iglv1"    "Jchain"   "Mzb1"

``` r
all(liver_hubDarkturquoise %in% liver_allDarkturquoise$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubDarkturquoise, file="liver_hubDarkturquoise.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###TAN module
# Retrieve each module's genes
liver_allTan <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="tan"])
colnames(liver_allTan) <- "SYMBOL" 
rownames(liver_allTan) <- liver_allTan$SYMBOL 
nrow(liver_allTan) #number of genes included in the module
```

    ## [1] 267

``` r
write.table(liver_allTan, file="liver_allTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_tan_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes3, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes3, liver_column3])>.8
# Check how many genes are hub genes 
table(liver_tan_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_tan_filterGenes
    ## FALSE  TRUE 
    ##   253    14

``` r
# Check which genes are hub genes
liver_hubTan <- rownames(liver_allTan)[liver_tan_filterGenes]
liver_hubTan
```

    ##  [1] "Babam1"  "Cope"    "Gps1"    "Grhpr"   "Mbl1"    "Mgst1"   "Mup1"   
    ##  [8] "Nme1"    "Psmc4"   "Psmd13"  "Psmd4"   "Psmd6"   "Rarres1" "Rexo2"

``` r
all(liver_hubTan %in% liver_allTan$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubTan, file="liver_hubTan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###RED module
# Retrieve each module's genes
liver_allRed <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="red"])
colnames(liver_allRed) <- "SYMBOL" 
rownames(liver_allRed) <- liver_allRed$SYMBOL 
nrow(liver_allRed) #number of genes included in the module
```

    ## [1] 536

``` r
write.table(liver_allRed, file="liver_allRed.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_red_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes4, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes4, liver_column4])>.8
# Check how many genes are hub genes 
table(liver_red_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_red_filterGenes
    ## FALSE  TRUE 
    ##   503    33

``` r
# Check which genes are hub genes
liver_hubRed <- rownames(liver_allRed)[liver_red_filterGenes]
liver_hubRed
```

    ##  [1] "Abcd1"         "Acad9"         "Acss3"         "Agmo"         
    ##  [5] "Akr1c20"       "Akr1d1"        "Aldh9a1"       "Apol7a"       
    ##  [9] "B630019A10Rik" "Bphl"          "Car5a"         "Ces1g"        
    ## [13] "Dhrs7"         "Dpys"          "Echs1"         "Erg28"        
    ## [17] "Gm4756"        "Gstt3"         "Hadh"          "Macrod1"      
    ## [21] "Mpc1"          "Mpc1-ps"       "Mup-ps16"      "Ndrg2"        
    ## [25] "Pecr"          "Plpp3"         "Plscr2"        "Pon1"         
    ## [29] "Serpinb1a"     "Shmt1"         "Slc47a1"       "Tox"          
    ## [33] "Vnn3"

``` r
all(liver_hubRed %in% liver_allRed$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubRed, file="liver_hubRed.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###DARKOLIVEGREEN module
# Retrieve each module's genes
liver_allDarkolivegreen <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="darkolivegreen"])
colnames(liver_allDarkolivegreen) <- "SYMBOL" 
rownames(liver_allDarkolivegreen) <- liver_allDarkolivegreen$SYMBOL 
nrow(liver_allDarkolivegreen) #number of genes included in the module
```

    ## [1] 70

``` r
write.table(liver_allDarkolivegreen, file="liver_allDarkolivegreen.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_darkolivegreen_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes5, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes5, liver_column5])>.8
# Check how many genes are hub genes 
table(liver_darkolivegreen_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_darkolivegreen_filterGenes
    ## FALSE  TRUE 
    ##    59    11

``` r
# Check which genes are hub genes
liver_hubDarkolivegreen <- rownames(liver_allDarkolivegreen)[liver_darkolivegreen_filterGenes]
liver_hubDarkolivegreen
```

    ##  [1] "Acss2"  "Dhcr7"  "Fdft1"  "Hmgcs1" "Mmab"   "Msmo1"  "Mvd"    "Mvk"   
    ##  [9] "Pmvk"   "Rdh11"  "Spns2"

``` r
all(liver_hubDarkolivegreen %in% liver_allDarkolivegreen$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubDarkolivegreen, file="liver_hubDarkolivegreen.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###DARKGREY module
# Retrieve each module's genes
liver_allDarkgrey <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="darkgrey"])
colnames(liver_allDarkgrey) <- "SYMBOL" 
rownames(liver_allDarkgrey) <- liver_allDarkgrey$SYMBOL 
nrow(liver_allDarkgrey) #number of genes included in the module
```

    ## [1] 116

``` r
write.table(liver_allDarkgrey, file="liver_allDarkgrey.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_darkgrey_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes6, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes6, liver_column6])>.8
# Check how many genes are hub genes 
table(liver_darkgrey_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_darkgrey_filterGenes
    ## FALSE  TRUE 
    ##    92    24

``` r
# Check which genes are hub genes
liver_hubDarkgrey <- rownames(liver_allDarkgrey)[liver_darkgrey_filterGenes]
liver_hubDarkgrey
```

    ##  [1] "9130409I23Rik" "Arsa"          "Chpt1"         "Cidec"        
    ##  [5] "Clstn3"        "Cox19"         "Cyp2u1"        "Dpy19l3"      
    ##  [9] "Fancl"         "Fitm1"         "Gpc1"          "Gprc5b"       
    ## [13] "Nat8"          "Ntrk1"         "Olig1"         "Osbpl3"       
    ## [17] "Pard3b"        "Rassf3"        "Snhg11"        "Unc119"       
    ## [21] "Uox"           "Zfp979"        "Zfp982"        "Zfp992"

``` r
all(liver_hubDarkgrey %in% liver_allDarkgrey$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubDarkgrey, file="liver_hubDarkgrey.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###CYAN module
# Retrieve each module's genes
liver_allCyan <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="cyan"])
colnames(liver_allCyan) <- "SYMBOL" 
rownames(liver_allCyan) <- liver_allCyan$SYMBOL 
nrow(liver_allCyan) #number of genes included in the module
```

    ## [1] 223

``` r
write.table(liver_allCyan, file="liver_allCyan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_cyan_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes7, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes7, liver_column7])>.8
# Check how many genes are hub genes 
table(liver_cyan_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_cyan_filterGenes
    ## FALSE  TRUE 
    ##   210    13

``` r
# Check which genes are hub genes
liver_hubCyan <- rownames(liver_allCyan)[liver_cyan_filterGenes]
liver_hubCyan
```

    ##  [1] "Arcn1"   "Copg1"   "Creld2"  "Hspa5"   "Iars"    "Manf"    "Sdf2l1" 
    ##  [8] "Sec22b"  "Sec24d"  "Sec61a1" "Serp1"   "Slc33a1" "Ssr1"

``` r
all(liver_hubCyan %in% liver_allCyan$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubCyan, file="liver_hubCyan.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###BROWN module
# Retrieve each module's genes
liver_allBrown <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="brown"])
colnames(liver_allBrown) <- "SYMBOL" 
rownames(liver_allBrown) <- liver_allBrown$SYMBOL 
nrow(liver_allBrown) #number of genes included in the module
```

    ## [1] 791

``` r
write.table(liver_allBrown, file="liver_allBrown.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_brown_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes8, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes8, liver_column8])>.8
# Check how many genes are hub genes 
table(liver_brown_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_brown_filterGenes
    ## FALSE  TRUE 
    ##   723    68

``` r
# Check which genes are hub genes
liver_hubBrown <- rownames(liver_allBrown)[liver_brown_filterGenes]
liver_hubBrown
```

    ##  [1] "A1bg"        "Acot3"       "Aldh3b3"     "Arrdc4"      "Atp6v0d2"   
    ##  [6] "Chic1"       "Cux2"        "Cyp17a1"     "Cyp2a22"     "Cyp2a4"     
    ## [11] "Cyp2b10"     "Cyp2b9"      "Cyp2c38"     "Cyp2c39"     "Cyp2c68"    
    ## [16] "Cyp2c69"     "Cyp2g1"      "Cyp3a41a"    "Dqx1"        "Echdc3"     
    ## [21] "Eci3"        "Esr1"        "Fmo1"        "Fmo2"        "Fmo3"       
    ## [26] "Fmo4"        "Gm11695"     "Gm37273"     "Gm42375"     "Gm6135"     
    ## [31] "Gypc"        "Hamp2"       "Hao2"        "Hexb"        "Hpd"        
    ## [36] "Ildr2"       "Kat6b-ps2"   "Klhl13"      "Maob"        "Nipal1"     
    ## [41] "Nt5e"        "Papss2"      "Prlr"        "Rdh16f2"     "Rtn4"       
    ## [46] "Sall1"       "Sh2d4a"      "Slc16a5"     "Slc22a26"    "Slc22a27"   
    ## [51] "Slco1a4"     "St3gal6"     "Sult1a1"     "Sult1d1"     "Sult2a-ps1" 
    ## [56] "Sult2a-ps2"  "Sult2a1"     "Sult2a2"     "Sult2a7"     "Sult3a1"    
    ## [61] "Tcn2"        "Tm6sf2"      "Tmem167-ps2" "Tmem98"      "Uba7"       
    ## [66] "Vldlr"       "Xist"        "Zbed4-ps2"

``` r
all(liver_hubBrown %in% liver_allBrown$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubBrown, file="liver_hubBrown.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

###BLUE module
# Retrieve each module's genes
liver_allBlue <- as.data.frame(names(liverDataExpr)[liver_moduleColors=="blue"])
colnames(liver_allBlue) <- "SYMBOL" 
rownames(liver_allBlue) <- liver_allBlue$SYMBOL 
nrow(liver_allBlue) #number of genes included in the module
```

    ## [1] 954

``` r
write.table(liver_allBlue, file="liver_allBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

# Retrieve each module's hub genes (GS > 0.2 and MM >0.8)
liver_blue_filterGenes = abs(liver_geneTraitSignificance_sex[liver_moduleGenes6, 1])>.2 & abs(liver_geneModuleMembership[liver_moduleGenes6, liver_column6])>.8
# Check how many genes are hub genes 
table(liver_blue_filterGenes) #TRUE indicates the number of hub genes
```

    ## liver_blue_filterGenes
    ## FALSE  TRUE 
    ##    92    24

``` r
# Check which genes are hub genes
liver_hubBlue <- rownames(liver_allBlue)[liver_blue_filterGenes]
liver_hubBlue
```

    ##   [1] "1110038B12Rik" "1810013L24Rik" "5730455P16Rik" "9130401M01Rik"
    ##   [5] "AI182371"      "AI463229"      "Abcb10"        "Abcg2"        
    ##   [9] "Acox1"         "Actr1b"        "Ankrd27"       "Ankrd52"      
    ##  [13] "Arfgap3"       "Arl2bp"        "Arsg"          "Asap3"        
    ##  [17] "Asb6"          "Bcap31"        "Bpnt1"         "C6"           
    ##  [21] "C730027H18Rik" "Cadm4"         "Caml"          "Capn7"        
    ##  [25] "Capza2"        "Cbfa2t2"       "Cdc34"         "Cdh15"        
    ##  [29] "Cdhr5"         "Cdipt"         "Ces2a"         "Ces2b"        
    ##  [33] "Chmp1b"        "Chrd"          "Cmtm6"         "Cnot7"        
    ##  [37] "Csad"          "Csnk1d"        "Ctr9"          "Cul2"         
    ##  [41] "Cxadr"         "Cyp4a12a"      "Dcakd"         "Dnase2b"      
    ##  [45] "Dnd1"          "Dus1l"         "Dusp8"         "Dync1h1"      
    ##  [49] "Dynll2"        "Eed"           "Eif5"          "Elovl2"       
    ##  [53] "Elovl3"        "Emc6"          "Ephx1"         "Ephx2"        
    ##  [57] "Eri1"          "Etfbkmt"       "Fcor"          "Fech"         
    ##  [61] "Ftl2-ps"       "Fyttd1"        "Galnt1"        "Garem1"       
    ##  [65] "Gas2l1"        "Gm11963"       "Gm15883"       "Gm31036"      
    ##  [69] "Gm34654"       "Gm40787"       "Gm42688"       "Gm44507"      
    ##  [73] "Gm45727"       "Gm49338"       "Gpat4"         "Gpr39"        
    ##  [77] "Gprin3"        "Gpsm2"         "Gspt1"         "Gsr"          
    ##  [81] "Gt(ROSA)26Sor" "Gtpbp4-ps1"    "Hras"          "Hsd17b12"     
    ##  [85] "Igsf5"         "Ikbke"         "Inhbc"         "Inpp4a"       
    ##  [89] "Ints9"         "Kdm4b"         "Klk1b4"        "Lcp1"         
    ##  [93] "Ldah"          "Lmo4"          "Lrrc3"         "Lrrc42"       
    ##  [97] "Lurap1l"       "Map1lc3b"      "Mbnl2"         "Mcm2"         
    ## [101] "Mcph1"         "Mctp2"         "Mecr"          "Med26"        
    ## [105] "Mfsd11"        "Mgrn1"         "Msrb3"         "Mta2"         
    ## [109] "Mup6"          "Mup9"          "Myo1c"         "Myo6"         
    ## [113] "Mzt1"          "Nfyc"          "Npepps"        "Ogdh"         
    ## [117] "Ola1"          "Pafah1b3"      "Paip2"         "Pak1ip1"      
    ## [121] "Paqr7"         "Pdilt"         "Pgs1"          "Phf20l1"      
    ## [125] "Phlda2"        "Phyh"          "Pip5k1a"       "Pisd"         
    ## [129] "Plekha1"       "Plk2"          "Ppp2ca"        "Ppp2r5c"      
    ## [133] "Psmd10"        "Psmd14"        "Psme3"         "Ptma"         
    ## [137] "Pum2"          "Rap2a"         "Rest"          "Rpf1"         
    ## [141] "Rpl26-ps6"     "Sae1"          "Sall2"         "Samd1"        
    ## [145] "Sars"          "Scp2"          "Selenbp2"      "Sephs1"       
    ## [149] "Serinc1"       "Serpina11"     "Serpina4-ps1"  "Serpine2"     
    ## [153] "Sfr1"          "Sgpp1"         "Slc35b3"       "Slc35e3"      
    ## [157] "Smpd1"         "Snhg3"         "Socs4"         "Sort1"        
    ## [161] "Spaca6"        "Stat6"         "Susd4"         "Tesk2"        
    ## [165] "Tex30"         "Thoc6"         "Thtpa"         "Thumpd3"      
    ## [169] "Tm9sf1"        "Tmem125"       "Tmigd1"        "Topors"       
    ## [173] "Tpmt"          "Tpp1"          "Traf2"         "Tram1"        
    ## [177] "Tsg101"        "Ttc33"         "Uap1l1"        "Ubac2"        
    ## [181] "Ufm1"          "Ugdh"          "Ugt2b38"       "Ugt2b5"       
    ## [185] "Uri1"          "Vkorc1l1"      "Vrk3"          "Yipf3"        
    ## [189] "Ythdf1"        "Zbtb42"        "Zc3h12a"       "Zc3h14"       
    ## [193] "Zdhhc16"       "Zfp125"        "Zfp706"        "Zfpm1"        
    ## [197] "Zfyve1"        "Zkscan3"       "Zranb2"        "Zyg11a"

``` r
all(liver_hubBlue %in% liver_allBlue$SYMBOL)
```

    ## [1] TRUE

``` r
write.table(liver_hubBlue, file="liver_hubBlue.txt", row.names = F, col.names = "Gene", quote=FALSE, sep="\t")

#create table with module membership and gene significance values comprising only hub genes
liver_hubTan_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubTan,]
liver_hubTan_geneSignificance$Genes <- rownames(liver_hubTan_geneSignificance)
liver_hubTan_geneSignificance$Module <- "Tan"

liver_hubSalmon_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubSalmon,]
liver_hubSalmon_geneSignificance$Genes <- rownames(liver_hubSalmon_geneSignificance)
liver_hubSalmon_geneSignificance$Module <- "Salmon"

liver_hubRed_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubRed,]
liver_hubRed_geneSignificance$Genes <- rownames(liver_hubRed_geneSignificance)
liver_hubRed_geneSignificance$Module <- "Red"

liver_hubDarkturquoise_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubDarkturquoise,]
liver_hubDarkturquoise_geneSignificance$Genes <- rownames(liver_hubDarkturquoise_geneSignificance)
liver_hubDarkturquoise_geneSignificance$Module <- "Darkturquoise"

liver_hubDarkolivegreen_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubDarkolivegreen,]
liver_hubDarkolivegreen_geneSignificance$Genes <- rownames(liver_hubDarkolivegreen_geneSignificance)
liver_hubDarkolivegreen_geneSignificance$Module <- "Darkolivegreen"

liver_hubDarkgrey_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubDarkgrey,]
liver_hubDarkgrey_geneSignificance$Genes <- rownames(liver_hubDarkgrey_geneSignificance)
liver_hubDarkgrey_geneSignificance$Module <- "Darkgrey"

liver_hubCyan_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubCyan,]
liver_hubCyan_geneSignificance$Genes <- rownames(liver_hubCyan_geneSignificance)
liver_hubCyan_geneSignificance$Module <- "Cyan"

liver_hubBrown_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubBrown,]
liver_hubBrown_geneSignificance$Genes <- rownames(liver_hubBrown_geneSignificance)
liver_hubBrown_geneSignificance$Module <- "Brown"

liver_hubBlue_geneSignificance <- liver_geneSignificance[rownames(liver_geneSignificance) %in% liver_hubBlue,]
liver_hubBlue_geneSignificance$Genes <- rownames(liver_hubBlue_geneSignificance)
liver_hubBlue_geneSignificance$Module <- "Blue"

liver_hub_geneSignificance <- rbind(liver_hubTan_geneSignificance,liver_hubSalmon_geneSignificance,liver_hubRed_geneSignificance,liver_hubDarkturquoise_geneSignificance,liver_hubDarkolivegreen_geneSignificance,liver_hubDarkgrey_geneSignificance,liver_hubCyan_geneSignificance,liver_hubBrown_geneSignificance,liver_hubBlue_geneSignificance)
write.table(liver_hub_geneSignificance, file="liver_hub_geneSignificance.txt", row.names = F, quote=FALSE, sep="\t")

liver_hubTan_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubTan,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubTan_moduleMembership$Genes <- rownames(liver_hubTan_moduleMembership)
liver_hubTan_moduleMembership$Module <- "Tan"

liver_hubSalmon_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubSalmon,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubSalmon_moduleMembership$Genes <- rownames(liver_hubSalmon_moduleMembership)
liver_hubSalmon_moduleMembership$Module <- "Salmon"

liver_hubRed_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubRed,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubRed_moduleMembership$Genes <- rownames(liver_hubRed_moduleMembership)
liver_hubRed_moduleMembership$Module <- "Red"

liver_hubDarkturquoise_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubDarkturquoise,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubDarkturquoise_moduleMembership$Genes <- rownames(liver_hubDarkturquoise_moduleMembership)
liver_hubDarkturquoise_moduleMembership$Module <- "Darkturquoise"

liver_hubDarkolivegreen_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubDarkolivegreen,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubDarkolivegreen_moduleMembership$Genes <- rownames(liver_hubDarkolivegreen_moduleMembership)
liver_hubDarkolivegreen_moduleMembership$Module <- "Darkolivegreen"

liver_hubDarkgrey_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubDarkgrey,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubDarkgrey_moduleMembership$Genes <- rownames(liver_hubDarkgrey_moduleMembership)
liver_hubDarkgrey_moduleMembership$Module <- "Darkgrey"

liver_hubCyan_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubCyan,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubCyan_moduleMembership$Genes <- rownames(liver_hubCyan_moduleMembership)
liver_hubCyan_moduleMembership$Module <- "Cyan"

liver_hubBrown_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubBrown,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubBrown_moduleMembership$Genes <- rownames(liver_hubBrown_moduleMembership)
liver_hubBrown_moduleMembership$Module <- "Brown"

liver_hubBlue_moduleMembership <- liver_moduleMembership[rownames(liver_moduleMembership) %in% liver_hubBlue,c("MM.tan","MM.salmon","MM.red","MM.darkturquoise","MM.darkolivegreen","MM.darkgrey","MM.cyan","MM.brown","MM.blue")]
liver_hubBlue_moduleMembership$Genes <- rownames(liver_hubBlue_moduleMembership)
liver_hubBlue_moduleMembership$Module <- "Blue"

liver_hub_moduleMembership <- rbind(liver_hubTan_moduleMembership,liver_hubSalmon_moduleMembership,liver_hubRed_moduleMembership,liver_hubDarkturquoise_moduleMembership,liver_hubDarkolivegreen_moduleMembership,liver_hubDarkgrey_moduleMembership,liver_hubCyan_moduleMembership,liver_hubBrown_moduleMembership,liver_hubBlue_moduleMembership)
write.table(liver_hub_moduleMembership, file="liver_hub_moduleMembership.txt", row.names = F, quote=FALSE, sep="\t")
```

### Pancreas

There are no modules significantly associated with age nor with sex.

# 4. Upset plots

## Fig 4/Fig S4

### Brain

``` r
# Import gene lists
brain_table_sig = read.csv(file = "brain_table_sig.csv", header = T); #top dynamic genes = top trendy
rownames(brain_table_sig) <- brain_table_sig$X
brain_table_sig$X <- NULL
sig_brain <- brain_table_sig$Feature

# Create list with sets of interest (i.e. top trendy, module genes and hub genes)
# TAN module
brain_trendyTan <- list(sig_brain, brain_allTan$SYMBOL, brain_hubTan)
names(brain_trendyTan) <- c('Trendy genes','Module genes','Hub genes')

# GREY60 module
brain_trendyGrey60 <- list(sig_brain, brain_allGrey60$SYMBOL, brain_hubGrey60)
names(brain_trendyGrey60) <- c('Trendy genes','Module genes','Hub genes')
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

# insert list of interest here
data=fromList(brain_trendyTan)
# data=fromList(brain_trendyGrey60)
sets = c('Trendy genes','Module genes','Hub genes');
order.by = "freq";
text.scale = c(2.5,2,2,2,2.5,3);
point.size = 8; 
line.size = 1; 
sets.bar.color= c("#74A089","#FDDDA0",  "#F8AFA8") # Trendy-Module-Hub
 

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
  
  # Trendy-Module-Hub
  for(i in 1:3) {
      j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
      if(length(j) > 0) Matrix_layout$color[j] <- c("#74A089","#FDDDA0",  "#F8AFA8")[i]
  }

Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #74A089   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  1 2     1 #74A089   1.0         2yes
    ## 5  2 2     1 #FDDDA0   1.0         2yes
    ## 6  3 2     0  gray83   0.5          6No
    ## 7  1 3     0  gray83   0.5          7No
    ## 8  2 3     1 #FDDDA0   1.0         3yes
    ## 9  3 3     0  gray83   0.5          9No
    ## 10 1 4     1 #74A089   1.0         4yes
    ## 11 2 4     1 #FDDDA0   1.0         4yes
    ## 12 3 4     1 #F8AFA8   1.0         4yes

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-29-1.png)

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-29-2.png)

``` r
# save plot here
tiff('Fig-4-brain-trendy-tan.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-brain-trendy-grey60.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
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
grid.text("Brain Tan Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.613, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.866, y=0.05)
# grid.text("Brain Grey60 Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.783, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.888, y=0.05)
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
### TAN module
list_brain_trendyTan <- overlapGroups(brain_trendyTan)
names(list_brain_trendyTan)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"

``` r
list1_brain_trendyTan <- attr(list_brain_trendyTan, "elements")[list_brain_trendyTan[[1]]]
list2_brain_trendyTan <- attr(list_brain_trendyTan, "elements")[list_brain_trendyTan[[2]]]
list3_brain_trendyTan <- attr(list_brain_trendyTan, "elements")[list_brain_trendyTan[[3]]]
list4_brain_trendyTan <- attr(list_brain_trendyTan, "elements")[list_brain_trendyTan[[4]]]

upsetList_brain_trendyTan <- list(list1_brain_trendyTan,list2_brain_trendyTan,list3_brain_trendyTan,list4_brain_trendyTan)

upsetList_brain_trendyTan <- as.data.frame(sapply(upsetList_brain_trendyTan, '[', seq(max(sapply(upsetList_brain_trendyTan, length)))))
colnames(upsetList_brain_trendyTan) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_brain_trendyTan, file="upsetList_brain_trendyTan.txt", row.names = F, quote=FALSE, sep="\t")

### GREY60 module
list_brain_trendyGrey60 <- overlapGroups(brain_trendyGrey60)
names(list_brain_trendyGrey60)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Module genes:Hub genes"              "Trendy genes:Module genes"          
    ## [5] "Trendy genes:Module genes:Hub genes"

``` r
list1_brain_trendyGrey60 <- attr(list_brain_trendyGrey60, "elements")[list_brain_trendyGrey60[[1]]]
list2_brain_trendyGrey60 <- attr(list_brain_trendyGrey60, "elements")[list_brain_trendyGrey60[[2]]]
list3_brain_trendyGrey60 <- attr(list_brain_trendyGrey60, "elements")[list_brain_trendyGrey60[[3]]]
list4_brain_trendyGrey60 <- attr(list_brain_trendyGrey60, "elements")[list_brain_trendyGrey60[[4]]]
list5_brain_trendyGrey60 <- attr(list_brain_trendyGrey60, "elements")[list_brain_trendyGrey60[[5]]]

upsetList_brain_trendyGrey60 <- list(list1_brain_trendyGrey60,list2_brain_trendyGrey60,list3_brain_trendyGrey60,list4_brain_trendyGrey60,list5_brain_trendyGrey60)

upsetList_brain_trendyGrey60 <- as.data.frame(sapply(upsetList_brain_trendyGrey60, '[', seq(max(sapply(upsetList_brain_trendyGrey60, length)))))
colnames(upsetList_brain_trendyGrey60) <- c("Trendy genes","Module genes","Module genes:Hub genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_brain_trendyGrey60, file="upsetList_brain_trendyGrey60.txt", row.names = F, quote=FALSE, sep="\t")
```

### Heart

``` r
# Import gene lists
heart_table_sig = read.csv(file = "heart_table_sig.csv", header = T); #top dynamic genes = top trendy
rownames(heart_table_sig) <- heart_table_sig$X
heart_table_sig$X <- NULL
sig_heart <- heart_table_sig$Feature

# Create list with sets of interest (i.e. top trendy, module genes and hub genes)
# TAN module
heart_trendyTan <- list(sig_heart, heart_allTan$SYMBOL, heart_hubTan)
names(heart_trendyTan) <- c('Trendy genes','Module genes','Hub genes')

# BLUE module
heart_trendyBlue <- list(sig_heart, heart_allBlue$SYMBOL, heart_hubBlue)
names(heart_trendyBlue) <- c('Trendy genes','Module genes','Hub genes')
```

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

# insert list of interest here
data=fromList(heart_trendyTan)
# data=fromList(heart_trendyBlue)
sets = c('Trendy genes','Module genes','Hub genes');
order.by = "freq";
text.scale = c(2.5,2,2,2,2.5,3);
point.size = 8; 
line.size = 1; 
sets.bar.color= c("#74A089","#FDDDA0","#F8AFA8") # #Trendy-Module-Hub

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

  # Trendy-Module-Hub
  for(i in 1:3) {
      j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
      if(length(j) > 0) Matrix_layout$color[j] <- c("#74A089","#FDDDA0","#F8AFA8")[i]
  }

Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #74A089   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  1 2     1 #74A089   1.0         2yes
    ## 5  2 2     1 #FDDDA0   1.0         2yes
    ## 6  3 2     0  gray83   0.5          6No
    ## 7  1 3     0  gray83   0.5          7No
    ## 8  2 3     1 #FDDDA0   1.0         3yes
    ## 9  3 3     0  gray83   0.5          9No
    ## 10 1 4     1 #74A089   1.0         4yes
    ## 11 2 4     1 #FDDDA0   1.0         4yes
    ## 12 3 4     1 #F8AFA8   1.0         4yes

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-33-1.png)

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-33-2.png)

``` r
# save plot here
tiff('Fig-4-heart-trendy-tan.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-4-heart-trendy-blue.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
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
grid.text("Heart Tan Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.613, y=0.05)
grid.text("*", gp=gpar(fontsize=30), x = 0.866, y=0.05)
# grid.text("Heart Blue Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.6775, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.782, y=0.05)
dev.off()
```

    ## png 
    ##   2

``` r
### TAN module
list_heart_trendyTan <- overlapGroups(heart_trendyTan)
names(list_heart_trendyTan)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"

``` r
list1_heart_trendyTan <- attr(list_heart_trendyTan, "elements")[list_heart_trendyTan[[1]]]
list2_heart_trendyTan <- attr(list_heart_trendyTan, "elements")[list_heart_trendyTan[[2]]]
list3_heart_trendyTan <- attr(list_heart_trendyTan, "elements")[list_heart_trendyTan[[3]]]
list4_heart_trendyTan <- attr(list_heart_trendyTan, "elements")[list_heart_trendyTan[[4]]]

upsetList_heart_trendyTan <- list(list1_heart_trendyTan,list2_heart_trendyTan,list3_heart_trendyTan,list4_heart_trendyTan)

upsetList_heart_trendyTan <- as.data.frame(sapply(upsetList_heart_trendyTan, '[', seq(max(sapply(upsetList_heart_trendyTan, length)))))
colnames(upsetList_heart_trendyTan) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_heart_trendyTan, file="upsetList_heart_trendyTan.txt", row.names = F, quote=FALSE, sep="\t")


### BLUE module
list_heart_trendyBlue <- overlapGroups(heart_trendyBlue)
names(list_heart_trendyBlue)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_heart_trendyBlue <- attr(list_heart_trendyBlue, "elements")[list_heart_trendyBlue[[1]]]
list2_heart_trendyBlue <- attr(list_heart_trendyBlue, "elements")[list_heart_trendyBlue[[2]]]
list3_heart_trendyBlue <- attr(list_heart_trendyBlue, "elements")[list_heart_trendyBlue[[3]]]
list4_heart_trendyBlue <- attr(list_heart_trendyBlue, "elements")[list_heart_trendyBlue[[4]]]
list5_heart_trendyBlue <- attr(list_heart_trendyBlue, "elements")[list_heart_trendyBlue[[5]]]

upsetList_heart_trendyBlue <- list(list1_heart_trendyBlue,list2_heart_trendyBlue,list3_heart_trendyBlue,list4_heart_trendyBlue,list5_heart_trendyBlue)

upsetList_heart_trendyBlue <- as.data.frame(sapply(upsetList_heart_trendyBlue, '[', seq(max(sapply(upsetList_heart_trendyBlue, length)))))
colnames(upsetList_heart_trendyBlue) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_heart_trendyBlue, file="upsetList_heart_trendyBlue.txt", row.names = F, quote=FALSE, sep="\t")
```

### Muscle

``` r
# Import gene lists
muscle_table_sig = read.csv(file = "muscle_table_sig.csv", header = T); #top dynamic genes = top trendy
rownames(muscle_table_sig) <- muscle_table_sig$X
muscle_table_sig$X <- NULL
sig_muscle <- muscle_table_sig$Feature

# Create list with sets of interest (i.e. top trendy, module genes and hub genes)
# MAGENTA module
muscle_trendyMagenta <- list(sig_muscle, muscle_allMagenta$SYMBOL, muscle_hubMagenta)
names(muscle_trendyMagenta) <- c('Trendy genes','Module genes','Hub genes')

# BROWN module
muscle_trendyBrown <- list(sig_muscle, muscle_allBrown$SYMBOL, muscle_hubBrown)
names(muscle_trendyBrown) <- c('Trendy genes','Module genes','Hub genes')

# RED module
muscle_trendyRed <- list(sig_muscle, muscle_allRed$SYMBOL, muscle_hubRed)
names(muscle_trendyRed) <- c('Trendy genes','Module genes','Hub genes')

# PURPLE module
muscle_trendyPurple <- list(sig_muscle, muscle_allPurple$SYMBOL, muscle_hubPurple)
names(muscle_trendyPurple) <- c('Trendy genes','Module genes','Hub genes')

# GREENYELLOW module
muscle_trendyGreenyellow <- list(sig_muscle, muscle_allGreenyellow$SYMBOL, muscle_hubGreenyellow)
names(muscle_trendyGreenyellow) <- c('Trendy genes','Module genes','Hub genes')

# BLUE module
muscle_trendyBlue <- list(sig_muscle, muscle_allBlue$SYMBOL, muscle_hubBlue)
names(muscle_trendyBlue) <- c('Trendy genes','Module genes','Hub genes')
```

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

# insert list of interest here
data=fromList(muscle_trendyMagenta)
# data=fromList(muscle_trendyBrown)
# data=fromList(muscle_trendyRed)
# data=fromList(muscle_trendyPurple)
# data=fromList(muscle_trendyGreenyellow)
# data=fromList(muscle_trendyBlue)
sets = c('Trendy genes','Module genes','Hub genes');
order.by = "freq";
text.scale = c(2.5,2,2,2,2.5,3);
point.size = 8; 
line.size = 1; 
sets.bar.color= c("#74A089","#FDDDA0", "#F8AFA8") #Trendy-Module-Hub
 

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
 
  # Trendy-Module-Hub
   for(i in 1:3) {
       j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
       if(length(j) > 0) Matrix_layout$color[j] <- c("#74A089","#FDDDA0", "#F8AFA8")[i]
   }

Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #74A089   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  1 2     0  gray83   0.5          4No
    ## 5  2 2     1 #FDDDA0   1.0         2yes
    ## 6  3 2     0  gray83   0.5          6No
    ## 7  1 3     1 #74A089   1.0         3yes
    ## 8  2 3     1 #FDDDA0   1.0         3yes
    ## 9  3 3     0  gray83   0.5          9No
    ## 10 1 4     1 #74A089   1.0         4yes
    ## 11 2 4     1 #FDDDA0   1.0         4yes
    ## 12 3 4     1 #F8AFA8   1.0         4yes

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-36-1.png)

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-36-2.png)

``` r
# save plot here
tiff('Fig-4-muscle-trendy-magenta.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-4-muscle-trendy-brown.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-muscle-trendy-red.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-muscle-trendy-purple.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-muscle-trendy-greenyellow.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-muscle-trendy-blue.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
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
grid.text("Muscle Magenta Module", gp=gpar(fontsize=25), x = 0.70, y=0.97)
grid.text("*", gp=gpar(fontsize=30), x = 0.866, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.738, y=0.05)
# grid.text("Muscle Brown Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.675, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.783, y=0.05)
# grid.text("Muscle Red Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.784, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.675, y=0.05)
# grid.text("Muscle Purple Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.678, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.782, y=0.05)
# grid.text("Muscle Greenyellow Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.887, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.678, y=0.05)
# grid.text("Muscle Blue Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.887, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.783, y=0.05)
dev.off()
```

    ## png 
    ##   2

``` r
### MAGENTA module
list_muscle_trendyMagenta <- overlapGroups(muscle_trendyMagenta)
names(list_muscle_trendyMagenta)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"

``` r
list1_muscle_trendyMagenta <- attr(list_muscle_trendyMagenta, "elements")[list_muscle_trendyMagenta[[1]]]
list2_muscle_trendyMagenta <- attr(list_muscle_trendyMagenta, "elements")[list_muscle_trendyMagenta[[2]]]
list3_muscle_trendyMagenta <- attr(list_muscle_trendyMagenta, "elements")[list_muscle_trendyMagenta[[3]]]
list4_muscle_trendyMagenta <- attr(list_muscle_trendyMagenta, "elements")[list_muscle_trendyMagenta[[4]]]

upsetList_muscle_trendyMagenta <- list(list1_muscle_trendyMagenta,list2_muscle_trendyMagenta,list3_muscle_trendyMagenta,list4_muscle_trendyMagenta)

upsetList_muscle_trendyMagenta <- as.data.frame(sapply(upsetList_muscle_trendyMagenta, '[', seq(max(sapply(upsetList_muscle_trendyMagenta, length)))))
colnames(upsetList_muscle_trendyMagenta) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_muscle_trendyMagenta, file="upsetList_muscle_trendyMagenta.txt", row.names = F, quote=FALSE, sep="\t")


### BROWN module
list_muscle_trendyBrown <- overlapGroups(muscle_trendyBrown)
names(list_muscle_trendyBrown)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_muscle_trendyBrown <- attr(list_muscle_trendyBrown, "elements")[list_muscle_trendyBrown[[1]]]
list2_muscle_trendyBrown <- attr(list_muscle_trendyBrown, "elements")[list_muscle_trendyBrown[[2]]]
list3_muscle_trendyBrown <- attr(list_muscle_trendyBrown, "elements")[list_muscle_trendyBrown[[3]]]
list4_muscle_trendyBrown <- attr(list_muscle_trendyBrown, "elements")[list_muscle_trendyBrown[[4]]]
list5_muscle_trendyBrown <- attr(list_muscle_trendyBrown, "elements")[list_muscle_trendyBrown[[5]]]

upsetList_muscle_trendyBrown <- list(list1_muscle_trendyBrown,list2_muscle_trendyBrown,list3_muscle_trendyBrown,list4_muscle_trendyBrown,list5_muscle_trendyBrown)

upsetList_muscle_trendyBrown <- as.data.frame(sapply(upsetList_muscle_trendyBrown, '[', seq(max(sapply(upsetList_muscle_trendyBrown, length)))))
colnames(upsetList_muscle_trendyBrown) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_muscle_trendyBrown, file="upsetList_muscle_trendyBrown.txt", row.names = F, quote=FALSE, sep="\t")

### RED module
list_muscle_trendyRed <- overlapGroups(muscle_trendyRed)
names(list_muscle_trendyRed)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_muscle_trendyRed <- attr(list_muscle_trendyRed, "elements")[list_muscle_trendyRed[[1]]]
list2_muscle_trendyRed <- attr(list_muscle_trendyRed, "elements")[list_muscle_trendyRed[[2]]]
list3_muscle_trendyRed <- attr(list_muscle_trendyRed, "elements")[list_muscle_trendyRed[[3]]]
list4_muscle_trendyRed <- attr(list_muscle_trendyRed, "elements")[list_muscle_trendyRed[[4]]]
list5_muscle_trendyRed <- attr(list_muscle_trendyRed, "elements")[list_muscle_trendyRed[[5]]]

upsetList_muscle_trendyRed <- list(list1_muscle_trendyRed,list2_muscle_trendyRed,list3_muscle_trendyRed,list4_muscle_trendyRed,list5_muscle_trendyRed)

upsetList_muscle_trendyRed <- as.data.frame(sapply(upsetList_muscle_trendyRed, '[', seq(max(sapply(upsetList_muscle_trendyRed, length)))))
colnames(upsetList_muscle_trendyRed) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_muscle_trendyRed, file="upsetList_muscle_trendyRed.txt", row.names = F, quote=FALSE, sep="\t")

### PURPLE module
list_muscle_trendyPurple <- overlapGroups(muscle_trendyPurple)
names(list_muscle_trendyPurple)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_muscle_trendyPurple <- attr(list_muscle_trendyPurple, "elements")[list_muscle_trendyPurple[[1]]]
list2_muscle_trendyPurple <- attr(list_muscle_trendyPurple, "elements")[list_muscle_trendyPurple[[2]]]
list3_muscle_trendyPurple <- attr(list_muscle_trendyPurple, "elements")[list_muscle_trendyPurple[[3]]]
list4_muscle_trendyPurple <- attr(list_muscle_trendyPurple, "elements")[list_muscle_trendyPurple[[4]]]
list5_muscle_trendyPurple <- attr(list_muscle_trendyPurple, "elements")[list_muscle_trendyPurple[[5]]]

upsetList_muscle_trendyPurple <- list(list1_muscle_trendyPurple,list2_muscle_trendyPurple,list3_muscle_trendyPurple,list4_muscle_trendyPurple,list5_muscle_trendyPurple)

upsetList_muscle_trendyPurple <- as.data.frame(sapply(upsetList_muscle_trendyPurple, '[', seq(max(sapply(upsetList_muscle_trendyPurple, length)))))
colnames(upsetList_muscle_trendyPurple) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes" )

write.table(upsetList_muscle_trendyPurple, file="upsetList_muscle_trendyPurple.txt", row.names = F, quote=FALSE, sep="\t")

### GREENYELLOW module
list_muscle_trendyGreenyellow <- overlapGroups(muscle_trendyGreenyellow)
names(list_muscle_trendyGreenyellow)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Module genes:Hub genes"             
    ## [5] "Trendy genes:Module genes:Hub genes"

``` r
list1_muscle_trendyGreenyellow <- attr(list_muscle_trendyGreenyellow, "elements")[list_muscle_trendyGreenyellow[[1]]]
list2_muscle_trendyGreenyellow <- attr(list_muscle_trendyGreenyellow, "elements")[list_muscle_trendyGreenyellow[[2]]]
list3_muscle_trendyGreenyellow <- attr(list_muscle_trendyGreenyellow, "elements")[list_muscle_trendyGreenyellow[[3]]]
list4_muscle_trendyGreenyellow <- attr(list_muscle_trendyGreenyellow, "elements")[list_muscle_trendyGreenyellow[[4]]]
list5_muscle_trendyGreenyellow <- attr(list_muscle_trendyGreenyellow, "elements")[list_muscle_trendyGreenyellow[[5]]]

upsetList_muscle_trendyGreenyellow <- list(list1_muscle_trendyGreenyellow,list2_muscle_trendyGreenyellow,list3_muscle_trendyGreenyellow,list4_muscle_trendyGreenyellow,list5_muscle_trendyGreenyellow)

upsetList_muscle_trendyGreenyellow <- as.data.frame(sapply(upsetList_muscle_trendyGreenyellow, '[', seq(max(sapply(upsetList_muscle_trendyGreenyellow, length)))))
colnames(upsetList_muscle_trendyGreenyellow) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Module genes:Hub genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_muscle_trendyGreenyellow, file="upsetList_muscle_trendyGreenyellow.txt", row.names = F, quote=FALSE, sep="\t")

### BLUE module
list_muscle_trendyBlue <- overlapGroups(muscle_trendyBlue)
names(list_muscle_trendyBlue)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Module genes:Hub genes"              "Trendy genes:Module genes"          
    ## [5] "Trendy genes:Module genes:Hub genes"

``` r
list1_muscle_trendyBlue <- attr(list_muscle_trendyBlue, "elements")[list_muscle_trendyBlue[[1]]]
list2_muscle_trendyBlue <- attr(list_muscle_trendyBlue, "elements")[list_muscle_trendyBlue[[2]]]
list3_muscle_trendyBlue <- attr(list_muscle_trendyBlue, "elements")[list_muscle_trendyBlue[[3]]]
list4_muscle_trendyBlue <- attr(list_muscle_trendyBlue, "elements")[list_muscle_trendyBlue[[4]]]
list5_muscle_trendyBlue <- attr(list_muscle_trendyBlue, "elements")[list_muscle_trendyBlue[[5]]]

upsetList_muscle_trendyBlue <- list(list1_muscle_trendyBlue,list2_muscle_trendyBlue,list3_muscle_trendyBlue,list4_muscle_trendyBlue,list5_muscle_trendyBlue)

upsetList_muscle_trendyBlue <- as.data.frame(sapply(upsetList_muscle_trendyBlue, '[', seq(max(sapply(upsetList_muscle_trendyBlue, length)))))
colnames(upsetList_muscle_trendyBlue) <- c("Trendy genes","Module genes","Module genes:Hub genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_muscle_trendyBlue, file="upsetList_muscle_trendyBlue.txt", row.names = F, quote=FALSE, sep="\t")
```

### Liver

``` r
# Import gene lists
liver_table_sig = read.csv(file = "liver_table_sig.csv", header = T); #top dynamic genes = top trendy
rownames(liver_table_sig) <- liver_table_sig$X
liver_table_sig$X <- NULL
sig_liver <- liver_table_sig$Feature

# Create list with sets of interest (i.e. top trendy, module genes and hub genes)
# SALMON module
liver_trendySalmon <- list(sig_liver, liver_allSalmon$SYMBOL, liver_hubSalmon)
names(liver_trendySalmon) <- c('Trendy genes','Module genes','Hub genes')

# DARKTURQUOISE module
liver_trendyDarkturquoise <- list(sig_liver, liver_allDarkturquoise$SYMBOL, liver_hubDarkturquoise)
names(liver_trendyDarkturquoise) <- c('Trendy genes','Module genes','Hub genes')

# TAN module
liver_trendyTan <- list(sig_liver, liver_allTan$SYMBOL, liver_hubTan)
names(liver_trendyTan) <- c('Trendy genes','Module genes','Hub genes')

# RED module
liver_trendyRed <- list(sig_liver, liver_allRed$SYMBOL, liver_hubRed)
names(liver_trendyRed) <- c('Trendy genes','Module genes','Hub genes')

# DARKOLIVEGREEN module
liver_trendyDarkolivegreen <- list(sig_liver, liver_allDarkolivegreen$SYMBOL, liver_hubDarkolivegreen)
names(liver_trendyDarkolivegreen) <- c('Trendy genes','Module genes','Hub genes')

# DARKGREY module
liver_trendyDarkgrey <- list(sig_liver, liver_allDarkgrey$SYMBOL, liver_hubDarkgrey)
names(liver_trendyDarkgrey) <- c('Trendy genes','Module genes','Hub genes')

# CYAN module
liver_trendyCyan <- list(sig_liver, liver_allCyan$SYMBOL, liver_hubCyan)
names(liver_trendyCyan) <- c('Trendy genes','Module genes','Hub genes')

# BROWN module
liver_trendyBrown <- list(sig_liver, liver_allBrown$SYMBOL, liver_hubBrown)
names(liver_trendyBrown) <- c('Trendy genes','Module genes','Hub genes')

# BLUE module
liver_trendyBlue <- list(sig_liver, liver_allBlue$SYMBOL, liver_hubBlue)
names(liver_trendyBlue) <- c('Trendy genes','Module genes','Hub genes')
```

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

# insert list of interest here
data=fromList(liver_trendySalmon);
# data=fromList(liver_trendyDarkturquoise)
# data=fromList(liver_trendyTan)
# data=fromList(liver_trendyRed)
# data=fromList(liver_trendyDarkolivegreen)
# data=fromList(liver_trendyDarkgrey)
# data=fromList(liver_trendyCyan)
# data=fromList(liver_trendyBrown)
# data=fromList(liver_trendyBlue)
sets = c('Trendy genes','Module genes','Hub genes');
order.by = "freq";
text.scale = c(2.5,2,2,2,2.5,3);
point.size = 8; 
line.size = 1; 
sets.bar.color= c("#74A089","#FDDDA0","#F8AFA8") #Trendy-Module-Hub
 

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
  
  # Trendy-Module-Hub
   for(i in 1:3) {
      j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
      if(length(j) > 0) Matrix_layout$color[j] <- c("#74A089","#FDDDA0", "#F8AFA8")[i]
   }

Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #74A089   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  1 2     1 #74A089   1.0         2yes
    ## 5  2 2     1 #FDDDA0   1.0         2yes
    ## 6  3 2     0  gray83   0.5          6No
    ## 7  1 3     0  gray83   0.5          7No
    ## 8  2 3     1 #FDDDA0   1.0         3yes
    ## 9  3 3     0  gray83   0.5          9No
    ## 10 1 4     1 #74A089   1.0         4yes
    ## 11 2 4     1 #FDDDA0   1.0         4yes
    ## 12 3 4     1 #F8AFA8   1.0         4yes

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-39-1.png)

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-39-2.png)

``` r
# save plot here
tiff('Fig-4-liver-trendy-salmon.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-4-liver-trendy-darkturquoise.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-tan.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-red.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-darkolivegreen.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-darkgrey.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-cyan.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-brown.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
# tiff('Fig-S4-liver-trendy-blue.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
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
grid.text("Liver Salmon Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
grid.text("*", gp=gpar(fontsize=30), x = 0.865, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.615, y=0.05)
# grid.text("Liver Darkturquoise Module", gp=gpar(fontsize=25), x = 0.70, y=0.98)
# grid.text("*", gp=gpar(fontsize=30), x = 0.677, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.782, y=0.05)
# grid.text("Liver Tan Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.57, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.78, y=0.05)
# grid.text("Liver Red Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.675, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.78, y=0.05)
# grid.text("Liver Darkolivegreen Module", gp=gpar(fontsize=25), x = 0.71, y=0.97)
# grid.text("*", gp=gpar(fontsize=30), x = 0.675, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.885, y=0.05)
# grid.text("Liver Darkgrey Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.784, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.573, y=0.05)
# grid.text("Liver Cyan Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.57, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.78, y=0.05)
# grid.text("Liver Brown Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.785, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.676, y=0.05)
# grid.text("Liver Blue Module", gp=gpar(fontsize=25), x = 0.70, y=0.95)
# grid.text("*", gp=gpar(fontsize=30), x = 0.887, y=0.05)
# grid.text("*", gp=gpar(fontsize=30), x = 0.675, y=0.05)
dev.off()
```

    ## png 
    ##   2

``` r
### SALMON module
list_liver_trendySalmon <- overlapGroups(liver_trendySalmon)
names(list_liver_trendySalmon)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"

``` r
list1_liver_trendySalmon <- attr(list_liver_trendySalmon, "elements")[list_liver_trendySalmon[[1]]]
list2_liver_trendySalmon <- attr(list_liver_trendySalmon, "elements")[list_liver_trendySalmon[[2]]]
list3_liver_trendySalmon <- attr(list_liver_trendySalmon, "elements")[list_liver_trendySalmon[[3]]]
list4_liver_trendySalmon <- attr(list_liver_trendySalmon, "elements")[list_liver_trendySalmon[[4]]]

upsetList_liver_trendySalmon <- list(list1_liver_trendySalmon,list2_liver_trendySalmon,list3_liver_trendySalmon,list4_liver_trendySalmon)

upsetList_liver_trendySalmon <- as.data.frame(sapply(upsetList_liver_trendySalmon, '[', seq(max(sapply(upsetList_liver_trendySalmon, length)))))
colnames(upsetList_liver_trendySalmon) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_liver_trendySalmon, file="upsetList_liver_trendySalmon.txt", row.names = F, quote=FALSE, sep="\t")

### DARKTURQUOISE module
list_liver_trendyDarkturquoise <- overlapGroups(liver_trendyDarkturquoise)
names(list_liver_trendyDarkturquoise)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyDarkturquoise <- attr(list_liver_trendyDarkturquoise, "elements")[list_liver_trendyDarkturquoise[[1]]]
list2_liver_trendyDarkturquoise <- attr(list_liver_trendyDarkturquoise, "elements")[list_liver_trendyDarkturquoise[[2]]]
list3_liver_trendyDarkturquoise <- attr(list_liver_trendyDarkturquoise, "elements")[list_liver_trendyDarkturquoise[[3]]]
list4_liver_trendyDarkturquoise <- attr(list_liver_trendyDarkturquoise, "elements")[list_liver_trendyDarkturquoise[[4]]]
list5_liver_trendyDarkturquoise <- attr(list_liver_trendyDarkturquoise, "elements")[list_liver_trendyDarkturquoise[[5]]]

upsetList_liver_trendyDarkturquoise <- list(list1_liver_trendyDarkturquoise,list2_liver_trendyDarkturquoise,list3_liver_trendyDarkturquoise,list4_liver_trendyDarkturquoise,list5_liver_trendyDarkturquoise)

upsetList_liver_trendyDarkturquoise <- as.data.frame(sapply(upsetList_liver_trendyDarkturquoise, '[', seq(max(sapply(upsetList_liver_trendyDarkturquoise, length)))))
colnames(upsetList_liver_trendyDarkturquoise) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyDarkturquoise, file="upsetList_liver_trendyDarkturquoise.txt", row.names = F, quote=FALSE, sep="\t")

### TAN module
list_liver_trendyTan <- overlapGroups(liver_trendyTan)
names(list_liver_trendyTan)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyTan <- attr(list_liver_trendyTan, "elements")[list_liver_trendyTan[[1]]]
list2_liver_trendyTan <- attr(list_liver_trendyTan, "elements")[list_liver_trendyTan[[2]]]
list3_liver_trendyTan <- attr(list_liver_trendyTan, "elements")[list_liver_trendyTan[[3]]]
list4_liver_trendyTan <- attr(list_liver_trendyTan, "elements")[list_liver_trendyTan[[4]]]
list5_liver_trendyTan <- attr(list_liver_trendyTan, "elements")[list_liver_trendyTan[[5]]]

upsetList_liver_trendyTan <- list(list1_liver_trendyTan,list2_liver_trendyTan,list3_liver_trendyTan,list4_liver_trendyTan,list5_liver_trendyTan)

upsetList_liver_trendyTan <- as.data.frame(sapply(upsetList_liver_trendyTan, '[', seq(max(sapply(upsetList_liver_trendyTan, length)))))
colnames(upsetList_liver_trendyTan) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyTan, file="upsetList_liver_trendyTan.txt", row.names = F, quote=FALSE, sep="\t")

### RED module
list_liver_trendyRed <- overlapGroups(liver_trendyRed)
names(list_liver_trendyRed)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyRed <- attr(list_liver_trendyRed, "elements")[list_liver_trendyRed[[1]]]
list2_liver_trendyRed <- attr(list_liver_trendyRed, "elements")[list_liver_trendyRed[[2]]]
list3_liver_trendyRed <- attr(list_liver_trendyRed, "elements")[list_liver_trendyRed[[3]]]
list4_liver_trendyRed <- attr(list_liver_trendyRed, "elements")[list_liver_trendyRed[[4]]]
list5_liver_trendyRed <- attr(list_liver_trendyRed, "elements")[list_liver_trendyRed[[5]]]

upsetList_liver_trendyRed <- list(list1_liver_trendyRed,list2_liver_trendyRed,list3_liver_trendyRed,list4_liver_trendyRed,list5_liver_trendyRed)

upsetList_liver_trendyRed <- as.data.frame(sapply(upsetList_liver_trendyRed, '[', seq(max(sapply(upsetList_liver_trendyRed, length)))))
colnames(upsetList_liver_trendyRed) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyRed, file="upsetList_liver_trendyRed.txt", row.names = F, quote=FALSE, sep="\t")

### DARKOLIVEGREEN module
list_liver_trendyDarkolivegreen <- overlapGroups(liver_trendyDarkolivegreen)
names(list_liver_trendyDarkolivegreen)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Module genes:Hub genes"             
    ## [5] "Trendy genes:Module genes:Hub genes"

``` r
list1_liver_trendyDarkolivegreen <- attr(list_liver_trendyDarkolivegreen, "elements")[list_liver_trendyDarkolivegreen[[1]]]
list2_liver_trendyDarkolivegreen <- attr(list_liver_trendyDarkolivegreen, "elements")[list_liver_trendyDarkolivegreen[[2]]]
list3_liver_trendyDarkolivegreen <- attr(list_liver_trendyDarkolivegreen, "elements")[list_liver_trendyDarkolivegreen[[3]]]
list4_liver_trendyDarkolivegreen <- attr(list_liver_trendyDarkolivegreen, "elements")[list_liver_trendyDarkolivegreen[[4]]]
list5_liver_trendyDarkolivegreen <- attr(list_liver_trendyDarkolivegreen, "elements")[list_liver_trendyDarkolivegreen[[5]]]

upsetList_liver_trendyDarkolivegreen <- list(list1_liver_trendyDarkolivegreen,list2_liver_trendyDarkolivegreen,list3_liver_trendyDarkolivegreen,list4_liver_trendyDarkolivegreen,list5_liver_trendyDarkolivegreen)

upsetList_liver_trendyDarkolivegreen <- as.data.frame(sapply(upsetList_liver_trendyDarkolivegreen, '[', seq(max(sapply(upsetList_liver_trendyDarkolivegreen, length)))))
colnames(upsetList_liver_trendyDarkolivegreen) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Module genes:Hub genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_liver_trendyDarkolivegreen, file="upsetList_liver_trendyDarkolivegreen.txt", row.names = F, quote=FALSE, sep="\t")

### DARKGREY module
list_liver_trendyDarkgrey <- overlapGroups(liver_trendyDarkgrey)
names(list_liver_trendyDarkgrey)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyDarkgrey <- attr(list_liver_trendyDarkgrey, "elements")[list_liver_trendyDarkgrey[[1]]]
list2_liver_trendyDarkgrey <- attr(list_liver_trendyDarkgrey, "elements")[list_liver_trendyDarkgrey[[2]]]
list3_liver_trendyDarkgrey <- attr(list_liver_trendyDarkgrey, "elements")[list_liver_trendyDarkgrey[[3]]]
list4_liver_trendyDarkgrey <- attr(list_liver_trendyDarkgrey, "elements")[list_liver_trendyDarkgrey[[4]]]
list5_liver_trendyDarkgrey <- attr(list_liver_trendyDarkgrey, "elements")[list_liver_trendyDarkgrey[[5]]]

upsetList_liver_trendyDarkgrey <- list(list1_liver_trendyDarkgrey,list2_liver_trendyDarkgrey,list3_liver_trendyDarkgrey,list4_liver_trendyDarkgrey,list5_liver_trendyDarkgrey)

upsetList_liver_trendyDarkgrey <- as.data.frame(sapply(upsetList_liver_trendyDarkgrey, '[', seq(max(sapply(upsetList_liver_trendyDarkgrey, length)))))
colnames(upsetList_liver_trendyDarkgrey) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyDarkgrey, file="upsetList_liver_trendyDarkgrey.txt", row.names = F, quote=FALSE, sep="\t")

### CYAN module
list_liver_trendyCyan <- overlapGroups(liver_trendyCyan)
names(list_liver_trendyCyan)
```

    ## [1] "Trendy genes"                        "Trendy genes:Module genes"          
    ## [3] "Module genes"                        "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyCyan <- attr(list_liver_trendyCyan, "elements")[list_liver_trendyCyan[[1]]]
list2_liver_trendyCyan <- attr(list_liver_trendyCyan, "elements")[list_liver_trendyCyan[[2]]]
list3_liver_trendyCyan <- attr(list_liver_trendyCyan, "elements")[list_liver_trendyCyan[[3]]]
list4_liver_trendyCyan <- attr(list_liver_trendyCyan, "elements")[list_liver_trendyCyan[[4]]]
list5_liver_trendyCyan <- attr(list_liver_trendyCyan, "elements")[list_liver_trendyCyan[[5]]]

upsetList_liver_trendyCyan <- list(list1_liver_trendyCyan,list2_liver_trendyCyan,list3_liver_trendyCyan,list4_liver_trendyCyan,list5_liver_trendyCyan)

upsetList_liver_trendyCyan <- as.data.frame(sapply(upsetList_liver_trendyCyan, '[', seq(max(sapply(upsetList_liver_trendyCyan, length)))))
colnames(upsetList_liver_trendyCyan) <- c("Trendy genes","Trendy genes:Module genes","Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyCyan, file="upsetList_liver_trendyCyan.txt", row.names = F, quote=FALSE, sep="\t")

### BROWN module
list_liver_trendyBrown <- overlapGroups(liver_trendyBrown)
names(list_liver_trendyBrown)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Trendy genes:Module genes:Hub genes"
    ## [5] "Module genes:Hub genes"

``` r
list1_liver_trendyBrown <- attr(list_liver_trendyBrown, "elements")[list_liver_trendyBrown[[1]]]
list2_liver_trendyBrown <- attr(list_liver_trendyBrown, "elements")[list_liver_trendyBrown[[2]]]
list3_liver_trendyBrown <- attr(list_liver_trendyBrown, "elements")[list_liver_trendyBrown[[3]]]
list4_liver_trendyBrown <- attr(list_liver_trendyBrown, "elements")[list_liver_trendyBrown[[4]]]
list5_liver_trendyBrown <- attr(list_liver_trendyBrown, "elements")[list_liver_trendyBrown[[5]]]

upsetList_liver_trendyBrown <- list(list1_liver_trendyBrown,list2_liver_trendyBrown,list3_liver_trendyBrown,list4_liver_trendyBrown,list5_liver_trendyBrown)

upsetList_liver_trendyBrown <- as.data.frame(sapply(upsetList_liver_trendyBrown, '[', seq(max(sapply(upsetList_liver_trendyBrown, length)))))
colnames(upsetList_liver_trendyBrown) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Trendy genes:Module genes:Hub genes","Module genes:Hub genes")

write.table(upsetList_liver_trendyBrown, file="upsetList_liver_trendyBrown.txt", row.names = F, quote=FALSE, sep="\t")

### BLUE module
list_liver_trendyBlue <- overlapGroups(liver_trendyBlue)
names(list_liver_trendyBlue)
```

    ## [1] "Trendy genes"                        "Module genes"                       
    ## [3] "Trendy genes:Module genes"           "Module genes:Hub genes"             
    ## [5] "Trendy genes:Module genes:Hub genes"

``` r
list1_liver_trendyBlue <- attr(list_liver_trendyBlue, "elements")[list_liver_trendyBlue[[1]]]
list2_liver_trendyBlue <- attr(list_liver_trendyBlue, "elements")[list_liver_trendyBlue[[2]]]
list3_liver_trendyBlue <- attr(list_liver_trendyBlue, "elements")[list_liver_trendyBlue[[3]]]
list4_liver_trendyBlue <- attr(list_liver_trendyBlue, "elements")[list_liver_trendyBlue[[4]]]
list5_liver_trendyBlue <- attr(list_liver_trendyBlue, "elements")[list_liver_trendyBlue[[5]]]

upsetList_liver_trendyBlue <- list(list1_liver_trendyBlue,list2_liver_trendyBlue,list3_liver_trendyBlue,list4_liver_trendyBlue,list5_liver_trendyBlue)

upsetList_liver_trendyBlue <- as.data.frame(sapply(upsetList_liver_trendyBlue, '[', seq(max(sapply(upsetList_liver_trendyBlue, length)))))
colnames(upsetList_liver_trendyBlue) <- c("Trendy genes","Module genes","Trendy genes:Module genes","Module genes:Hub genes","Trendy genes:Module genes:Hub genes")

write.table(upsetList_liver_trendyBlue, file="upsetList_liver_trendyBlue.txt", row.names = F, quote=FALSE, sep="\t")
```

## Fig 5 - left

``` r
# create list comprising each tissue's TRENDY-MODULE overlap; in tissues with more than one module, we merged the vectors of each module and extracted the unique elements 


allTissues_hub <- list(
  na.omit(unique(c(upsetList_brain_trendyTan$`Trendy genes:Module genes:Hub genes`))),
  
  na.omit(unique(c(upsetList_heart_trendyTan$`Trendy genes:Module genes:Hub genes`,upsetList_heart_trendyBlue$`Trendy genes:Module genes:Hub genes`))),
  
  na.omit(unique(c(upsetList_liver_trendySalmon$`Trendy genes:Module genes:Hub genes`,
                   upsetList_liver_trendyDarkturquoise$`Trendy genes:Module genes:Hub genes`,
                   upsetList_liver_trendyTan$`Trendy genes:Module genes:Hub genes`,
                   upsetList_liver_trendyBlue$`Trendy genes:Module genes:Hub genes`,
                   upsetList_liver_trendyCyan$`Trendy genes:Module genes:Hub genes`,
                   upsetList_liver_trendyDarkgrey$`Trendy genes:Module genes:Hub genes`))),
  
  na.omit(unique(c(upsetList_muscle_trendyMagenta$`Trendy genes:Module genes:Hub genes`,
                   upsetList_muscle_trendyBrown$`Trendy genes:Module genes:Hub genes`,
                   upsetList_muscle_trendyBlue$`Trendy genes:Module genes:Hub genes`)))
  )
names(allTissues_hub) <- c('Brain','Heart','Liver','Muscle')
```

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
sets.bar.color= c("#EEDD88","#EE8866","#FFAABB","#77AADD")#Liver-Heart-Muscle-Brain


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

  
     #Liver-Heart-Muscle-Brain
   for(i in 1:4) {
       j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
       if(length(j) > 0) Matrix_layout$color[j] <- c("#EEDD88","#EE8866","#FFAABB","#77AADD")[i]
   }
  

  
Matrix_layout
```

    ##    y x value   color alpha Intersection
    ## 1  1 1     1 #EEDD88   1.0         1yes
    ## 2  2 1     0  gray83   0.5          2No
    ## 3  3 1     0  gray83   0.5          3No
    ## 4  4 1     0  gray83   0.5          4No
    ## 5  1 2     0  gray83   0.5          5No
    ## 6  2 2     1 #EE8866   1.0         2yes
    ## 7  3 2     0  gray83   0.5          7No
    ## 8  4 2     0  gray83   0.5          8No
    ## 9  1 3     0  gray83   0.5          9No
    ## 10 2 3     0  gray83   0.5         10No
    ## 11 3 3     1 #FFAABB   1.0         3yes
    ## 12 4 3     0  gray83   0.5         12No
    ## 13 1 4     0  gray83   0.5         13No
    ## 14 2 4     0  gray83   0.5         14No
    ## 15 3 4     0  gray83   0.5         15No
    ## 16 4 4     1 #77AADD   1.0         4yes
    ## 17 1 5     1 #EEDD88   1.0         5yes
    ## 18 2 5     1 #EE8866   1.0         5yes
    ## 19 3 5     0  gray83   0.5         19No
    ## 20 4 5     0  gray83   0.5         20No

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-42-1.png)

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

![](WGCNA_2_files/figure-markdown_github/unnamed-chunk-42-2.png)

``` r
tiff('Fig-5-left_hub.tiff', units="cm", width=30, height=20, res=300, compression = 'lzw')
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
grid.text("Overlap Trendy-Hub genes", gp=gpar(fontsize=25), x = 0.67, y=0.93)
dev.off()
```

    ## png 
    ##   2

``` r
#age - hub
list_allTissues_hub <- overlapGroups(allTissues_hub)
names(list_allTissues_hub)
```

    ## [1] "Liver"       "Heart"       "Muscle"      "Brain"       "Heart:Liver"

``` r
list1_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[1]]]
list2_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[2]]]
list3_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[3]]]
list4_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[4]]]
list5_allTissues_hub <- attr(list_allTissues_hub, "elements")[list_allTissues_hub[[5]]]

upsetList_allTissues_hub <- list(list1_allTissues_hub,list2_allTissues_hub,list3_allTissues_hub,list4_allTissues_hub,list5_allTissues_hub)

upsetList_allTissues_hub <- as.data.frame(sapply(upsetList_allTissues_hub, '[', seq(max(sapply(upsetList_allTissues_hub, length)))))
colnames(upsetList_allTissues_hub) <- c("Liver","Heart","Muscle","Brain","Heart:Liver")

write.table(upsetList_allTissues_hub, file="upsetList_allTissues_hub_genes.txt", row.names = F, quote=FALSE, sep="\t")
```

# Session Info

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
    ## [1] UpSetR_1.4.0          ggplot2_3.3.2         WGCNA_1.69           
    ## [4] fastcluster_1.1.25    dynamicTreeCut_1.63-1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Biobase_2.50.0        bit64_4.0.5           splines_4.0.3        
    ##  [4] foreach_1.5.1         Formula_1.2-4         stats4_4.0.3         
    ##  [7] latticeExtra_0.6-29   blob_1.2.1            yaml_2.2.1           
    ## [10] impute_1.62.0         pillar_1.4.7          RSQLite_2.2.1        
    ## [13] backports_1.2.0       lattice_0.20-41       glue_1.4.2           
    ## [16] digest_0.6.27         RColorBrewer_1.1-2    checkmate_2.0.0      
    ## [19] colorspace_2.0-0      htmltools_0.5.0       preprocessCore_1.50.0
    ## [22] Matrix_1.2-18         plyr_1.8.6            pkgconfig_2.0.3      
    ## [25] purrr_0.3.4           GO.db_3.12.1          scales_1.1.1         
    ## [28] jpeg_0.1-8.1          tibble_3.0.4          htmlTable_2.1.0      
    ## [31] generics_0.1.0        farver_2.0.3          IRanges_2.24.0       
    ## [34] ellipsis_0.3.1        withr_2.3.0           nnet_7.3-14          
    ## [37] BiocGenerics_0.36.0   survival_3.2-7        magrittr_2.0.1       
    ## [40] crayon_1.3.4          memoise_1.1.0         evaluate_0.14        
    ## [43] doParallel_1.0.16     foreign_0.8-80        tools_4.0.3          
    ## [46] data.table_1.13.4     lifecycle_0.2.0       matrixStats_0.57.0   
    ## [49] stringr_1.4.0         S4Vectors_0.28.0      munsell_0.5.0        
    ## [52] cluster_2.1.0         AnnotationDbi_1.52.0  compiler_4.0.3       
    ## [55] rlang_0.4.9           iterators_1.0.13      rstudioapi_0.13      
    ## [58] htmlwidgets_1.5.2     labeling_0.4.2        base64enc_0.1-3      
    ## [61] rmarkdown_2.5         gtable_0.3.0          codetools_0.2-16     
    ## [64] DBI_1.1.0             R6_2.5.0              gridExtra_2.3        
    ## [67] knitr_1.30            dplyr_1.0.2           bit_4.0.4            
    ## [70] Hmisc_4.4-1           stringi_1.5.3         parallel_4.0.3       
    ## [73] Rcpp_1.0.5            vctrs_0.3.5           rpart_4.1-15         
    ## [76] png_0.1-7             tidyselect_1.1.0      xfun_0.19
