![](./logo.png)
------ 
## Table of Content

* [Quick start](#quick-start)
* [Introduction](#introduction)
* [How to install](#how-to-install)
* [Inputs](#inputs)
* [Parameters](#parameters)
* [Data preparation for running coGSEA](#data-preparation-for-running-cogsea)
   * [MicroArray data](#microarray-data)

   * [RNAseq data](#rnaseq-data)


# Quick start

``` r
coGSEA(ElistObject = elist, contrastMatrix = contrast,
ENTREZGenesId = elist$genes$ENTREZ, geneSetCollection = "H",
specie = "Mus musculus", directoryPath = "/path/to/existing/dir")

```
# Introduction
**coGSEA** (**co**mparative **G**ene **S**et **E**nrichment **A**nalysis) is a **R** package developped to help you compare and combine up to 12 different methods of Gene Set Enrichment Analsysis.
In the current version, those 12 methods include :
- [camera](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
- [gage](https://bioconductor.org/packages/release/bioc/html/gage.html)
- [globaltest](https://bioconductor.org/packages/release/bioc/html/globaltest.html)
- [gsva](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
- [ssgsea](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
- [zscore](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
- [plage](https://bioconductor.org/packages/release/bioc/html/GSVA.html)
- [ora](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3134237/)
- [padog](http://bioconductor.org/packages/release/bioc/html/PADOG.html)
- [roast](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
- [safe](http://bioconductor.org/packages/release/bioc/html/safe.html)
- [setrank](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1571-6)

**Disclaimer** : This tool is largely inspired by the [eGSEA](http://bioconductor.org/packages/release/bioc/html/EGSEA.html) R package (and contains some of its code)

# How to install

``` r
devtools::install_github("maxibor/coGSEA",build_vignettes = TRUE)
```



# Inputs
- **ElistObject** : a Elist object. Your Elist object must be either a *limma-voom* object in case of RNAseq data, or have a estimated dispersion in case of MicroArray data. Dispersion estimates can be performed using the `estimateDisp()` function from the *edgeR* package.

- **contrastMatrix** : the contrast matrix object for your experiment. For some help to make a contrast matrix, you can check the vignette of the limma package.

- **ENTREZGenesIds** : the Elist object containing the ENTREZ gene identifiers. For example `elist$genes$ENTREZ`

# Parameters

- **geneSetCollection** : A geneset collection from the [MSigDB](http://software.broadinstitute.org/gsea/msigdb) database. Currently only three are supported : `H` (Hallmark) and two subset of `C2_KEGG` and `C2_REACTOME`(Kegg and Reactome). You can also use your own geneset collection. In this case, you need to provide an list object with genesets and genes. Such an example file can be found [here](./exampleData/geneset.rds)

- **Specie** : The organism from which the data were extracted. Currently only `Homo sapiens` and `Mus musculus` are supported. Default : `Mus musculus`

- **directoryPath** : a path to an existing directory where coGSEA results and plots will be saved. Default : `"./"`  
Example :  `"~/coGSEA_results"`

- **alpha** : Alpha *p* value threshold. Default : `0.05`   
example : `0.05`

- **pvalAdjMethod** : *p* value adjustment method for multiple testing. Default : `BH`.  
To select among the following methods :
	- `holm`
	- `hochberg`
	- `hommel`
	- `bonferroni`
	- `BH` (Benjamini Hochberg)
	- `BY` (Benjamini Yekutieli)
	- `fdr` (False Discovery rate)
	- `none`  


- **pvalCombMethod** : the method to combine all the p.values of one geneset accros the different GSEA methods. Default : `sumlog`   
To select among the following methods :
	- `sumz` (sum z method)
	- `votep` (vote counting method)
	- `minimump` (Wilkinson's method)
	- `sumlog` (Fisher's method) - selected by default
	- `sump` (sum p method)
	- `logitp` (logit method)
	- `meanp` (mean p method)
	- `maximump` (Wilkinson's method)


- **min.intersection.size** : Graphical Parameter to select the minimum size of an intersection of genesets between different methods. Default : `1`  
Example : `2`

- **GSEA.Methods** : between 1 and 12 methods to select from the methods listed in the introduction. Default : ` c("camera", "gage","globaltest", "gsva", "ssgsea", "zscore", "plage", "ora", "padog", "roast","safe")`    
Example : ` c("camera", "gage","globaltest")`  

- **num.workers** : number of thread for multithreading to decrease execution time. Default : `4`  
Example : `4`

- **shinyMode** : Boolean value. Shouldn't be changed. Only used when running coGSEA is running in the background of a shiny application. Default : `FALSE`

# Data preparation for running coGSEA

## MicroArray data

This example dataset is from a [article](https://www.ncbi.nlm.nih.gov/pubmed/21836821) about Astrocymotas tumors published in 2011 by Liu et al.  

You can download it on the Gene Expression Omnibus database with the accession number [GSE19728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19728)  


#### Loading the necessary packages for this analysis

``` r
library(affy)
library(hgu133plus2.db)
library(edgeR)
library(gage)
library(coGSEA)
```

#### Loading the CEL FILES and Normalizing them usin `rma`

``` r
setwd("~/GitLab/GSE19728/data/")
celfiles = ReadAffy()
celfiles = rma(celfiles)
```

#### Getting the expression data, and the Probe Names

``` r

intensity = exprs(celfiles)
intensity = cbind(rownames(intensity), intensity)
colnames(intensity)[1] = "PROBEID"
intensity = as.data.frame(intensity)
intensity$PROBEID = as.character(intensity$PROBEID)
```

#### Getting the annotations of Probes IDS to ENTREZ accession numbers

``` r
annots = select(hgu133plus2.db, intensity$PROBEID, "ENTREZID", "PROBEID")
```

#### Merge expression matrix and annotation

``` r
res = merge(intensity, annots, by= "PROBEID")

```


#### Getting rid of PROBE ids and type casting

``` r
resmin = res[,2:ncol(res)]
cname = colnames(resmin)
resmin = apply(resmin, 2, as.numeric)
colnames(resmin)= cname
resmin = as.data.frame(resmin)
```

#### Aggregating the PROBES matching the same ENTREZ accession number by averaging them and applying a log transformation

``` r
result = aggregate(. ~ ENTREZID, resmin, mean)
result$ENTREZID = levels(as.factor(as.character(resmin$ENTREZID)))
rownames(result) = result$ENTREZID
result = result[,-1]
result = log(result)
```

#### Visualizing in boxplot to check normalization

``` r
boxplot(result, las = 2)
```

#### Changing column names to remove the `.CEL` filename extension

``` r
colnames(result) = gsub(".CEL","", colnames(result))
print(colnames(result))
```


#### Selecting only the samples we are interested in

``` r
result2 = cbind(result$GSM492649_Astrocytomas_N, result$GSM525014, result$GSM525015, result$GSM525016, result$`GSM492662_Astrocytomas_T4-1`, result$`GSM492663_Astrocytomas_T4-2` , result$`GSM492664_Astrocytomas_T4-3`, result$`GSM492665_Astrocytomas_T4-4`, result$`GSM492666_Astrocytomas_T4-5`)
colnames(result2) = c("GSM492649_Astrocytomas_N", "GSM525014", "GSM525015", "GSM525016","GSM492662_Astrocytomas_T4-1", "GSM492663_Astrocytomas_T4-2", "GSM492664_Astrocytomas_T4-3", "GSM492665_Astrocytomas_T4-4", "GSM492666_Astrocytomas_T4-5" )
rownames(result2) = rownames(result)
```

#### Preparing the design matrix

``` r
Normal = c(rep(1,4),rep(0,5))
Tumor = c(rep(0,4),rep(1,5))
design = cbind(Normal, Tumor)
rownames(design) = colnames(result2)
```

#### Preparing the contrast matrix

``` r
contr.matrix = makeContrasts(NormalVSTumor = Normal - Tumor, levels = design)
```


#### Preparing expression list object

``` r
temp = new("EList")
temp$design = design
temp$E = as.matrix(result2)
rownames(temp$E) = as.numeric(rownames(temp$E))
temp$genes$ENTREZ = rownames(result2)
temp$common.dispersion = estimateDisp(temp$E, design = temp$design)$common.dispersion
temp$samples = colnames(result2)
```

#### Preparing gene set collection

``` r
gs = gage::kegg.gsets(species = "hsa", id.type = "entrez")
geneset = gs$kg.sets
```

Some GSEA methods do not work properly with exotic gene names, so we need to simplify them  

#### Function for simplifying gene sets names

``` r
nameshorter = function(names){
  namemod = c()
  for (i in seq(1,length(names))){
    namemod[i] = paste(strsplit(names[i], " ")[[1]][-1], sep = "", collapse = " ")
    namemod[i] = gsub("/","", names[i])
    namemod[i] = gsub(" ","_", names[i])
  }
  return(namemod)
}
```

#### Simplifying gene sets names

``` r
names(geneset) = nameshorter(names(geneset))
names(geneset) = gsub("/","_",names(geneset))
```

#### Saving necessary objects to RDS files

``` r
saveRDS(contr.matrix, "~/GitLab/GSE19728/contrast.rds")
saveRDS(temp, "~/GitLab/GSE19728/elist.rds")
saveRDS(geneset, "~/GitLab/GSE19728/geneset.rds")
```


#### Reading necessary objects (generated above) from RDS files

``` r
elist = readRDS("~/GitLab/GSE19728/elist.rds")
contrast = readRDS("~/GitLab/GSE19728/contrast.rds")
geneset = readRDS("~/GitLab/GSE19728/geneset.rds")
```



#### Running coGSEA analysis

``` r
coGSEA(ElistObject = elist, contrastMatrix = contrast, ENTREZGenesIds = elist$genes$ENTREZ, geneSetCollection = geneset,specie = "Homo sapiens", directoryPath = "~/GitLab/GSE19728/results", alpha = 0.05, pvalAdjMethod = "BH", pvalCombMethod = "sumlog",min.intersection.size = 1, GSEA.Methods = c("camera", "gage","globaltest", "gsva", "ssgsea", "zscore", "ora", "padog", "roast","safe"), num.workers = 4, shinyMode = FALSE)

```

## RNAseq data

You can download this example dataset on the Gene Expression Omnibus database with the accession number [GSE63310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310).    
This dataset was analyzed in a very detailed [article](https://f1000research.com/articles/5-1408/v2) on how to do differential expression analysis that we strongly advise you to read.  
The file you're looking for is : `GSE63310_RAW.tar`  

#### Loading necessary packages

``` r
library(edgeR)
library(limma)
library(Mus.musculus)
library(coGSEA)
```

#### Reading files

``` r
setwd("~/GitLab/GSE63310/data/")
files <- c(
"GSM1545535_10_6_5_11.txt",
"GSM1545536_9_6_5_11.txt",   
"GSM1545538_purep53.txt",
"GSM1545539_JMS8-2.txt",
"GSM1545540_JMS8-3.txt",
"GSM1545541_JMS8-4.txt",
"GSM1545542_JMS8-5.txt",
"GSM1545544_JMS9-P7c.txt",
"GSM1545545_JMS9-P8c.txt")
x <- readDGE(files, columns=c(1,3))
```

#### Simplifying file names

``` r
colnames(x) = substring(colnames(x), 12, nchar(colnames(x)))
```

#### Grouping by sample condition
``` r
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
x$samples$group <- group
```

#### Grouping by lane
``` r
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples
```

#### Annotation
``` r
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),keytype="ENTREZID")
dim(genes)
head(genes)
```

#### Getting rid of duplicated annoations by keeping only the first one
``` r
genes <- genes[!duplicated(genes$ENTREZID),]
```

#### Count per million of reads
``` r
cpm <- cpm(x)
```


#### Removing genes lowly expressed  (genes with 0 expression across all samples)
``` r
#Removing genes lowly expressed (genes with 0 expression across all samples)
table(rowSums(x$counts==0)==9)

keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

```

#### Normlization with edgeR
``` r
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
```

#### Making the design matrix
``` r
design <- model.matrix( ~ 0 + group + lane)
colnames(design) <- gsub("group", "", colnames(design))
rownames(design) = colnames(x)
```

#### Making the contrast matrix
``` r

contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design))

contr.matrix

```

#### Applying voom transformation
``` r
v <- voom(x, design, plot=F)
v$genes$ENTREZ = rownames(v$E)
```


#### Saving objects
``` r
saveRDS(v, "~/GitLab/GSE63310/elist.rds")
saveRDS(contr.matrix, "~/GitLab/GSE63310/contrast.rds")
```


#### Reading RDS objects (previously genereated above)
``` r
elist = readRDS("~/GitLab/GSE63310/elist.rds")
contrast = readRDS("~/GitLab/GSE63310/contrast.rds")
```


#### Running coGSEA
``` r
coGSEA(ElistObject = elist, contrastMatrix = contrast, ENTREZGenesIds = elist$genes$ENTREZ, geneSetCollection = "C2_KEGG",specie = "Mus musculus", directoryPath = "~/GitLab/GSE63310/results", alpha = 0.05, pvalAdjMethod = "BH", pvalCombMethod = "sumlog",min.intersection.size = 1, GSEA.Methods = c("camera", "gage","globaltest", "gsva", "ssgsea", "zscore", "ora", "padog", "roast","safe"), num.workers = 4, shinyMode = FALSE)

```
