# coGSEA
**co**mparative **G**ene **S**et **E**nrichment **Analysis**

## Information

This version doesn't include the *SetRank* method for reason of package installation speeed optimization. To include it, you can use the 	[coGSEA_SetRank](https://gitlab.pasteur.fr/mborry/coGSEA_SetRank) version of this package.

## How to install

```
require(devtools)
devtools::install_git('https://mborry@gitlab.pasteur.fr/mborry/coGSEA.git')
```

## Input files

- An [elist](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/EList.html) object. This object is very similar to a R list. [Example data file](./exampleData/elist.rds)
- A contrast matrix. [Example data](./exampleData/contrast.rds)
- Optionally, a custome geneset collection can be added as well if the collection is not H, C2_KEGG, C2_REACTOME. This object is of list class. [Example data](./exampleData/geneset.rds)

More information about the method can be found in the [method introduction](./intro.md)




#### Disclaimer :
This tool is largely inspired by the [eGSEA](http://bioconductor.org/packages/release/bioc/html/EGSEA.html) R package (and contains some of its code)
