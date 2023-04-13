# 🐶 hlabud: tidy HLA amino acid analysis

hlabud is an R package that provides functions to download and analyze
human leukocyte antigen (HLA) genotypes from
[IMGTHLA](https://github.com/ANHIG/IMGTHLA) in a tidy R workflow.

For example, let’s consider a simple question about two HLA genotypes
DRB1\*04:174 and DRB1\*15:152. What amino acid positions are different
between these two genotypes?

```r
library(hlabud)
a <- hla_alignments("DRB1")
amino_dosage(c("DRB1*04:174", "DRB1*15:152"), a$onehot)
```

```
            P9_E P9_W
DRB1*04:174    1    0
DRB1*15:152    0    1
```

From this output, we can see that the two genotypes are nearly
identical, but they differ at position 9:

-   DRB1\*04:174 has P9\_E (Glu)
-   DRB1\*15:152 has P9\_W (Trp)

Installation
============

The quickest way to get hlabud is to install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("slowkow/hlabud")
```

Examples
========

See the [usage examples](vignettes/examples.md) to get some ideas for
how to use hlabud in your analyses.

-   [Get a one-hot encoded matrix for all HLA-DRB1
    alleles](articles/examples.html#get-a-one-hot-encoded-matrix-for-all-hla-drb1-alleles)
-   [Convert genotypes to a dosage
    matrix](articles/examples.html#convert-genotypes-to-a-dosage-matrix)
-   [Logistic regression association for amino acid
    positions](articles/examples.html#logistic-regression-association-for-amino-acid-positions)
-   [UMAP embedding of 3,486 HLA-DRB1
    alleles](articles/examples.html#umap-embedding-of-3486-hla-drb1-alleles)
-   [Download and unpack all data from the latest IMGTHLA
    release](articles/examples.html#download-and-unpack-all-data-from-the-latest-imgthla-release)

<a href="articles/examples.html">
<img width="49%" src="articles/examples_files/figure-html/glm-volcano-1.png">
<img width="49%" src="articles/examples_files/figure-html/umap1-1.png">
</a>

Related work
============

I recommend this article for anyone new to HLA, because the beautiful
figures help to build intuition:

-   La Gruta NL, Gras S, Daley SR, Thomas PG, Rossjohn J. [Understanding
    the drivers of MHC restriction of T cell
    receptors.](https://www.ncbi.nlm.nih.gov/pubmed/29636542) Nat Rev
    Immunol. 2018;18: 467–478.

Learn about the conventions for HLA nomenclature:

-   Marsh SGE, Albert ED, Bodmer WF, Bontrop RE, Dupont B, Erlich HA, et
    al. [Nomenclature for factors of the HLA
    system, 2010.](https://www.ncbi.nlm.nih.gov/pubmed/20356336) Tissue
    Antigens. 2010;75: 291–455.

For case-control analysis of HLA genotype data, consider the
[BIGDAWG](https://CRAN.R-project.org/package=BIGDAWG) R package
available on CRAN. Here is the related article:

-   Pappas DJ, Marin W, Hollenbach JA, Mack SJ. [Bridging ImmunoGenomic
    Data Analysis Workflow Gaps (BIGDAWG): An integrated case-control
    analysis pipeline.](https://pubmed.ncbi.nlm.nih.gov/26708359) Hum
    Immunol. 2016;77: 283–287.
