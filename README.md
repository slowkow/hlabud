-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Examples](#examples)
-   [Related work](#related-work)

Introduction
============

hlabud is an R package that provides functions to download and analyze
human leukocyte antigen (HLA) genotypes from
[IMGTHLA](https://github.com/ANHIG/IMGTHLA) in a tidy R workflow.

What amino acid positions are different between the DRB1\*04:174 and
DRB1\*15:152 genotypes?

``` r
library(hlabud)
a <- hla_alignments("DRB1")
amino_dosage(c("DRB1*04:174", "DRB1*15:152"), a$onehot)
```

    ##             P9_E P9_W
    ## DRB1*04:174    1    0
    ## DRB1*15:152    0    1

The two genotypes are nearly identical, but they differ at position 9:

-   E (Glu) in DRB1\*04:174
-   W (Trp) in DRB1\*15:152

Installation
============

The quickest way to get hlabud is to install from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("slowkow/hlabud")
```

Examples
========

See the [usage examples](vignettes/examples.Rmd) to see how hlabud can
be used in an analysis.

Related work
============

I recommend this article for anyone new to HLA, because the beautiful
figures help to build intuition:

-   La Gruta NL, Gras S, Daley SR, Thomas PG, Rossjohn J. [Understanding
    the drivers of MHC restriction of T cell
    receptors.](https://www.ncbi.nlm.nih.gov/pubmed/29636542) Nat Rev
    Immunol. 2018;18: 467–478.
    <a href="doi:10.1038/s41577-018-0007-5" class="uri">doi:10.1038/s41577-018-0007-5</a>

Learn about the conventions for HLA nomenclature:

-   Marsh SGE, Albert ED, Bodmer WF, Bontrop RE, Dupont B, Erlich HA, et
    al. [Nomenclature for factors of the HLA
    system, 2010.](https://www.ncbi.nlm.nih.gov/pubmed/20356336) Tissue
    Antigens. 2010;75: 291–455.
    <a href="doi:10.1111/j.1399-0039.2010.01466.x" class="uri">doi:10.1111/j.1399-0039.2010.01466.x</a>

For case-control analysis of HLA genotype data, please consider the
[BIGDAWG](https://CRAN.R-project.org/package=BIGDAWG) R package
available on CRAN, as cite the associated article:

-   Pappas DJ, Marin W, Hollenbach JA, Mack SJ. [Bridging ImmunoGenomic
    Data Analysis Workflow Gaps (BIGDAWG): An integrated case-control
    analysis pipeline.](https://pubmed.ncbi.nlm.nih.gov/26708359) Hum
    Immunol. 2016;77: 283–287.
    <a href="doi:10.1016/j.humimm.2015.12.006" class="uri">doi:10.1016/j.humimm.2015.12.006</a>
