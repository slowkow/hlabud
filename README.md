# hlabud <img width="25%" align="right" src="https://github.com/slowkow/hlabud/assets/209714/b39a3f04-c9a8-4867-a3e0-9434f0f9ef20"></img>

[![R-CMD-check](https://github.com/slowkow/hlabud/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/slowkow/hlabud/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11093557.svg)](https://doi.org/10.5281/zenodo.11093557)

hlabud provides methods to retrieve sequence alignment data from
[IMGTHLA](https://github.com/ANHIG/IMGTHLA) and convert the data into
convenient R matrices ready for downstream analysis. See the [usage
examples](https://slowkow.github.io/hlabud/articles/examples.html) to
learn how to use the data with logistic regression and dimensionality
reduction. We also share tips on how to [visualize the 3D molecular
structure](https://slowkow.github.io/hlabud/articles/visualize-hla-structure.html)
of HLA proteins and highlight specific amino acid residues.

For example, let’s consider a simple question about two HLA genotypes.

What amino acid positions are different between two genotypes?

``` r
library(hlabud)
a <- hla_alignments("DRB1")
a$release
```

    ## [1] "3.59.0"

``` r
dosage(a$onehot, c("DRB1*03:01:05", "DRB1*03:02:03"))
```

    ##               F26 Y26 D28 E28 F47 Y47 G86 V86
    ## DRB1*03:01:05   0   1   1   0   1   0   0   1
    ## DRB1*03:02:03   1   0   0   1   0   1   1   0

What nucleotides are different?

``` r
n <- hla_alignments("DRB1", type = "nuc")
n$release
```

    ## [1] "3.59.0"

``` r
dosage(n$onehot, c("DRB1*03:01:05", "DRB1*03:02:03"))
```

    ##               A164 T164 C171 G171 A227 T227 A240 G240 G344 T344 G345 T345 A357 G357
    ## DRB1*03:01:05    1    0    1    0    0    1    1    0    0    1    1    0    1    0
    ## DRB1*03:02:03    0    1    0    1    1    0    0    1    1    0    0    1    0    1

# Installation

The quickest way to get hlabud is to install from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("slowkow/hlabud")
```

# Examples

See the [usage
examples](https://slowkow.github.io/hlabud/articles/examples.html) to
get some ideas for how to use hlabud in your analyses.

-   [Get a one-hot encoded matrix for all HLA-DRB1
    alleles](https://slowkow.github.io/hlabud/articles/examples.html#get-a-one-hot-encoded-matrix-for-all-hla-drb1-alleles)

-   [Convert genotypes to a dosage
    matrix](https://slowkow.github.io/hlabud/articles/examples.html#convert-genotypes-to-a-dosage-matrix)

-   [Logistic regression association for amino acid
    positions](https://slowkow.github.io/hlabud/articles/examples.html#logistic-regression-association-for-amino-acid-positions)

-   [UMAP embedding of 3,516 HLA-DRB1
    alleles](https://slowkow.github.io/hlabud/articles/examples.html#umap-embedding-of-3516-hla-drb1-alleles)

-   [Get HLA allele frequencies from Allele Frequency Net Database
    (AFND)](https://slowkow.github.io/hlabud/articles/examples.html#get-hla-allele-frequencies-from-allele-frequency-net-database-afnd)

-   [Compute HLA divergence with the Grantham distance
    matrix](https://slowkow.github.io/hlabud/articles/examples.html#compute-hla-divergence-with-the-grantham-distance-matrix)

-   [Download and unpack all data from the latest IMGTHLA
    release](https://slowkow.github.io/hlabud/articles/examples.html#download-and-unpack-all-data-from-the-latest-imgthla-release)

<a href="https://slowkow.github.io/hlabud/articles/examples.html#logistic-regression-association-for-amino-acid-positions">
<img width="49%" src="vignettes/articles/examples_files/figure-html/glm-volcano-1.png">
</a>
<a href="https://slowkow.github.io/hlabud/articles/examples.html#umap-embedding-of-3516-hla-drb1-alleles">
<img width="49%" src="vignettes/articles/examples_files/figure-html/umap-2digit-1.png">
</a>
<a href="https://slowkow.github.io/hlabud/articles/examples.html#get-hla-allele-frequencies-from-allele-frequency-net-database-afnd">
<img width="49%" src="vignettes/articles/examples_files/figure-html/afnd_dqb1_02_01-1.png">
</a>
<a href="https://slowkow.github.io/hlabud/articles/visualize-hla-structure.html">
<img width="49%" src="https://github.com/slowkow/ggrepel/assets/209714/4843a850-a4fd-4832-9600-0b8e9c1bb904">
</a>

# Citation

`hlabud` provides access to the data in IMGT/HLA database. Therefore, if
you use `hlabud` then please cite the IMGT/HLA paper:

-   Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE.
    [IPD-IMGT/HLA Database.](https://pubmed.ncbi.nlm.nih.gov/31667505/)
    Nucleic Acids Res. 2020;48: D948–D955.
    <https://doi.org/10.1093/nar/gkz950>

`hlabud` also provides access to the data in Allele Frequency Net
Database (AFND). Therefore, if you use `hlabud::hla_frequencies()` then
please cite the AFND paper:

-   Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L,
    Ortega-Rivera ND, et al. [Allele frequency net database (AFND) 2020
    update: gold-standard data classification, open access genotype data
    and new query tools.](https://pubmed.ncbi.nlm.nih.gov/31722398)
    Nucleic Acids Res. 2020;48: D783–D788.
    <https://doi.org/10.1093/nar/gkz1029>

Additionally, you can also cite the `hlabud` package like this:

-   Slowikowski K. hlabud: HLA analysis in R. Zenodo.
    <https://doi.org/10.5281/zenodo.11093557>

# Learn more

I recommend this article for anyone new to HLA, because the beautiful
figures help to build intuition:

-   La Gruta NL, Gras S, Daley SR, Thomas PG, Rossjohn J. [Understanding
    the drivers of MHC restriction of T cell
    receptors.](https://pubmed.ncbi.nlm.nih.gov/29636542/) Nat Rev
    Immunol. 2018;18: 467–478.

Learn about the conventions for HLA nomenclature:

-   Marsh SGE, Albert ED, Bodmer WF, Bontrop RE, Dupont B, Erlich HA, et
    al. [Nomenclature for factors of the HLA system,
    2010.](https://pubmed.ncbi.nlm.nih.gov/20356336/) Tissue Antigens.
    2010;75: 291–455.

# Related work

[HLAtools](https://github.com/sjmack/HLAtools) is an R package that also
makes IPD-IMGT/HLA resources available for analysis, and also works with
BIGDAWG data formats.

For case-control analysis of HLA genotype data, consider the
[BIGDAWG](https://CRAN.R-project.org/package=BIGDAWG) R package
available on CRAN. Here is the related article:

-   Pappas DJ, Marin W, Hollenbach JA, Mack SJ. [Bridging ImmunoGenomic
    Data Analysis Workflow Gaps (BIGDAWG): An integrated case-control
    analysis pipeline.](https://pubmed.ncbi.nlm.nih.gov/26708359) Hum
    Immunol. 2016;77: 283–287.

[HATK](https://github.com/WansonChoi/HATK) is set of Python scripts for
processing and analyzing IMGT-HLA data. Here is the related article:

-   Choi W, Luo Y, Raychaudhuri S, Han B. [HATK: HLA analysis
    toolkit.](https://pubmed.ncbi.nlm.nih.gov/32735319) Bioinformatics.
    2021;37: 416–418. <doi:10.1093/bioinformatics/btaa684>

The HLA divergence code in hlabud is a translation of the [original Perl
code](https://sourceforge.net/projects/granthamdist) by [Tobias
Lenz](https://orcid.org/0000-0002-7203-0044), which was published in:

-   Pierini F, Lenz TL. [Divergent allele advantage at human MHC genes:
    signatures of past and ongoing
    selection.](https://pubmed.ncbi.nlm.nih.gov/29893875) Mol Biol
    Evol. 2018. <doi:10.1093/molbev/msy116>

[HLAdivR](https://github.com/rbentham/HLAdivR/) is another R package for
calculating HLA divergence.
