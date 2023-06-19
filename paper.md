---
title: 'hlabud: An R package for analysis of HLA genotypes'
tags:
  - R
  - genetics
  - immunology
  - bioinformatics
authors:
  - name: Kamil Slowikowski
	orcid: 0000-0002-2843-6370
	corresponding: true
    affiliation: "1, 2, 3, 4"
  - name: Alexandra-Chloé Villani
	orcid: 0000-0001-7461-0408
    affiliation: "1, 2, 3, 4, 5"
affiliations:
  - name: Center for Immunology and Inflammatory Diseases, Department of Medicine, Massachusetts General Hospital, Boston, MA, USA 
    index: 1
  - name: Massachusetts General Hospital, Cancer Center, Boston, MA, USA 
    index: 2
  - name: Broad Institute of Massachusetts Institute of Technology and Harvard, Cambridge, MA, USA
    index: 3
  - name: Harvard Medical School, Boston, MA, USA
	index: 4
  - name: Division of Rheumatology, North Shore Physicians Group, Department of Medicine, Mass General Brigham Healthcare Center, Lynn, MA, USA
	index: 5
date: 12 June 2023
bibliography: paper.bib

---

# Summary

Human leukocyte antigen (HLA) genes encode the proteins that display antigens for the immune system to recognize pathogens like bacteria and viruses.
Genes in the HLA locus on chromosome 6 in the human genome have thousands of different alleles in the human population.
The single-nucleotide polymorphisms (SNPs) encoding different amino acids in HLA genes are the genetic variants with the largest effect sizes, and they are associated with risk of developing autoimmune disease [@Kennedy2017].
The fields of immunology and genomics aim to discover molecular factors, such as HLA genotypes, that explain the functions of the human immune system in health and disease.
Analysis of this genotype data requires computational methods for managing collections of genetic data and transforming data into different encodings for downstream analyses [@Sakaue2022].

-![HLA-DRB1 genotypes embedded with UMAP](vignettes/examples_files/figure-html/umap1-1.png)


# Statement of need

`hlabud` is an R package that simplifies the tasks of downloading and parsing data from the IMGT/HLA database of HLA genotypes and sequence alignments [@Robinson2020].
The R programming language has a comprehensive repository of open-source libraries for statistical modeling and data visualization that can be applied to any data analysis.
The API for `hlabud` is designed to provide functions that output simple lists of matrices and tables that facilitate seamless integration with all other R packages.
HLA genotype data is lazily downloaded (as-needed) from the IMGT-HLA GitHub repository [@imgthla] and automatically cached in a user-configurable directory.
The documentation includes usage examples for analysis of the one-hot encoding of amino acid positions such as association analysis with logistic regression and low dimensional embedding with UMAP.

The BIGDAWG R package provides functions for chi-squared Hardy-Weinberg and case-control association tests of highly polymorphic genetic data like HLA genotypes [@Pappas2016].
In contrast, `hlabud` was designed to be minimal and flexible, to support any downstream HLA analyses, including analyses that have yet to be imagined.
`hlabud` can be used by biomedical researchers and also by students in courses that teach immunology, genetics, and bioinformatics.


# References


