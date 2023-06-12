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
    affiliation: "1, 2, 3, 4, 5"
  - name: Alexandra-Chloe Villani
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
  - name: Harvard Medical School, Boston, MA, USA
	index: 5
date: 12 June 2023
bibliography: paper.bib

---

# Summary

Human leukocyte antigen (HLA) genes encode the proteins that display peptides from pathogens such as bacteria and viruses.
Genes in the HLA locus on chromosome 6 in the human genome have thousands of different alleles in the human population.
The single-nucleotide polymorphisms (SNPs) encoding different amino acids in HLA genes are the genetic variants with the largest effect sizes, and they are associated with risk of developing autoimmune disease.
The fields of immunology and genomics aim to discover molecular factors, such as HLA genotypes, that explain the functions of the human immune system in health and disease.
Analysis of this genotype data requires computational tools for managing collections of genetic data and transforming data into different encodings for downstream analyses.

-![HLA-DRB1 genotypes embedded with UMAP](vignettes/examples_files/figure-html/umap1-1.png)


# Statement of need

`hlabud` is an R package for conveniently downloading and analyzing the data from the IMGT/HLA database of HLA genotypes [@Robinson2020].
R provides a rich ecosystem of libraries for statistical modeling and data visualization that can be applied to any data analysis.
The API for `hlabud` is designed to provide simple functions that output lists of matrices and data frames to facilitate seamless integration with any other packages.
The documentation includes usage examples such as logistic regression association analysis with amino acid positions and UMAP [@McInnes2018] analysis of the one-hot encoded amino acid data.
HLA genotype data is downloaded lazily, as-needed, and cached in a user-configurable directory.

`hlabud` was designed to be used by biomedical researchers and by students in courses on genetics and bioinformatics.


# References


