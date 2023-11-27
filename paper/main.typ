//#import "template.typ": *

#import "lapreprint.typ": template

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: template.with(
  title: "hlabud: HLA genotype analysis in R",
  short-title: "hlabud",
  venue: [bio#text(fill: red.darken(20%))[R]$chi$iv],
  // You can make all dates optional, however, `date` is by default `datetime.today()`
  //date: (
  //  (title: "Published", date: datetime(year: 2023, month: 08, day: 21)),
  //  (title: "Accepted", date: datetime(year: 2022, month: 12, day: 10)),
  //  (title: "Submitted", date: datetime(year: 2022, month: 12, day: 10)),
  //),
  logo: "logo.png",
  theme: blue.darken(20%),
  font-face: "Noto Sans",
  paper-size: "us-letter",
  authors: (
    (
      name: "Kamil Slowikowski",
      email: "kslowikowski@mgh.harvard.edu",
      affiliations: "1,2,3,4",
      orcid: "0000-0002-2843-6370"
    ),
    (
      name: "Alexandra-Chloe Villani",
      email: "avillani@mgh.harvard.edu",
      affiliations: "1,2,3,4",
      orcid: "0000-0001-7461-0408"
    ),
  ),
  affiliations: (
    (id: "1", name: "Center for Immunology and Inflammatory Diseases, Division of Rheumatology, Allergy an Immunology, Department of Medicine, Massachusetts General Hospital"),
    (id: "2", name: "Cancer Center, Massachusetts General Hospital"),
    (id: "3", name: "Broad Institute"),
    (id: "4", name: "Harvard Medical School")
  ),
  kind: "Pre-Print",
  // Insert your abstract after the colon, wrapped in brackets.
  abstract: (
    (title: "Summary", content: [The human leukocyte antigen (HLA) genes have thousands of different alleles in the human population, and have more associations with human diseases than any other genes. Data for all known HLA genotypes are curated in the international ImMunoGeneTics (IMGT) database in versioned releases on #link("https://github.com/ANHIG/IMGTHLA")[GitHub]. Here, we introduce _hlabud_, an R package that provides access to data from the IMGT/HLA database and the Allele Frequency Net Database (AFND), functions to encode the data in different formats, and tutorials for association analysis, embedding, and HLA divergence.]),
    (title: "Availability", content: [Source code and documentation are available at *#link("https://github.com/slowkow/hlabud")[github.com/slowkow/hlabud]*]),
    (title: "Contact", content: [#link("mailto:kslowikowski@mgh.harvard.edu")[kslowikowski\@mgh.harvard.edu]])
  ),
  keywords: ("immunoinformatics", "genetics", "immunology", "HLA"),
  open-access: true,
)

= Introduction

Human leukocyte antigen (HLA) genes encode the proteins that display antigens so the immune system can recognize pathogens such as bacteria and viruses.
Geneticists have identified thousands of variants (e.g. single nucleotide polymorphisms) in the human genome that are associated with hundreds of different disease and phenotypes @Kennedy2017.

The HLA genes encode a protein complex that presents antigens to other cells.

To facilitate HLA genotype analysis, we developed _hlabud_, a free and open-source software package that downloads information from the IMGT/HLA database of HLA genotypes and sequence alignments @Robinson2020 directly in the R programming environment.
The _hlabud_ package provides functions that return convenient lists of items, where each item is either a matrix or a data frame.
The simple design makes _hlabud_ easy to integrate with any downstream R packages for data analysis or visualization.

_hlabud_ downloads HLA genotype data from the IMGT-HLA GitHub repository @imgthla and automatically caches it in a user-configurable folder.
Functionality includes parsing the custom IMGT/HLA file format for multiple sequence alignments, converting sequence alignments to a one-hot matrix, and calculating the Grantham divergence between HLA alleles @Pierini2018.

The documentation includes tutorials for analysis of the one-hot encoding of amino acid positions, including association analysis with logistic regression and low-dimensional embedding with UMAP @McInnes2018.
_hlabud_ also provides direct access to the allele frequencies for all HLA genes from the Allele Frequency Net Database (AFND) @Gonzalez-Galarza2020.

= Description

Comprehensive HLA genotype data is curated in the IMGT/HLA database, and the data is archived in a GitHub repository (#link("https://github.com/ANHIG/IMGTHLA")[github.com/ANHIG/IMGTHLA]).
We can use _hlabud_ to download the sequence alignment data, read it into R, and automatically encode the data as a one-hot matrix like this:

```R
a <- hla_alignments("DRB1")
```

When the user runs this line of code, _hlabud_ will:

- Download data from the IMGT/HLA Github repository.

- Cache data files in a local folder that supports multiple releases of the data.

- Read the data into data frames and matrices for downstream analysis.

- Create a one-hot encoding of the multiple sequence alignment data.

Many amino acid residues at specific loci have been associated with human diseases and blood protein levels @Krishna2023.
Researchers have developed software tools for calling HLA genotypes with high accuracy from DNA-seq or RNA-seq next-generation sequencing reads @Claeys2023, so there are opportunities to use that data for association studies.

Once we have a list of genotypes for each individual (e.g. `"DRB1*04:01,DRB1*05:01"`), we can use _hlabud_ to prepare data for regression analysis to find which amino acid positions are associated with a phenotype in a sample of individuals. We call `dosage(genotypes, a$onehot)` where `genotypes` is a vector of genotypes and `a$onehot` is a one-hot matrix representation of HLA alleles (from the example above). The `dosage()` function returns the number of copies of each amino acid at each position for each individual, which can then be used for omnibus regression @Sakaue2023 or single-position testing (@fig1\A).

UMAP accepts the one-hot matrix of HLA alleles as input, and it can be used to visualize the dataset in a latent space with reduced dimensionality (@fig1\B).

_hlabud_ provides direct access to the allele frequencies HLA genes reported in the Allele Frequency Net Database (AFND) (#link("http://allelefrequences.net")) (@fig1\C).

Each HLA allele binds a specific set of peptides.
So, an individual with two highly dissimilar alleles can bind a greater number of different peptides than a homozygous individual @Wakeland1990.
_hlabud_ implements the Grantham divergence calculations based on the original Perl code @Pierini2018:

```R
my_genos <- c("A*23:01:12,A*24:550", "A*25:12N,A*11:27", "A*24:381,A*33:85")
hla_divergence(my_genos, method = "grantham")
#> A*23:01:12,A*24:550    A*25:12N,A*11:27    A*24:381,A*33:85 
#>           0.4924242           3.3333333           4.9015152
```

#figure(
  image("figure-examples.png", width: 71%),
  caption: [(*A*) Association between amino acid positions and simulated case-control status. The x-axis represents the odds ratio and the y-axis represents $-log_10 P$ from a logistic regression analysis in R.
  (*B*) 3,516 HLA-DRB1 alleles represented as dots in a two-dimensional embedding computed by UMAP from a one-hot encoding of amino acids.
  (*C*) Allele frequencies for HLA-DQB1*02:01 in the AFND.],
) <fig1>


= Installation and documentation

_hlabud_ can be installed in an R session with:

```R
remotes::install_github("slowkow/hlabud")
```

Each function is documented extensively, and the complete manual can be viewed on the _hlabud_ website at #link("https://slowkow.github.io/hlabud"). _hlabud_ has been tested on Linux/Unix, Mac OS (Darwin) and Windows.

= Discussion

Our open-source R package _hlabud_ enables easy access to HLA data from two public databases, and provides functions to enable HLA divergence calculations, regression analysis, and low-dimensional embedding. We hope that _hlabud_ will raise awareness of the IMGT/HLA and AFND databases and influence other developers to share more open-source tools for HLA analysis. We envision that _hlabud_ will be used by biomedical researchers, and also by teachers and students who study genetics and bioinformatics.

= Acknowledgments

This work was supported by a NIAID grant T32AR007258 (to K.S.) and the National Institute of Health Director’s New Innovator Award (DP2CA247831; to A.C.V.) Thanks to Sreekar Mantena for reporting issues with the code. 

= Competing Interests

No competing interest is declared.

= Author contributions statement

K.S. wrote the software and the manuscript. A.C.V. reviewed the manuscript.

= Related Work

BIGDAWG is an R package that provides functions for chi-squared Hardy-Weinberg and case-control association tests of highly polymorphic genetic data like HLA genotypes @Pappas2016. HATK is set of Python scripts for processing and analyzing IMGT-HLA data @Choi2020.


#bibliography("references.bib")
