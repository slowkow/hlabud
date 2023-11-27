//#import "template.typ": *

#import "lapreprint.typ": template

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: template.with(
  title: [hlabud: HLA genotype analysis in R],
  short-title: [_hlabud_],
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
    (
      title: "Summary",
      content: [The human leukocyte antigen (HLA) genes have thousands of different alleles in the human population, and have more associations with human diseases than any other genes. Data for all known HLA genotypes are curated in the international ImMunoGeneTics (IMGT) database, and the Allele Frequency Net Database (AFND) provides allele frequencies for each HLA allele across human populations. Our open-source R package _hlabud_ facilitates access to HLA data from IMGT/HLA and AFND, and provides functions for HLA divergence calculations, fine-mapping analysis of amino acid (or nucleotide) positions, and low-dimensional embedding.]
    ),
    (title: "Availability", content: [Source code and documentation are available at *#link("https://github.com/slowkow/hlabud")[github.com/slowkow/hlabud]*]),
    (title: "Contact", content: [#link("mailto:kslowikowski@mgh.harvard.edu")[kslowikowski\@mgh.harvard.edu]])
  ),
  keywords: ("immunoinformatics", "genetics", "immunology", "HLA"),
  open-access: true,
)

= Introduction

Human leukocyte antigen (HLA) genes encode the proteins that enable cells to display antigens to other cells, so the immune system can recognize pathogens such as bacteria and viruses.
Geneticists have identified thousands of variants (e.g. single nucleotide polymorphisms) in the human genome that are associated with hundreds of different diseases and phenotypes @Kennedy2017.

HLA nomenclature consists of allele names like _HLA*01:01_ to indicate the genotype of each individual in a study.
Each allele name corresponds to multiple mutations at different positions throughout the gene's sequence, so it is difficult to estimate the similarity of two alleles solely from the allele names.
This ambiguity about specific amino acid positions means that allele names are not ideal for statistical analysis.

Researchers have developed software tools for calling HLA genotypes (@diagram) with high accuracy from DNA-seq or RNA-seq next-generation sequencing reads @Claeys2023, so there may be opportunities to use this type of data for HLA association studies.
Most software tools report allele names, not genotypes at specific nucleotide positions.
Commercial providers of HLA typing services also report genotypes with the traditional HLA allele names (i.e. _HLA*01:01_) instead of reporting alleles at specific nucleotide positions (@diagram).

#figure(
  image("diagram.png", width: 130%),
  caption: [_hlabud_ converts HLA genotypes to amino acid position matrices.]
) <diagram>

In contrast, fine-mapping analysis involves associating a phenotype with each amino acid position.
Many amino acid residues at specific loci have been associated with human diseases and blood protein levels @Krishna2023.
Published associations at specific amino acid positions have created opportunities for experimental validation that might advance our understanding of disease-associated mechanisms related to HLA proteins.

Fine-mapping can be more sensitive than allele-level analysis, and the results can be interpreted in the context of the protein structures that are affected by the associated amino acid positions.
For example, we might have different ideas about the function of a mutation in the peptide binding groove than a mutation in the interior region of the protein.

To facilitate HLA fine-mapping analysis, we developed _hlabud_, a free and open-source R package that downloads data from the IMGT/HLA database @Robinson2020 and automatically creates amino acid (or nucleotide) position matrices that are ready for analysis (@diagram).
_hlabud_ functions return simple lists, where each item in the list is a matrix or a data frame.
This simple design makes _hlabud_ easy to integrate with any downstream R packages for data analysis or visualization.



= Examples

Curated HLA genotype data is provided by the IMGT/HLA database at GitHub (#link("https://github.com/ANHIG/IMGTHLA")[github.com/ANHIG/IMGTHLA]).
In the example below, we use _hlabud_ to download the sequence alignment data for _HLA-DRB1_, read it into R, and encode it as a one-hot matrix:

```R
a <- hla_alignments("DRB1")
```

With one line of code, _hlabud_ will:

- Download data from the IMGT/HLA Github repository.

- Cache files in a local folder that supports multiple data releases.

- Read the data into matrices and dataframes for downstream analysis.

- Create a one-hot encoding of the multiple sequence alignment data.

Once we have obtained a list of genotypes for each individual (e.g. `"DRB1*04:01,DRB1*05:01"`), we can use _hlabud_ to prepare data for fine-mapping regression analysis that will reveal which amino acid positions are associated with a phenotype in a sample of individuals. To calculate the number of copies of each amino acid at each position for each individual, we can run:

```R
dosage(genotypes, a$onehot)
```

where `genotypes` is a vector of _HLA-DRB1_ genotypes and `a$onehot` is a one-hot matrix representation of _HLA-DRB1_ alleles.
The dosage matrix can then be used for omnibus regression @Sakaue2023 or fine-mapping (i.e. regression with each single position) (@figexamples\A).

Visualizing data in a two-dimensional embedding with algorithms like UMAP @McInnes2018 can help to build intuition about the relationship between all objects in a dataset.
UMAP accepts the one-hot matrix of HLA alleles as input, and the resulting embedding can be used to visualize the dataset for exploratory data analysis (@figexamples\B).

_hlabud_ provides direct access to the allele frequencies of HLA genes in the Allele Frequency Net Database (AFND) @Gonzalez-Galarza2020 (#link("http://allelefrequences.net")) (@figexamples\C).

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
) <figexamples>


= Installation and documentation

The easiest way to install _hlabud_ is to run this command in an R session:

```R
remotes::install_github("slowkow/hlabud")
```

The complete manual is available at #link("https://slowkow.github.io/hlabud"). _hlabud_ has been tested on Linux/Unix, Mac OS (Darwin) and Windows.

= Discussion

Our open-source R package _hlabud_ enables easy access to HLA data from two public databases and provides functions for HLA divergence calculations, amino acid or nucleotide fine-mapping analysis, and low-dimensional embedding. 
_hlabud_ downloads HLA genotype data from the IMGT-HLA GitHub repository @imgthla, caches it in a user-configurable folder, and prepares the data for downstream analysis in R.

We provide tutorials for HLA divergence calculation, fine-mapping association analysis with logistic regression, and embedding with UMAP.
_hlabud_ also provides direct access to the allele frequencies for all HLA genes from the Allele Frequency Net Database (AFND) @Gonzalez-Galarza2020.

= Acknowledgments

This work was supported by a NIAID grant T32AR007258 (to K.S.) and the National Institute of Health Director’s New Innovator Award (DP2CA247831; to A.C.V.) Thanks to Sreekar Mantena for reporting issues with the code. 

= Competing Interests

No competing interest is declared.

= Author contributions statement

K.S. wrote the software and the manuscript. A.C.V. reviewed the manuscript.

= Related Work

BIGDAWG is an R package that provides functions for chi-squared Hardy-Weinberg and case-control association tests of highly polymorphic genetic data like HLA genotypes @Pappas2016. HATK is set of Python scripts for processing and analyzing IMGT-HLA data @Choi2020.


#bibliography("references.bib")
