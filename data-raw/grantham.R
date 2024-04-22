
grantham <- read.table(header = TRUE, text = "amino c p v
Ser	1.42	9.2	32
Arg	0.65	10.5	124
Leu	0	4.9	111
Pro	0.39	8	32.5
Thr	0.71	8.6	61
Ala	0	8.1	31
Val	0	5.9	84
Gly	0.74	9	3
Ile	0	5.2	111
Phe	0	5.2	132
Tyr	0.2	6.2	136
Cys	2.75	5.5	55
His	0.58	10.4	96
Gln	0.89	10.5	85
Asn	1.33	11.6	56
Lys	0.33	11.3	119
Asp	1.38	13	54
Glu	0.92	12.3	83
Met	0	5.7	105
Trp	0.13	5.4	170
")

usethis::use_data(grantham, overwrite = TRUE)
