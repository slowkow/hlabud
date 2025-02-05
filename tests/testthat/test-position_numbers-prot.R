
library(stringr)
library(glue)

my_release <- "3.59.0"
my_type <- "prot"
my_genes <- unique(sort(hla_genes(release = my_release)$gene))

test_that(glue('positions are correct for HLA-A {my_type} {my_release}'), {

  a <- hla_alignments(gene = "A", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNM", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 24:1), as.character(1:67))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})

test_that(glue('positions are correct for HLA-B {my_type} {my_release}'), {

  a <- hla_alignments(gene = "B", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEY", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 24:1), as.character(1:59))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})

test_that(glue('positions are correct for HLA-C {my_type} {my_release}'), {

  a <- hla_alignments(gene = "C", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MRVMAPRTLILLLSGALALTETWACSHSMKYFFTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQ", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 24:1), as.character(1:65))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})

test_that(glue('positions are correct for HLA-DMA {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DMA", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MGHEQNQGAALLQMLPLLWLLPHSWAVPEAPTPMWPDDLQNHTFLHTVYCQDGSPSVGLSEAYDEDQLFFFDFSQNTRVPRLPEFADWAQEQGDAP", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 26:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DMB {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DMB", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MITFLPLLLGLSLGCTGAGGFVAHVESTCLLDDAGTPKDFTYCISFNKDLLTCWDPEENKMAPCEFGVLNSLANVLSQHLNQKDTLMQRLRNGLQNCA", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 18:1), as.character(1:80))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DOA {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DOA", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAH", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 25:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DOB {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DOB", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MGSGWVPWVVALLVNLTRLDSSMTQGTDSPEDFVIQAKADCYFTNGTEKVQFVVRFIFNLEEYVRFDSDVGMFVALTKLGQPDAEQWNSRLDLLER", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 26:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DPA1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DPA1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MRPEDRMFHIRAVILRALSLAFLLSLRGAGAIKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGL", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 31:1), as.character(1:60))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DPB1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DPB1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MMVLQVSAAPRTVALTALLMVLLTSVVQGRATPENYVYQGRQECYAFNGTQRFLERYIYNREEYARFDSDVGEFRAVTELGRPAAEYWNSQK", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 29:1), as.character(1:63))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DQA1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DQA1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVA", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 23:1), as.character(1:69))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})

test_that(glue('positions are correct for HLA-DQA2 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DQA2", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MILNKALLLGALALTAVMSPCGGEDIVADHVASYGVNFYQSHGPSGQYTHEFDGDEEFYVDLETKETVWQLPMFSKFISFDPQSALRNMAVGK", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 23:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DQB1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DQB1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MSWKKSLRIPGDLRVATVTLMLAILSSSLAEGRDSPEDFVYQFKGLCYFTNGTERVRGVTRHIYNREEYVRFDSDVGVYRAVTPQGRPVAEY", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 32:1), as.character(1:60))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DQB2 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DQB2", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MSWKMALQIPGGFWAAAVTVMLVMLSTPVAEARDFPKDFLVQFKGMCYFTNGTERVRGVARYIYNREEYGRFDSDVGEFQAVTELGRSIEDW", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 32:1), as.character(1:60))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DRA {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DRA", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MAISGVPVLGFFIIAVLMSAQESWAIKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANL", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 25:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DRB1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DRB1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVCLKLPGGSCMTALTVTLMVLSSPLALAGDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLL", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 29:1), as.character(1:68))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DRB3 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DRB3", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVCLKLPGGSSLAALTVTLMVLSSRLAFAGDTRPRFLELRKSECHFFNGTERVRYLDRYFHNQEEFLRFDSDVGEYRAVTELGRPVAESWNSQKDL", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 29:1), as.character(1:67))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DRB4 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DRB4", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVCLKLPGGSCMAALTVTLTVLSSPLALAGDTQPRFLEQAKCECHFLNGTERVWNLIRYIYNQEEYARYNSDLGEYQAVTELGRPDAEYWNSQKDLLER", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 29:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-DRB5 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "DRB5", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVCLKLPGGSYMAKLTVTLMVLSSPLALAGDTRPRFLQQDKYECHFFNGTERVRFLHRDIYNQEEDLRFDSDVGEYRAVTELGRPDAEYWNSQKDFLED", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 29:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-E {my_type} {my_release}'), {

  a <- hla_alignments(gene = "E", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVDGTLLLLLSEALALTQTWAGSHSLKYFHTSVSRPGRGEPRFISVGYVDDTQFVRFDNDAASPRMVPRAPWMEQEGSEYWDRETRSARDT", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 21:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-F {my_type} {my_release}'), {

  a <- hla_alignments(gene = "F", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MAPRSLLLLLSGALALTDTWAGSHSLRYFSTAVSRPGRGEPRYIAVEYVDDTQFLRFDSDAAIPRMEPREPWVEQEGPQYWEWTTGYAKAN", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 21:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-G {my_type} {my_release}'), {

  a <- hla_alignments(gene = "G", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MVVMAPRTLFLLLSGALTLTETWAGSHSMRYFSAAVSRPGRGEPRFIAMGYVDDTQFVRFDSDSACPRMEPRAPWVEQEGPEYWEEETRNTKAH", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 24:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-HFE {my_type} {my_release}'), {

  a <- hla_alignments(gene = "HFE", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MGPRARPALLLLMLLQTAVLQGRLLRSHSLHYLFMGASEQDLGLSLFEALGYVDDQLFVFYDHESRRVEPRTPWVSSRISSQMWLQLSQSL", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 21:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-MICA {my_type} {my_release}'), {

  a <- hla_alignments(gene = "MICA", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MGLGPVFLLLAGIFPFAPPGAAAEPHSLRYNLTVLSWDGSVQSGFLTEVHLDGQPFLRCDRQKCRAKPQGQWAEDVLGNKTWDRETRDLTGNG", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 23:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-MICB {my_type} {my_release}'), {

  a <- hla_alignments(gene = "MICB", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MGLGRVLLFLAVAFPFAPPAAAAEPHSLRYNLMVLSQDESVQSGFLAEGHLDGQPFLRYDRQKRRAKPQGQWAEDVLGAKTWDTETEDLTENG", "")[[1]]
  names(aminos) <- c(sprintf("n%s", 23:1), as.character(1:70))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-TAP1 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "TAP1", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MASSRCPAPRGCRCLPGASLAWLGTVLLLLADWVLLRTALPRIFSLLVPTALPLLRVWAVGLSRWAVLWLGACGVLRATVGSKSENAGAQGWLAALKPLA", "")[[1]]
  names(aminos) <- as.character(1:length(aminos))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


test_that(glue('positions are correct for HLA-TAP2 {my_type} {my_release}'), {

  a <- hla_alignments(gene = "TAP2", type = my_type, release = my_release, verbose = FALSE)
  expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))

  # Manually copied from the 3.59.0 release file
  aminos <- strsplit("MRLPDLRPWTSLLLVDAALLWLLQGPLGTLLPQGLPGLWLEGTLRLGGLWGLLKLRGLLGFVGTLLLPLCLATPLTVSLRALVAGASRAPPARVASAPWS", "")[[1]]
  names(aminos) <- as.character(1:length(aminos))
  expect_equal(aminos, a$alleles[1,names(aminos)])

})


