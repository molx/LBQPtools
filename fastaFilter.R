#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

library(Biostrings)

f <- readAAStringSet("2seqFasta.fasta")


f1 <- f[[1]]

query <- "cytochrome"

ff <- f[grepl(query, names(f))]
ff

writeXStringSet(x = ff, filepath = "output.fasta")


##### Jaques

library(stringi)

csv <- read.csv("Raiz grandiflora BANCO TRANSCRIPTOMA 410 ptns.csv", skip = 2, stringsAsFactors = FALSE)

query <- csv$Accession

fasta <- readAAStringSet("proteins grandiflora.fasta")

test <- sapply(names(fasta), function(i) any(sapply(query, function(j) grepl(j, i, fixed = TRUE))))

sum(test)

out <- fasta[test]

writeXStringSet(x = out, filepath = "JaquesOut.fasta")






