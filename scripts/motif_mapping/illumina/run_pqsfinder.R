#!/usr/bin/env Rscript
# pqsfinder -- exhaustive, imperfection-tolerant G-quadruplex prediction.
# Confirmed API directly from the installed package: pqsfinder(subject, ...)
# takes a single DNAString, returns a PQSViews object whose coordinates come
# from as.data.frame() and whose score/strand/run-length metadata comes from
# S4Vectors::mcols() -- verified interactively before writing this, not
# guessed from documentation alone.
# Usage: Rscript run_pqsfinder.R <SEQUENCES_FASTA> <OUT_GFF3>
suppressMessages({
  library(pqsfinder)
  library(Biostrings)
  library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
seqs_file <- args[1]
out_file <- args[2]

seqs <- readDNAStringSet(seqs_file)

all_gff <- character(0)
for (i in seq_along(seqs)) {
  pqs <- pqsfinder(seqs[[i]], verbose = FALSE)
  if (length(pqs) == 0) next
  coords <- as.data.frame(pqs)
  meta <- as.data.frame(mcols(pqs))
  for (j in seq_len(nrow(coords))) {
    line <- paste(
      names(seqs)[i], "pqsfinder", "G_quadruplex",
      coords$start[j], coords$end[j], meta$score[j], meta$strand[j], ".",
      sprintf("ID=G4_%d_%d", i, j),
      sep = "\t"
    )
    all_gff <- c(all_gff, line)
  }
}
writeLines(all_gff, out_file)
cat(sprintf("Wrote %d G-quadruplex predictions to %s\n", length(all_gff), out_file))
