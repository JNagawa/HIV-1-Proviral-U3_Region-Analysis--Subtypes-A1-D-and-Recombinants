#!/usr/bin/env Rscript
# TFBSTools -- PWM motif scanning against the 6 core-TF JASPAR set, using
# TFBSTools::readJASPARMatrix on the raw JASPAR-format flat file (not the
# MEME-format one FIMO uses). Confirmed API from TFBSTools documentation
# (searchSeq/toPWM/PWMatrixList), not guessed.
# Usage: Rscript run_tfbstools.R <SEQUENCES_FASTA> <JASPAR_FLATFILE> <OUT_GFF3>
suppressMessages({
  library(TFBSTools)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
seqs_file <- args[1]
jaspar_file <- args[2]
out_file <- args[3]

pfms <- readJASPARMatrix(jaspar_file, matrixClass = "PFM")
pwms <- toPWM(pfms)

seqs <- readDNAStringSet(seqs_file)

all_gff <- character(0)
for (i in seq_along(seqs)) {
  hits <- searchSeq(pwms, seqs[[i]], seqname = names(seqs)[i], min.score = "80%", strand = "*")
  gff <- writeGFF3(hits)
  if (nrow(gff) > 0) {
    all_gff <- c(all_gff, capture.output(write.table(gff, sep = "\t", quote = FALSE,
                                                       row.names = FALSE, col.names = FALSE)))
  }
}
writeLines(all_gff, out_file)
cat(sprintf("Wrote %d hits to %s\n", length(all_gff), out_file))
