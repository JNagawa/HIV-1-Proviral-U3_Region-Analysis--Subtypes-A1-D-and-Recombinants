#!/usr/bin/env Rscript
# gquad -- regex-based G-quadruplex prediction. Confirmed API directly from
# the installed package (?gquad, and an interactive test on a toy
# multi-sequence FASTA before writing this, not guessed from the CRAN PDF
# alone): gquad(x, xformat="fasta") takes a FASTA file PATH directly (not a
# loaded DNAStringSet) and returns a data frame with columns input_ID,
# sequence_position, sequence, sequence_length, likeliness ("*"/"**").
# Usage: Rscript run_gquad.R <SEQUENCES_FASTA> <OUT_GFF3>
suppressMessages(library(gquad))

args <- commandArgs(trailingOnly = TRUE)
seqs_file <- args[1]
out_file <- args[2]

# gquad's own input_ID is just a 1-based row index into the FASTA, not the
# record's actual header -- read headers separately so the GFF3 output
# names the real sequence, not "1"/"2".
headers <- sub("^>", "", grep("^>", readLines(seqs_file), value = TRUE))
headers <- sub("[ \t].*", "", headers)

res <- gquad(seqs_file, xformat = "fasta")

all_gff <- character(0)
if (is.data.frame(res) && nrow(res) > 0) {
  for (i in seq_len(nrow(res))) {
    seqname <- headers[as.integer(res$input_ID[i])]
    start <- as.integer(res$sequence_position[i])
    end <- start + as.integer(res$sequence_length[i]) - 1
    line <- paste(
      seqname, "gquad", "G_quadruplex",
      start, end, ".", "+", ".",
      sprintf("ID=G4_%d;likeliness=%s", i, res$likeliness[i]),
      sep = "\t"
    )
    all_gff <- c(all_gff, line)
  }
}
writeLines(all_gff, out_file)
cat(sprintf("Wrote %d G-quadruplex predictions to %s\n", length(all_gff), out_file))
