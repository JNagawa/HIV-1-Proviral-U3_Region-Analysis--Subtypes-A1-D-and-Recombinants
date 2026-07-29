#!/usr/bin/env bash

################################################################################

# Note that for boolean variables, only the exact value "true" (all lower case)
# will be interpreted as true, anything else is taken to mean false.

# What do you have to type into the command line to make these commands execute?
# (If the binary file lives in a directory that is not included in your $PATH
# variable, you will need to include the path here.)
python='python3'                                     # command used to invoke Python
BlastDBcommand='makeblastdb'                          # command to build a blast database
BlastNcommand='blastn'                                # command to run nucleotide blast
# smalt mapper command (used only if mapper=smalt)
smalt='smalt'
bwa='bwa'                                             # bwa mapper command (used only if mapper=bwa)
# bowtie2 mapper command (used only if mapper=bowtie)
bowtie2='bowtie2'
bowtie2_build='bowtie2-build'                         # command to build a bowtie2 index
samtools='samtools'                                   # samtools command for BAM/pileup handling
# mafft command for multiple-sequence alignment
mafft='mafft'
fastaq='fastaq'                                       # fastaq command for FASTA/FASTQ manipulation
# If shiver is installed with conda, you can run trimmomatic simply typing
# 'trimmomatic' at the command line. Otherwise, to be able to do that,
# 1. in a file named 'trimmomatic' (no file extension) copy the next three lines:
# ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# TheBinary=$(ls "$ThisDir"/trimmomatic-*.jar)
# java -jar "$TheBinary" "$@" || { echo "Problem running Trimmomatic."; exit 1; }
# and remove the # character at the start of each line
# 2. make that file executable, e.g. running the 'chmod u+x' command on it from
# the command line
# 3. move that file to the same directory that contains the trimmomatic java
# file you have downloaded, usually named like trimmomatic-XXX.jar (with numbers
# instead of XXX)
# 4. Add that directory to your PATH variable (Google how to do this if needed).
# After those four steps, you can leave the variable below set to 'trimmomatic'.
trimmomatic='trimmomatic'                             # command to run Trimmomatic for read trimming
# If you leave 'GiveHXB2coords', below, as 'true', we'll do pairwise alignment
# of the mapping reference with HXB2. You may as well use mafft options to make
# it more accurate (though slower).
# mafft options for the accurate pairwise reference-to-HXB2 alignment
MafftArgsForPairwise='--maxiterate 1000 --localpair'

# Minimum contig length: contigs shorter than this will be discarded at the
# start. In addition, when contigs are blasted against the existing reference
# set, we will only keep hits for which the length of the hit multipled by its
# identity to the reference is at least this length.
# drop contigs (and weak blast hits) shorter than 300bp
MinContigLength=300
# A contig will be split/cut if it has multiple blast hits. After alignment to
# a set of references, it may be split again if it contains a large gap: this
# parameter sets the gap size that will result in such a splitting...
# split an aligned contig at internal gaps this long or longer
MinGapSizeToSplitGontig=160
# ...and as such splitting occasionally results in cutting off a small bit of
# a contig into a new separate contig that you might not want to bother keeping,
# we have a second length threshold (which you could set to equal the
# MinContigLength parameter above but by default we are more permissive).
# keep post-split contig fragments at least this long
MinContigFragmentLength=80

# After aligning the contigs to the input existing references, by default we
# trim off any contig sequence that overhangs the whole reference alignment,
# i.e. in the alignment the beginning of the contig is to the left of the
# beginning of all references, or the end of the contig is to the right of the
# end of all references) Set the variable below to false to switch this
# trimming off.
# trim contig ends that overhang the reference alignment
TrimToKnownGenome=true

# Different blast 'tasks' (modes) to try when blasting the contigs. We will run
# each of the different tasks specified here (separated by whitespace) and merge
# the results. With this default, we try only -task megablast; if
# 'megablast blastn' were specified instead, we would also try -task blastn (and
# merge results).
# blast task(s) used when matching contigs to references
BlastTasks='megablast'

# Options to give blast when blasting the contigs; run your blastn command with
# -help to investigate possibilities.
# extra blastn options for the contig-vs-reference search
ContigBlastArgs="-max_target_seqs 1 -word_size 17"

# When two blast hits for the same contig have a fractional overlap (defined as
# the length of the part of the contig spanned by both hits divided by the
# length of the shorter of the two hits) equal to or greater than this value,
# we will merge them into a single hit. When the fractional overlap is less than
# this value, the two hits will be kept separate, resulting in the contig being
# split into two parts (one corresponding to each hit) to be aligned separately.
# A value of 1 or greater means partially overlapping hits are never merged
# (which is how shiver has always behaved). A value between 0 and 1 means they
# may or may not be merged, depending on how strongly they overlap.
# fractional overlap at/above which two blast hits are merged
ContigMinBlastOverlapToMerge='0.8'

# If you have a more recent mafft installation that includes the --addfragments
# option, we will use both --addfragments and --add to align the contigs to the
# existing reference alignment, then automatically choose one to keep. There are
# two available strategies for doing this: one is to use the alignment with the
# shortest length ("MinAlnLength"), the other is to calculate the fractional gap
# content of each contig after alignment, find the maximum over all contigs, and
# use the alignment with the smaller maximum ("MinMaxGappiness").
# rule for choosing between the --add and --addfragments alignments
MafftTestingStrategy="MinAlnLength"

# Shall we trim adapaters and low quality bases from reads, using trimmomatic?
# enable Trimmomatic adapter/quality trimming of reads
TrimReadsForAdaptersAndQual=true
# The trimmomatic manual explains at length the parameters controlling read
# trimming; the reader is referred to it for explanations of the following
# variables and other options not used here:
# Trimmomatic ILLUMINACLIP adapter-clipping parameters
IlluminaClipParams='2:10:7:1:true'
# Trimmomatic quality/length trimming parameters
BaseQualityParams='MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20'
# How many threads Trimmomatic should use (it sometimes multithreads unless told
# not to, which can be problematic on clusters).
# pin Trimmomatic to 1 thread to stay well-behaved on clusters
NumThreadsTrimmomatic=1

# Shall we trim exact matches to PCR primers from the end of reads using fastaq?
# enable trimming of exact PCR-primer matches from read ends
TrimReadsForPrimers=true
# Shall we also trim matches to the PCR primers that differ by a single base
# change? (This slows down the trimming step a lot.)
# do not also trim near-primer matches (1 SNP off) -- too slow
TrimPrimerWithOneSNP=false

# Shall we clean (remove read pairs that look like contaminants)?
# remove read pairs flagged as likely contaminants
CleanReads=true

# Which mapper to use? "smalt", "bwa" or "bowtie"? You can ignore the options
# for a mapper you're not using, and it doesn't need to be installed.
# Changed from the stock "smalt" to "bwa": bwa is already a required
# dependency of the production pipeline (Step 7), so this avoids installing
# a 4th mapper just for this comparison.
# use bwa (already a pipeline dependency) as the read mapper
mapper="bwa"

# Check the smalt documentation for a full explanation of options,
# including those not used by default here.
# The default options listed use the -x, -i, and -j options, which are needed for
# paired read data but should not be used for unpaired reads.
# A summary of the index options used below:
# -k sets the word (kmer) length, -s the sampling step size (i.e. is every word
#  hashed, or every second word, or one word in every 3, ...), when a hash table
# is made using the reference.
# A summary of the mapping options used below:
# -x means a read and its mate are mapped independently (not constraining them
# to be close), -y sets the minimum fraction of identical nucleotides a read
# must have to its reference before it is considered mapped, -j is the minimum
# insert size and -i the maximum insert size: outside of this range, the read
# pair is still mapped, but flagged as improperly paired.
# smalt index options (kmer length, sampling step)
smaltIndexOptions="-k 15 -s 3"
# smalt mapping options (independent mates, identity, insert range)
smaltMapOptions="-x -y 0.7 -j 0 -i 2000"

# Check the bowtie2 documentation for a full explanation of options,
# including those not used by default here.
# The default options listed use the options --maxins and --no-discordant, which
# are needed for paired read data but should not be use for unpaired reads.
# A summary of the options used below:
# --local means bowtie might soft-clip read ends if doing so maximizes the
# alignment score.
# --maxins 2000 means the maximum allowed insert size is 2000
# --no-discordant stops bowtie from looking for discordant alignments of mates
# in a pair (incorrectly oriented or exceeding the specified maximum insert
# size).
# --no-unal keeps unaligned reads out of the output (see also shiver's
# samtoolsReadFlags option below).
# --quiet means "Print nothing besides alignments and serious errors".
# bowtie2 mapping options (used only if mapper=bowtie)
bowtieOptions="--local --maxins 2000 --no-discordant --no-unal --quiet"

# Check the bwa mem documentation for a full explanation of options,
# including those not used by default here.
# A summary of the options used below:
# -v 2 sets the verbosity to "warnings and errors" but not "normal messages".
# bwa mem options (quiet: warnings/errors only)
bwaOptions='-v 2'

# After mapping, the choice of what kinds of reads should be kept is specified
# with SAM format flags, whose documentation is here:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# Flags are combined in a bitwise manner, which is fiddly. This page
# https://broadinstitute.github.io/picard/explain-flags.html
# gives a more user-friendly correspondance between SAM flags and kinds of
# reads.
# The flags used below mean unmapped reads are excluded (-F 4) and only properly
# aligned pairs are kept (-f 3). The '-f 3' should be removed for unpaired data.
# keep only properly paired, mapped reads (drop unmapped)
samtoolsReadFlags='-f 3 -F 4'

# See http://www.htslib.org/doc/samtools.html for a description of samtools
# mpileup options. Those used below mean that: the base alignment quality ('BAQ')
# calculation (described at https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtr076)
# is turned off, as seems to be appropriate for HIV
# (https://tinyurl.com/noBAQnoCry); the minimum quality for a base to be
# retained is 5 (for backward/historical consistency), and only the first
# 1000000 reads mapped to each point will be considered (NB a limit must be
# provided; the default is 250).
# Important note for data with overlapping read pairs: samtools mpileup avoids
# double counting sequence in the overlap of a read pair by setting the quality
# of all bases in the overlap to zero, for one of the two reads in the pair;
# then with any value of --min-BQ strictly greater than zero, these bases are
# effectively deleted, such that the overlap sequence is counted only once.
# If you set --min-BQ equal to zero, these bases will be counted (which is
# generally undesirable).
# mpileup options: BAQ off, min base quality 5, high depth cap
mpileupOptions='--no-BAQ --min-BQ 5 --max-depth 1000000'

# Parameters for calling the consensus base at each position:
# The minimum coverage (number of reads) to call a base instead of a '?'
# positions with fewer than 15 reads become '?' in the consensus
MinCov1=15
# The minimum coverage to use upper case for the base (to signal increased
# confidence)
# positions with >=30 reads are written upper-case (higher confidence)
MinCov2=30
# The minimum fraction of reads at a position before we call that base (or those
# bases, when one base alone does not reach that threshold fraction; e.g. say
# you have 60% A, 30% C and 10% G: if you set this fraction to 0.6 or lower we
# call an A, if you set it to 0.6-0.9 we call an M for "A or C", if you set it
# to 0.9-1 we call a V for "A, C or G".). Alternatively, if you choose a
# negative value, we always call the single most common base regardless of its
# fraction, unless two or more bases are equally (most) common, then we call the
# ambiguity code for those bases.
# negative: always call the single most common base (no ambiguity codes)
MinBaseFrac=-1

# Shall we remove read pairs marked as duplicates? i.e. using picard, for each
# set of pairs sharing the same mapped coordinates (start & end of each mate),
# keep only one pair and discard the rest? This can cause loss of diversity in
# the reads due to true biological variation as well sequencing error. We
# suggest you use this only if you understand duplication in your sequencing
# data...
# do not remove PCR/optical duplicate read pairs
deduplicate=false
# Desired command (note that MarkDuplicatesWithMateCigar exists, which may be
# better, however it still seems to have beta status; also note that you can
# include options in this command, such as a non-default
# DUPLICATE_SCORING_STRATEGY, but do not include options relating to file-naming
# or the associated shiver commands will break):
DeduplicationCommand="picard MarkDuplicates"          # command used if deduplicate=true

# Shall we remap to the consensus? (For remapping, gaps in coverage in the
# consensus will filled in by the corresponding part of the orginal reference,
# and ambiguity codes will simplified to just one of the bases they represent.
# Because of this, if remapping to the consensus, you are strongly advised to
# set the MinBaseFrac parameter above to any negative value.)
# remap reads to the first consensus for a refined final consensus
remap=true

# Shall we map contaminant reads to the reference (separately), to see which
# reads would have contaminanted our final bam file had they not been removed?
# do not separately map the removed contaminant reads
MapContaminantReads=false

# Shall we generate a version of the base frequencies file that also includes
# HXB2 coordinates (by aligning the reference used for mapping to HXB2)? Useful
# for HIV, clearly inappropriate for other viruses.
# If you are using HXB2 as the reference for mapping (instead of a reference
# constructed out of contigs as is normal for shiver), set this to false or
# there will be a problem with two identically named sequences.
# add HXB2 coordinates to the base-frequencies output (HIV-specific)
GiveHXB2coords=true

# Shall we align the contigs to the consensus, for comparison?
# do not produce a contigs-vs-consensus alignment
AlignContigsToConsensus=false

# With the default value of false, the reads in their state just before mapping
# (after any trimming of primers or adapters or low-quality bases, and after
# removal of suspected contaminant reads) will have 'temp_' prepended to their
# filenames so that they removed by the "rm temp*" command that you probably
# want to run after shiver to get rid of temporary files. Changing the value to
# true means the reads in that state don't have 'temp_' prepended to their
# filenames - handy if you want to keep them. (By request of shiver-pro Tanya!)
# let the pre-mapping reads keep the temp_ prefix (so they get cleaned up)
KeepPreMappingReads=false

# Finally, these two options are only needed for the deprecated 'fully automatic'
# version of shiver (bin/deprecated/shiver_full_auto.sh): the maximum allowed
# percentage of gaps inside contigs when aligned to their closest reference (too
# much gap content indicates misalignment, rather than deletions), and the
# minimum fraction of a contig's length that blasts to HIV.
# (deprecated auto mode) max gap fraction allowed inside a contig
MaxContigGappiness=0.05
# (deprecated auto mode) min fraction of a contig that must blast to HIV
MinContigHitFrac=0.9

# Suffixes we'll append to the sample ID for output files.
# If you change the extension (whatever follows the dot) you might break
# something.
OutputRefSuffix='_ref.fasta'                          # suffix for the per-sample mapping reference
DeduplicationStatsSuffix='_DedupStats.txt'            # suffix for the deduplication stats file
PreDeduplicationBamSuffix='_PreDedup'                 # suffix for the pre-deduplication BAM
# suffix for the mapped-contaminant-reads output
MappedContaminantReadsSuffix='_ContaminantReads'
BaseFreqsSuffix='_BaseFreqs.csv'                       # suffix for the base-frequencies table
# suffix for base freqs with global-alignment coordinates
BaseFreqsWGlobalSuffix='_BaseFreqs_ForGlobalAln.csv'
BaseFreqsWHXB2Suffix='_BaseFreqs_WithHXB2.csv'        # suffix for base freqs with HXB2 coordinates
# suffix for the insert-size distribution table
InsertSizeCountsSuffix='_InsertSizeCounts.csv'
# suffix for the coordinate-translation dictionary
CoordsDictSuffix='_coords.csv'
BlastSuffix='.blast'                                  # suffix for the contig blast output
MergedBlastSuffix='_MergedHits.blast'                 # suffix for the merged-blast-hits output
ReadsPreMapping1Suffix='_PreMapping_1.fastq'          # suffix for the pre-mapping R1 reads
ReadsPreMapping2Suffix='_PreMapping_2.fastq'          # suffix for the pre-mapping R2 reads
# suffix for the consensus prepared for global alignment
GlobalAlnSuffix='_ForGlobalAln.fasta'
BestContigToRefAlignmentSuffix='_ContigsAndBestRef.fasta' # only for fully auto.
################################################################################
# The names of temporary files we'll create in the working directory.
# If you change the extension, you may well break something.
RawContigFile1='temp_HIVcontigs_uncut1.fasta'         # temp: uncut HIV contigs (pass 1)
RawContigFile2='temp_HIVcontigs_uncut2.fasta'         # temp: uncut HIV contigs (pass 2)
# temp: contigs after cutting at blast-hit boundaries
CutContigFile='temp_HIVcontigs_cut.fasta'
# temp: contigs+refs aligned with mafft --add
TempContigAlignment1='temp_HIVcontigs_wRefs_MafftAdd.fasta'
# temp: contigs+refs aligned with mafft --addfragments
TempContigAlignment2='temp_HIVcontigs_wRefs_MafftAddFrags.fasta'
TempContigAlignment3='temp_HIVcontigs_wRefs_3.fasta'  # temp: intermediate contigs+refs alignment
TempRefAlignment='temp_RefAlignment.fasta'            # temp: reference-only alignment
# temp: reference with gaps plus extra sequence
GappyRefWithExtraSeq='temp_GappyRefWithExtraSeq.fasta'
# temp: contigs flattened to a single reference frame
FlattenedContigs='temp_FlattenedContigs.fasta'
AllContigsList='temp_AllContigsList.txt'              # temp: list of all contig names
HIVContigsListOrig='temp_HIVContigsListOrig.txt'      # temp: original list of HIV-matching contigs
HIVContigsListUser='temp_HIVContigsListUser.txt'      # temp: user-editable list of HIV contigs
ContaminantContigsList='temp_ContaminantContigsList.txt'  # temp: list of contaminant contigs
RefAndContaminantContigs='temp_RefAndContaminantContigs.fasta' # no whitespace!
BlastDB='temp_BlastDB' # no whitespace!
BadReadsBaseName='temp_ContaminantReads'              # temp: base name for contaminant read files
smaltIndex='temp_smaltRefIndex'                       # temp: smalt reference index
bowtieIndex='temp_bowtieRefIndex'                     # temp: bowtie2 reference index
# temp: contaminant reads incl. unmapped (SAM)
AllMappedContaminantReads='temp_ContaminantReads_IncUnmapped.sam'
RefFromAlignment='temp_RefFromAlignment.fasta'        # temp: reference extracted from the alignment
AllSeqsInAln='temp_AllSeqsInAln.txt'                  # temp: list of all sequences in the alignment
reads1asFasta='temp_reads1.fasta'                     # temp: R1 reads converted to FASTA
reads2asFasta='temp_reads2.fasta'                     # temp: R2 reads converted to FASTA
reads1blast1='temp_reads1_1.blast'                    # temp: R1 blast output (pass 1)
reads2blast1='temp_reads2_1.blast'                    # temp: R2 blast output (pass 1)
reads1blast2='temp_reads1_2.blast'                    # temp: R1 blast output (pass 2)
reads2blast2='temp_reads2_2.blast'                    # temp: R2 blast output (pass 2)
reads1sorted='temp_1_sorted.fastq'                    # temp: sorted R1 reads
reads2sorted='temp_2_sorted.fastq'                    # temp: sorted R2 reads
MapOutAsSam='temp_MapOut.sam'                          # temp: mapping output in SAM format
MapOutConversion1='temp_MapOutStep1'                  # temp: mapping-output conversion stage 1
MapOutConversion2='temp_MapOutStep2'                  # temp: mapping-output conversion stage 2
MapOutConversion3='temp_MapOutStep3'                  # temp: mapping-output conversion stage 3
InsertSizes1='temp_InsertSizes.txt'                   # temp: insert sizes (pass 1)
InsertSizes2='temp_InsertSizes2.txt'                  # temp: insert sizes (pass 2)
PileupFile='temp_MapOut.pileup'                       # temp: samtools pileup of the mapping
RefWithGaps='temp_RefWithGaps.fasta'                  # temp: reference with gaps re-inserted
reads1trim1='temp_reads1trim1.fastq'                  # temp: R1 reads after trim pass 1
reads1trim2='temp_reads1trim2.fastq'                  # temp: R1 reads after trim pass 2
reads2trim1='temp_reads2trim1.fastq'                  # temp: R2 reads after trim pass 1
reads2trim2='temp_reads2trim2.fastq'                  # temp: R2 reads after trim pass 2
reads1trimmings='temp_trimmings1.fastq'               # temp: sequence trimmed off R1
reads2trimmings='temp_trimmings2.fastq'               # temp: sequence trimmed off R2
# temp: alignment used when testing mafft strategies
AlignmentForTesting='temp_test.fasta'
# temp: contigs aligned with a single reference
ContigsWith1ref='temp_ContigsWith1ref.fasta'
RefMatchLog='temp_RefMatchLog.txt'                    # temp: log of contig-to-reference matches
ContigAlignmentsToRefsDir='temp_ContigAlignmentsToRefsDir' # no whitespace!
SamtoolsSortFile='temp_SamtoolsSortFile'              # temp: samtools sort scratch prefix
RefWHXB2unaln='temp_RefWHXB2unaln.fasta'              # temp: reference + HXB2, unaligned
RefWHXB2aln='temp_RefWHXB2aln.fasta'                  # temp: reference + HXB2, aligned
ContigsNoShortOnes='temp_contigs_NoShortOnes.fasta'   # temp: contigs after dropping short ones
DiscardedContigNames='temp_DiscardedContigNames.txt'  # temp: names of contigs that were discarded
