#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00			# adjust based on expected runtime
#SBATCH --ntasks=1			# single multi-threaded task
#SBATCH --cpus-per-task=8		# number of threads for bwa/samtools/bcftools
#SBATCH --mem=32G			# adjust based on genome size. (memory per thread x threads) + buffer. -> (1-2GB /thread) /
					# Safe  for viral genomes (1-2GB) X 8 -> 8-16GB + 16GB buffer

##-------DESCRIPTION--------##
## This bash script automates comparative variant analysis workflow for two strains of a selected microorganism.
#It performs all key steps, including, downloading reference genomes and raw sequencing reads, quality checks,
# aligning reads, calling variants, annotating them, and identifying unique and shared missense SNPs between strains.
## NOTE: Ensure this script is executed on a system with Anaconda, Miniconda or an equivalent Python distribution. / 
# Create an environment using the provided .yml file to ensure that all required channels and dependencies are installed /
# before running this script.

# GROUP 14
# No.	NAME                       STUDENT NO.		REG NO.
# 1.	Nagawa Jovita              2500726007		2025/HD07/26007U
# 2.	Uwimana Yves               2500726028		2025/HDO7/26028X
# 3.	Kirunda Jeremy Menya       2500725995		2025/HD07/25995U

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"  # load conda's init script so `conda activate` works
conda activate variant_calling  # activate the env holding bwa/samtools/bcftools/snpEff

if [ "$#" -ne 1 ]; then  # require exactly one argument: the accession-list file
    echo "Usage: $0 <file with strain accessions>"  # print correct usage on wrong invocation
    exit 1  # bail out so we don't run with bad input
fi

# Setting up to capture errors in commands, functions and pipelines
set -eEuo pipefail  # strict mode: exit on error/unset var/pipe failure, inherit ERR trap
trap 'echo "Pipeline finished at $(date)"' EXIT  # always print a finish line whenever the script exits
trap 'echo "ERROR at $(date) in ${FUNCNAME[0]:-main}: \"$BASH_COMMAND\" at line ${LINENO}" >&2; exit 1' ERR  # on any error, log the location and command, then exit 1

# logging errors in log files inside logs directory
mkdir -p logs  # ensure logs dir exists before we redirect output into it

# Define log files
TIMESTAMP=$(date +%Y%m%d_%H%M%S)  # unique timestamp so each run gets its own log files
STDOUT_LOG="logs/pipeline_stdout_${TIMESTAMP}.log"  # path for this run's stdout log
STDERR_LOG="logs/pipeline_stderr_${TIMESTAMP}.log"  # path for this run's stderr log
COMBINED_LOG="logs/pipeline_combined_${TIMESTAMP}.log"  # path for the merged stdout+stderr log

# Redirect stdout and stderr to separate log files and to a single combined log file
exec > >(tee -a "$STDOUT_LOG" "$COMBINED_LOG") 2> >(tee -a "$STDERR_LOG" "$COMBINED_LOG" >&2)  # mirror stdout/stderr to the console and to the log files


# Indicate reference organism name and accession number below
REF_ORG_NAME="Orthoebolavirus zairense"  # reference organism name (used in the snpEff config)
REF_ACC="NC_002549.1"  # reference accession to download, index and annotate against

# Store path to reference genome as it will be needed several times
REF_GENOME="data/reference/${REF_ACC}.fasta"  # local path to the reference FASTA, reused throughout

# Create a list array for strain accessions from the file with accs passed as cmd line args so we can be able to loop through them
mapfile -t STRAINS < $1  # read the strain accessions from the input file into an array

# Creating variant calling directory where the project will sit.
[ -d variant_calling ] || mkdir -p variant_calling  # create the project working dir if it doesn't exist
cd variant_calling  # run everything below from inside the project dir

# Store path to QC directory for the strains that are being worked with
QC_DIR="qc/$(IFS=_; echo "${STRAINS[*]}")_qc"  # QC output dir named after the strain set being processed

# Define number of threads to use: If the job is run under SLURM, SLURM_CPUS_PER_TASK set from the SBATCH directives will be used /
# otherwise it defaults to one CPU per task.
THREADS=${SLURM_CPUS_PER_TASK:-1}  # use the SLURM CPU allocation, else default to 1 thread

echo "Using $THREADS thread(s) for multi-threaded tools"  # report the chosen thread count

# Create file structure
mkdir -p alignments/sam alignments/unsorted_bam alignments/sorted_bam annotation checkpoints comparisons \
	 data/fastq_files data/reference data/sra_files ${QC_DIR} variants  # create the whole project directory tree in one go

	# Creating qc directories for each strain
	for STRAIN_ACC in ${STRAINS[@]}; do  # make a QC subdir for each strain
		mkdir -p ${QC_DIR}/${STRAIN_ACC}  # per-strain QC output dir
	done

# Define functions
download_reference () {  # function: download the reference genome FASTA
	# check if reference was already downloaded and exists in reference directory. If it doesn't, then it is downloaded. /
	# If it exists then the download is skipped
	local CHECKPOINT="checkpoints/${REF_ACC}_reference_download.done"  # checkpoint marking the reference as already downloaded

	if [ ! -f "${CHECKPOINT}" ]; then  # only download if the checkpoint is absent
		echo "Downloading ${REF_ACC} reference sequence..."  # progress marker
		wget -O "${REF_GENOME}" \
		"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${REF_ACC}&db=nuccore&report=fasta"  # fetch the reference FASTA from NCBI

		touch "${CHECKPOINT}"  # record completion so future runs skip this

		echo "REFERENCE SEQUENCE DOWNLOADED✅"  # done marker
	fi
}

download_fastq_files () {  # function: prefetch and convert reads for every strain
	# First we download the sra file, then convert locally to avoid timeouts incase of network interruptions
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain accession

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_fastq_download.done"  # per-strain download checkpoint

		if [ ! -f "${CHECKPOINT}" ]; then  # skip if this strain was already downloaded
			echo "Downloading ${STRAIN_ACC} sra files..."  # progress marker
			prefetch ${STRAIN_ACC} -O data/sra_files  # download the SRA object first (avoids mid-convert timeouts)
			echo "Download of ${STRAIN_ACC} sra files complete. Starting conversion to fastq"  # progress marker
			fasterq-dump --split-files ${STRAIN_ACC} -O data/fastq_files/${STRAIN_ACC}  # convert the SRA to paired FASTQs locally

			touch "${CHECKPOINT}"  # mark this strain done
		else
			echo "${STRAIN_ACC} reads already downloaded. Skipping..."  # skip message for an already-downloaded strain
		fi
	done
	echo "FASTQ FILES DOWNLOADED✅"  # phase done marker
}

qc () {  # function: FastQC per strain then aggregate with MultiQC

	local CHECKPOINT="checkpoints/qc_$(IFS=_; echo "${STRAINS[*]}").done"  # checkpoint for the whole QC step

	if [ ! -f ${CHECKPOINT} ]; then  # skip if QC was already run
		for STRAIN_ACC in ${STRAINS[@]}; do  # QC each strain in turn
			echo "Performing quality checks on ${STRAIN_ACC} reads..."  # progress marker
			fastqc data/fastq_files/${STRAIN_ACC}/*.fastq -o ${QC_DIR}/${STRAIN_ACC}  # run FastQC on this strain's FASTQs
			echo "Quality checks complete."  # progress marker
		done

		echo "Merging the quality checks results into a single html report..."  # progress marker
		multiqc ${QC_DIR} -o ${QC_DIR}  # merge all FastQC reports into one HTML
		echo "Merging quality checks complete."  # progress marker

		touch "${CHECKPOINT}"  # mark QC done

		echo "QC STEPS COMPLETE✅"  # phase done marker
	else
		echo "qc already done. Skipping..."  # skip message when QC is already done
	fi
}

index_ref_genome () {  # function: build BWA and samtools indexes on the reference
	# Indexing  reference genome

	echo "Indexing reference genome..."  # progress marker
	local CHECKPOINT="checkpoints/${REF_ACC}_ref_index.done"  # checkpoint marking the reference as indexed

	if [ ! -f "${CHECKPOINT}" ]; then  # skip if already indexed
		bwa index ${REF_GENOME} # Creating index files for alignment step
		samtools faidx ${REF_GENOME}  # build the .fai FASTA index (needed by bcftools)

		touch "${CHECKPOINT}"  # mark indexing done

		echo "INDEXING REFERENCE GENOME COMPLETE✅"  # done marker
	else
		echo "Reference genome already indexed. Skipping..."  # skip message when already indexed
	fi
}


align_reads () {  # function: align each strain's reads to the reference with BWA-MEM
	# Aligning reads for both strains by looping through list of strains
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

		echo "Aligning ${STRAIN_ACC} reads..."  # progress marker

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_alignment.done"  # per-strain alignment checkpoint

		if [ ! -f "${CHECKPOINT}" ]; then  # skip if already aligned
			bwa mem -t ${THREADS} ${REF_GENOME} \
			data/fastq_files/${STRAIN_ACC}/${STRAIN_ACC}_1.fastq \
			data/fastq_files/${STRAIN_ACC}/${STRAIN_ACC}_2.fastq \
			> alignments/sam/${STRAIN_ACC}.sam  # align the paired FASTQs to the reference, writing a SAM

			touch "${CHECKPOINT}"  # mark this strain aligned
		else
			echo "${STRAIN_ACC} reads already aligned. Skipping..."  # skip message for an already-aligned strain
		fi
	done
	echo "ALIGNING READS COMPLETE✅"  # phase done marker
}

sam_to_indexed_bam () {  # function: SAM -> BAM, optional merge, then sort and index
    # Converting SAM files to BAM, optionally merging, then sorting and indexing
    for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

        echo "Processing ${STRAIN_ACC} SAM to indexed BAM..."  # progress marker

        CHECKPOINT="checkpoints/${STRAIN_ACC}_indexed_bam.done"  # per-strain checkpoint for this conversion

        if [ ! -f "${CHECKPOINT}" ]; then  # skip if already done
            # Convert SAM to BAM
            samtools view -@ ${THREADS} -b alignments/sam/${STRAIN_ACC}.sam -o \
                alignments/unsorted_bam/${STRAIN_ACC}.bam  # convert the SAM to an unsorted BAM

            # Conditional merge (if more than one BAM exists)
            BAM_COUNT=$(ls alignments/unsorted_bam/${STRAIN_ACC}*.bam 2>/dev/null | wc -l)  # count BAMs to decide whether a merge is needed
            if [ "$BAM_COUNT" -gt 1 ]; then  # more than one BAM present for this strain?
                echo "Merging BAM files for ${STRAIN_ACC}..."  # progress marker
                samtools merge alignments/unsorted_bam/${STRAIN_ACC}_merged.bam \
                    alignments/unsorted_bam/${STRAIN_ACC}*.bam  # merge the multiple BAMs into one
                merged_bam=alignments/unsorted_bam/${STRAIN_ACC}_merged.bam  # point downstream steps at the merged BAM
            else
                echo "Only one BAM file for ${STRAIN_ACC}. Using it directly."  # single-BAM message
                merged_bam=alignments/unsorted_bam/${STRAIN_ACC}.bam  # use the single BAM directly
            fi

            # Sort and index
            echo "Sorting ${merged_bam}..."  # progress marker
            samtools sort -@ ${THREADS} "$merged_bam" -o alignments/sorted_bam/${STRAIN_ACC}.sorted.bam  # coordinate-sort the BAM

            echo "Indexing alignments/sorted_bam/${STRAIN_ACC}.sorted.bam..."  # progress marker
            samtools index alignments/sorted_bam/${STRAIN_ACC}.sorted.bam  # build the .bai index for random access

            touch "${CHECKPOINT}"  # mark this strain done
        else
            echo "SAM to indexed BAM for ${STRAIN_ACC} already done. Skipping..."  # skip message when already done
        fi
    done

    echo "SAM TO INDEXED_BAM CONVERSION COMPLETE ✅"  # phase done marker
}

variant_call () {  # function: call variants per strain with bcftools
	# Variant calling for each strain
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

		echo "Calling ${STRAIN_ACC} variants..."  # progress marker

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_calling.done"  # per-strain variant-calling checkpoint
		
		echo "Present working directory: $(pwd)"  # debug: show the current working directory
		echo "Folders in current working directory: $(ls)"  # debug: list the current directory contents
		echo "Processing strain: ${STRAIN_ACC}"  # debug: which strain is being processed
		echo "Reference genome: ${REF_GENOME}"  # debug: reference genome path
		echo "BAM file: alignments/sorted_bam/${STRAIN_ACC}.sorted.bam"  # debug: input BAM path
		echo "Threads: ${THREADS}"  # debug: thread count
		if [ ! -f "${CHECKPOINT}" ]; then  # skip if variants were already called
			bcftools mpileup \
				-f ${REF_GENOME} \
				--threads ${THREADS} \
				-a FORMAT/DP,AD \
				-Ou \
				alignments/sorted_bam/${STRAIN_ACC}.sorted.bam | \
			bcftools call \
				-mv \
				-Oz \
				-o variants/${STRAIN_ACC}.vcf.gz  # pile up reads and call variants into a compressed VCF

			touch "${CHECKPOINT}"  # mark this strain done
		else
			echo "${STRAIN_ACC} variants already called. Skipping..."  # skip message when already called
		fi
	done
	echo "VARIANT CALLING COMPLETE✅"  # phase done marker
}

index_variants () {  # function: index each strain's VCF
	# Indexing variants
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

		echo "Indexing ${STRAIN_ACC} variants..."  # progress marker

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_indexing.done"  # per-strain indexing checkpoint

		if [ ! -f "${CHECKPOINT}" ]; then  # skip if already indexed
			bcftools index variants/${STRAIN_ACC}.vcf.gz  # build the VCF index (.csi)

			touch "${CHECKPOINT}"  # mark this strain done
		else
			echo "Indexing ${STRAIN_ACC} variants already done. Skipping..."  # skip message when already indexed
		fi
	done
	echo "INDEXING VARIANTS COMPLETE✅"  # phase done marker
}

build_variant_annotation_database () {  # function: stage files and build the snpEff annotation database
	# Create snpEff_data/ref_acc directory for building variant annotation database to be used in effect prediction
	# Copy reference genome file into snpEff/ref_acc directory
	# Rename 'ref_acc.fasta' to 'sequences.fa' as that's what snpEff will look for when building the database and annotating

	local CHECKPOINT1="checkpoints/${REF_ACC}_snpEff_sequences.fa_creation.done"  # checkpoint for staging sequences.fa

	echo "Creating data/snpEff_data/${REF_ACC} directory and copying ${REF_GENOME} as sequences.fa to it..."  # progress marker
	if [ ! -f "${CHECKPOINT1}" ]; then  # skip if already staged
		mkdir -p data/snpEff_data/${REF_ACC}  # snpEff expects a per-genome data directory
		cp ${REF_GENOME} data/snpEff_data/${REF_ACC}  # copy the reference FASTA into it
		mv data/snpEff_data/${REF_ACC}/${REF_ACC}.fasta data/snpEff_data/${REF_ACC}/sequences.fa  # rename to sequences.fa (the name snpEff looks for)

		touch "${CHECKPOINT1}"  # mark this sub-step done
	else
		echo "sequences.fa already exists in data/snpEff_data/${REF_ACC} directory. Skipping..."  # skip message when already staged
	fi

	# Download reference general feature file (gff) into snpEff/ref_acc for building the database and annotating

	local CHECKPOINT2="checkpoints/${REF_ACC}_snpEff_genes.gff_download.done"  # checkpoint for the genes.gff download

	echo "Downloading ${REF_ACC} genes.gff file to data/snpEff_data/${REF_ACC} directory..."  # progress marker
	if [ ! -f "${CHECKPOINT2}" ]; then  # skip if already downloaded
		wget -O data/snpEff_data/${REF_ACC}/genes.gff \
		"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${REF_ACC}&db=nuccore&report=gff3"  # fetch the GFF3 annotation snpEff needs

		touch "${CHECKPOINT2}"  # mark this sub-step done
	else
		echo "genes.gff file already exists in data/snpEff_data/${REF_ACC} directory. Skipping..."  # skip message when already present
	fi

	# Creating and populating snpEff.config file
	
	echo "Creating and populating ${REF_ACC} snpEff.config file..."  # progress marker

	local CHECKPOINT3="checkpoints/${REF_ACC}_snpEff.config_creation.done"  # checkpoint for creating the config

	if [ ! -f "${CHECKPOINT3}" ]; then  # skip if the config already exists
		cat <<EOL > ./snpEff.config  # write the genome entry so snpEff can find this DB
        	${REF_ACC}.genome : ${REF_ORG_NAME}
        	${REF_ACC}.reference : ${REF_ACC}
EOL

		touch "${CHECKPOINT3}"  # mark this sub-step done
	else
		echo "${REF_ACC} snpEff.config file already exists. Skipping..."  # skip message when config exists
	fi

	# Build the reference database for annotation. /
	# Turn of protein and cds checks because protein.fa and cds.fa are not available /
	# Usually unavailable for viruses.

	echo "Building ${REF_ACC} variant annotation database..."  # progress marker

	local CHECKPOINT4="checkpoints/${REF_ACC}_ref_variant_annotation_database_build.done"  # checkpoint for the database build

	if [ ! -f "${CHECKPOINT4}" ]; then  # skip if already built
		snpEff build -gff3 -v -noCheckProtein -noCheckCds -dataDir data/snpEff_data ${REF_ACC}  # build the annotation DB (protein/CDS checks off for viruses)

		touch "${CHECKPOINT4}"  # mark this sub-step done
	else
		echo "${REF_ACC} variant annotation database already built. Skipping..."  # skip message when already built
	fi
	echo "BUILDING ${REF_ACC} snpEff DATABASE FOR VARIANT ANNOTATION COMPLETE✅"  # phase done marker
}

annotate_variants () {  # function: annotate each strain's variants with snpEff
	# Annotate variants and store annotations in annotation directory
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

		echo "Annotating ${STRAIN_ACC} variants..."  # progress marker
		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_annotation.done"  # per-strain annotation checkpoint

		if [ ! -f "${CHECKPOINT}" ]; then  # skip if already annotated
			snpEff -dataDir data/snpEff_data -v ${REF_ACC} variants/${STRAIN_ACC}.vcf.gz \
			> annotation/${STRAIN_ACC}.ann.vcf  # annotate variants with predicted functional effects

			touch "${CHECKPOINT}"  # mark this strain done
		else
			echo "${STRAIN_ACC} variants already annotated. Skipping..."  # skip message when already annotated
		fi
	done
	echo "VARIANT ANNOTATION COMPLETE✅"  # phase done marker
}

extract_missense_variants () {  # function: extract missense SNPs from each strain's annotations
	# Extract missense variants
	for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain

		echo "Extracting ${STRAIN_ACC} missense variants..."  # progress marker
		local CHECKPOINT="checkpoints/${STRAIN_ACC}_missense_variants_extraction.done"  # per-strain extraction checkpoint

		if [ ! -f "${CHECKPOINT}" ]; then  # skip if already extracted
			bcftools view -i 'ANN~"missense_variant"' annotation/${STRAIN_ACC}.ann.vcf \
			> annotation/${STRAIN_ACC}.missense.vcf  # keep only records annotated as missense variants

			touch "${CHECKPOINT}"  # mark this strain done
		else
			echo "Extracting ${STRAIN_ACC} missense variants already done. Skipping..."  # skip message when already extracted
		fi
	done
	echo "MISSENSE VARIANT EXTRACTION COMPLETE✅"  # phase done marker
}

compress_and_index_missense_variant_files () {  # function: bgzip and tabix-index the missense VCFs
	
	local CHECKPOINT="checkpoints/compressing_$(IFS=', '; echo "${STRAINS[*]}")_missense_vcfs.done"  # checkpoint for the compress/index step

	if [ ! -f ${CHECKPOINT} ]; then  # skip if already done
		# Compressing and indexing missense variant files because bcftools isec requires bgzip compressed files as input
		# Compress
		for STRAIN_ACC in ${STRAINS[@]}; do  # loop over each strain
			echo "Compressing variant files..."  # progress marker
			bgzip annotation/${STRAIN_ACC}.missense.vcf  # bgzip the VCF (bcftools isec needs bgzipped input)
	
			# Index
			echo "Indexing missense variants..."  # progress marker
			tabix -p vcf annotation/${STRAIN_ACC}.missense.vcf.gz  # build the tabix index for the bgzipped VCF
		done

		touch "${CHECKPOINT}"  # mark this step done
	else
		echo "Compressing and indexing missense variants already done. Skipping..."  # skip message when already done
	fi

		echo "COMPRESSING AND INDEXING MISSENSE VARIANTS COMPLETE✅"  # phase done marker
}


compare_missense_variants () {  # function: intersect strains' missense VCFs with bcftools isec

	local CHECKPOINT="checkpoints/$(IFS=_; echo "${STRAINS[*]}")_missense_variant_comparison.done"  # checkpoint for the comparison step

	echo "Comparing $(IFS=', '; echo "${STRAINS[*]}") missesnse variants..."  # progress marker
	if [ ! -f ${CHECKPOINT} ]; then  # skip if already compared

		local OUT_DIR="$1/$(IFS=_; echo "${STRAINS[*]}")"  # output dir under the passed-in base directory ($1)
		mkdir -p "$OUT_DIR"  # ensure the output dir exists

		vcfs=()  # collect the per-strain VCF paths
		for STRAIN_ACC in "${STRAINS[@]}"; do  # build the list of VCFs to intersect
			vcfs+=("annotation/${STRAIN_ACC}.missense.vcf.gz")  # add this strain's missense VCF
		done

		if (( ${#vcfs[@]} < 2 )); then  # bcftools isec needs at least two inputs
			echo "Error: need at least two VCF files to compare" >&2  # error message on stderr
			return 1  # abort the comparison
		fi

		bcftools isec -p "${OUT_DIR}" "${vcfs[@]}"  # intersect the VCFs into shared/unique record sets

		touch "${CHECKPOINT}"  # mark this step done
	else
		echo "$(IFS=', '; echo "${STRAINS[*]}") missense variants already compared. Skipping..."  # skip message when already compared
	fi
	echo "MISSENSE VARIANT COMPARISON COMPLETE✅"  # phase done marker
}

#===============MAIN VARIANT CALLING EXECUTION==================
# Download reference genome
download_reference  # step 1: fetch the reference genome

# Download fastq files
download_fastq_files  # step 2: download and convert the reads

# Perform quality checks
qc  # step 3: quality control

# Index reference genome
index_ref_genome  # step 4: index the reference

# Align reads to obtain a Sequence Alignment Map (SAM)
align_reads  # step 5: align reads -> SAM

# Convert SAM to BAM, merge, sort and index the BAM
sam_to_indexed_bam  # step 6: SAM -> sorted, indexed BAM

# Call variants
variant_call  # step 7: call variants

# Index variants
index_variants  # step 8: index the VCFs

# Build reference variant annotation database
build_variant_annotation_database  # step 9: build the snpEff annotation database

# Annotate variants
annotate_variants  # step 10: annotate the variants

# Extract missense variants
extract_missense_variants  # step 11: extract missense SNPs

# Compress and index missense variants
compress_and_index_missense_variant_files  # step 12: bgzip and index the missense VCFs

# Compare missense snps between strains
compare_missense_variants comparisons  # step 13: compare strains (results under comparisons/)

echo "READ ALIGNMENT AND, VARIANT CALLING, ANNOTATION & EFFECT PREDICTION PIPELINE COMPLETE✅✅✅"  # whole pipeline done marker
