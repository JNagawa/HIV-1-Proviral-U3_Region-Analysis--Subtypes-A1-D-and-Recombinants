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
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate variant_calling

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file with strain accessions>"
    exit 1
fi

# Setting up to capture errors in commands, functions and pipelines
set -eEuo pipefail
trap 'echo "Pipeline finished at $(date)"' EXIT
trap 'echo "ERROR at $(date) in ${FUNCNAME[0]:-main}: \"$BASH_COMMAND\" at line ${LINENO}" >&2; exit 1' ERR

# logging errors in log files inside logs directory
mkdir -p logs

# Define log files
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
STDOUT_LOG="logs/pipeline_stdout_${TIMESTAMP}.log"
STDERR_LOG="logs/pipeline_stderr_${TIMESTAMP}.log"
COMBINED_LOG="logs/pipeline_combined_${TIMESTAMP}.log"

# Redirect stdout and stderr to separate log files and to a single combined log file
exec > >(tee -a "$STDOUT_LOG" "$COMBINED_LOG") 2> >(tee -a "$STDERR_LOG" "$COMBINED_LOG" >&2)


# Indicate reference organism name and accession number below
REF_ORG_NAME="Orthoebolavirus zairense"
REF_ACC="NC_002549.1"

# Store path to reference genome as it will be needed several times
REF_GENOME="data/reference/${REF_ACC}.fasta"

# Create a list array for strain accessions from the file with accs passed as cmd line args so we can be able to loop through them
mapfile -t STRAINS < $1

# Creating variant calling directory where the project will sit.
[ -d variant_calling ] || mkdir -p variant_calling
cd variant_calling

# Store path to QC directory for the strains that are being worked with
QC_DIR="qc/$(IFS=_; echo "${STRAINS[*]}")_qc"

# Define number of threads to use: If the job is run under SLURM, SLURM_CPUS_PER_TASK set from the SBATCH directives will be used /
# otherwise it defaults to one CPU per task.
THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Using $THREADS thread(s) for multi-threaded tools"

# Create file structure
mkdir -p alignments/sam alignments/unsorted_bam alignments/sorted_bam annotation checkpoints comparisons \
	 data/fastq_files data/reference data/sra_files ${QC_DIR} variants

	# Creating qc directories for each strain
	for STRAIN_ACC in ${STRAINS[@]}; do
		mkdir -p ${QC_DIR}/${STRAIN_ACC}
	done

# Define functions
download_reference () {
	# check if reference was already downloaded and exists in reference directory. If it doesn't, then it is downloaded. /
	# If it exists then the download is skipped
	local CHECKPOINT="checkpoints/${REF_ACC}_reference_download.done"

	if [ ! -f "${CHECKPOINT}" ]; then
		echo "Downloading ${REF_ACC} reference sequence..."
		wget -O "${REF_GENOME}" \
		"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${REF_ACC}&db=nuccore&report=fasta"

		touch "${CHECKPOINT}"

		echo "REFERENCE SEQUENCE DOWNLOADED✅"
	fi
}

download_fastq_files () {
	# First we download the sra file, then convert locally to avoid timeouts incase of network interruptions
	for STRAIN_ACC in ${STRAINS[@]}; do

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_fastq_download.done"

		if [ ! -f "${CHECKPOINT}" ]; then
			echo "Downloading ${STRAIN_ACC} sra files..."
			prefetch ${STRAIN_ACC} -O data/sra_files
			echo "Download of ${STRAIN_ACC} sra files complete. Starting conversion to fastq"
			fasterq-dump --split-files ${STRAIN_ACC} -O data/fastq_files/${STRAIN_ACC}

			touch "${CHECKPOINT}"
		else
			echo "${STRAIN_ACC} reads already downloaded. Skipping..."
		fi
	done
	echo "FASTQ FILES DOWNLOADED✅"
}

qc () {

	local CHECKPOINT="checkpoints/qc_$(IFS=_; echo "${STRAINS[*]}").done"

	if [ ! -f ${CHECKPOINT} ]; then
		for STRAIN_ACC in ${STRAINS[@]}; do
			echo "Performing quality checks on ${STRAIN_ACC} reads..."
			fastqc data/fastq_files/${STRAIN_ACC}/*.fastq -o ${QC_DIR}/${STRAIN_ACC}
			echo "Quality checks complete."
		done

		echo "Merging the quality checks results into a single html report..."
		multiqc ${QC_DIR} -o ${QC_DIR}
		echo "Merging quality checks complete."

		touch "${CHECKPOINT}"

		echo "QC STEPS COMPLETE✅"
	else
		echo "qc already done. Skipping..."
	fi
}

index_ref_genome () {
	# Indexing  reference genome

	echo "Indexing reference genome..."
	local CHECKPOINT="checkpoints/${REF_ACC}_ref_index.done"

	if [ ! -f "${CHECKPOINT}" ]; then
		bwa index ${REF_GENOME} # Creating index files for alignment step
		samtools faidx ${REF_GENOME}

		touch "${CHECKPOINT}"

		echo "INDEXING REFERENCE GENOME COMPLETE✅"
	else
		echo "Reference genome already indexed. Skipping..."
	fi
}


align_reads () {
	# Aligning reads for both strains by looping through list of strains
	for STRAIN_ACC in ${STRAINS[@]}; do

		echo "Aligning ${STRAIN_ACC} reads..."

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_alignment.done"

		if [ ! -f "${CHECKPOINT}" ]; then
			bwa mem -t ${THREADS} ${REF_GENOME} \
			data/fastq_files/${STRAIN_ACC}/${STRAIN_ACC}_1.fastq \
			data/fastq_files/${STRAIN_ACC}/${STRAIN_ACC}_2.fastq \
			> alignments/sam/${STRAIN_ACC}.sam

			touch "${CHECKPOINT}"
		else
			echo "${STRAIN_ACC} reads already aligned. Skipping..."
		fi
	done
	echo "ALIGNING READS COMPLETE✅"
}

sam_to_indexed_bam () {
    # Converting SAM files to BAM, optionally merging, then sorting and indexing
    for STRAIN_ACC in ${STRAINS[@]}; do

        echo "Processing ${STRAIN_ACC} SAM to indexed BAM..."

        CHECKPOINT="checkpoints/${STRAIN_ACC}_indexed_bam.done"

        if [ ! -f "${CHECKPOINT}" ]; then
            # Convert SAM to BAM
            samtools view -@ ${THREADS} -b alignments/sam/${STRAIN_ACC}.sam -o \
                alignments/unsorted_bam/${STRAIN_ACC}.bam

            # Conditional merge (if more than one BAM exists)
            BAM_COUNT=$(ls alignments/unsorted_bam/${STRAIN_ACC}*.bam 2>/dev/null | wc -l)
            if [ "$BAM_COUNT" -gt 1 ]; then
                echo "Merging BAM files for ${STRAIN_ACC}..."
                samtools merge alignments/unsorted_bam/${STRAIN_ACC}_merged.bam \
                    alignments/unsorted_bam/${STRAIN_ACC}*.bam
                merged_bam=alignments/unsorted_bam/${STRAIN_ACC}_merged.bam
            else
                echo "Only one BAM file for ${STRAIN_ACC}. Using it directly."
                merged_bam=alignments/unsorted_bam/${STRAIN_ACC}.bam
            fi

            # Sort and index
            echo "Sorting ${merged_bam}..."
            samtools sort -@ ${THREADS} "$merged_bam" -o alignments/sorted_bam/${STRAIN_ACC}.sorted.bam

            echo "Indexing alignments/sorted_bam/${STRAIN_ACC}.sorted.bam..."
            samtools index alignments/sorted_bam/${STRAIN_ACC}.sorted.bam

            touch "${CHECKPOINT}"
        else
            echo "SAM to indexed BAM for ${STRAIN_ACC} already done. Skipping..."
        fi
    done

    echo "SAM TO INDEXED_BAM CONVERSION COMPLETE ✅"
}

variant_call () {
	# Variant calling for each strain
	for STRAIN_ACC in ${STRAINS[@]}; do

		echo "Calling ${STRAIN_ACC} variants..."

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_calling.done"
		
		echo "Present working directory: $(pwd)"
		echo "Folders in current working directory: $(ls)"
		echo "Processing strain: ${STRAIN_ACC}"
		echo "Reference genome: ${REF_GENOME}"
		echo "BAM file: alignments/sorted_bam/${STRAIN_ACC}.sorted.bam"
		echo "Threads: ${THREADS}"
		if [ ! -f "${CHECKPOINT}" ]; then
			bcftools mpileup \
				-f ${REF_GENOME} \
				--threads ${THREADS} \
				-a FORMAT/DP,AD \
				-Ou \
				alignments/sorted_bam/${STRAIN_ACC}.sorted.bam | \
			bcftools call \
				-mv \
				-Oz \
				-o variants/${STRAIN_ACC}.vcf.gz

			touch "${CHECKPOINT}"
		else
			echo "${STRAIN_ACC} variants already called. Skipping..."
		fi
	done
	echo "VARIANT CALLING COMPLETE✅"
}

index_variants () {
	# Indexing variants
	for STRAIN_ACC in ${STRAINS[@]}; do

		echo "Indexing ${STRAIN_ACC} variants..."

		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_indexing.done"

		if [ ! -f "${CHECKPOINT}" ]; then
			bcftools index variants/${STRAIN_ACC}.vcf.gz

			touch "${CHECKPOINT}"
		else
			echo "Indexing ${STRAIN_ACC} variants already done. Skipping..."
		fi
	done
	echo "INDEXING VARIANTS COMPLETE✅"
}

build_variant_annotation_database () {
	# Create snpEff_data/ref_acc directory for building variant annotation database to be used in effect prediction
	# Copy reference genome file into snpEff/ref_acc directory
	# Rename 'ref_acc.fasta' to 'sequences.fa' as that's what snpEff will look for when building the database and annotating

	local CHECKPOINT1="checkpoints/${REF_ACC}_snpEff_sequences.fa_creation.done"

	echo "Creating data/snpEff_data/${REF_ACC} directory and copying ${REF_GENOME} as sequences.fa to it..."
	if [ ! -f "${CHECKPOINT1}" ]; then
		mkdir -p data/snpEff_data/${REF_ACC}
		cp ${REF_GENOME} data/snpEff_data/${REF_ACC}
		mv data/snpEff_data/${REF_ACC}/${REF_ACC}.fasta data/snpEff_data/${REF_ACC}/sequences.fa

		touch "${CHECKPOINT1}"
	else
		echo "sequences.fa already exists in data/snpEff_data/${REF_ACC} directory. Skipping..."
	fi

	# Download reference general feature file (gff) into snpEff/ref_acc for building the database and annotating

	local CHECKPOINT2="checkpoints/${REF_ACC}_snpEff_genes.gff_download.done"

	echo "Downloading ${REF_ACC} genes.gff file to data/snpEff_data/${REF_ACC} directory..."
	if [ ! -f "${CHECKPOINT2}" ]; then
		wget -O data/snpEff_data/${REF_ACC}/genes.gff \
		"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${REF_ACC}&db=nuccore&report=gff3"

		touch "${CHECKPOINT2}"
	else
		echo "genes.gff file already exists in data/snpEff_data/${REF_ACC} directory. Skipping..."
	fi

	# Creating and populating snpEff.config file
	
	echo "Creating and populating ${REF_ACC} snpEff.config file..."

	local CHECKPOINT3="checkpoints/${REF_ACC}_snpEff.config_creation.done"

	if [ ! -f "${CHECKPOINT3}" ]; then
		cat <<EOL > ./snpEff.config
        	${REF_ACC}.genome : ${REF_ORG_NAME}
        	${REF_ACC}.reference : ${REF_ACC}
EOL

		touch "${CHECKPOINT3}"
	else
		echo "${REF_ACC} snpEff.config file already exists. Skipping..."
	fi

	# Build the reference database for annotation. /
	# Turn of protein and cds checks because protein.fa and cds.fa are not available /
	# Usually unavailable for viruses.

	echo "Building ${REF_ACC} variant annotation database..."

	local CHECKPOINT4="checkpoints/${REF_ACC}_ref_variant_annotation_database_build.done"

	if [ ! -f "${CHECKPOINT4}" ]; then
		snpEff build -gff3 -v -noCheckProtein -noCheckCds -dataDir data/snpEff_data ${REF_ACC}

		touch "${CHECKPOINT4}"
	else
		echo "${REF_ACC} variant annotation database already built. Skipping..."
	fi
	echo "BUILDING ${REF_ACC} snpEff DATABASE FOR VARIANT ANNOTATION COMPLETE✅"
}

annotate_variants () {
	# Annotate variants and store annotations in annotation directory
	for STRAIN_ACC in ${STRAINS[@]}; do

		echo "Annotating ${STRAIN_ACC} variants..."
		local CHECKPOINT="checkpoints/${STRAIN_ACC}_variant_annotation.done"

		if [ ! -f "${CHECKPOINT}" ]; then
			snpEff -dataDir data/snpEff_data -v ${REF_ACC} variants/${STRAIN_ACC}.vcf.gz \
			> annotation/${STRAIN_ACC}.ann.vcf

			touch "${CHECKPOINT}"
		else
			echo "${STRAIN_ACC} variants already annotated. Skipping..."
		fi
	done
	echo "VARIANT ANNOTATION COMPLETE✅"
}

extract_missense_variants () {
	# Extract missense variants
	for STRAIN_ACC in ${STRAINS[@]}; do

		echo "Extracting ${STRAIN_ACC} missense variants..."
		local CHECKPOINT="checkpoints/${STRAIN_ACC}_missense_variants_extraction.done"

		if [ ! -f "${CHECKPOINT}" ]; then
			bcftools view -i 'ANN~"missense_variant"' annotation/${STRAIN_ACC}.ann.vcf \
			> annotation/${STRAIN_ACC}.missense.vcf

			touch "${CHECKPOINT}"
		else
			echo "Extracting ${STRAIN_ACC} missense variants already done. Skipping..."
		fi
	done
	echo "MISSENSE VARIANT EXTRACTION COMPLETE✅"
}

compress_and_index_missense_variant_files () {
	
	local CHECKPOINT="checkpoints/compressing_$(IFS=', '; echo "${STRAINS[*]}")_missense_vcfs.done"

	if [ ! -f ${CHECKPOINT} ]; then
		# Compressing and indexing missense variant files because bcftools isec requires bgzip compressed files as input
		# Compress
		for STRAIN_ACC in ${STRAINS[@]}; do
			echo "Compressing variant files..."
			bgzip annotation/${STRAIN_ACC}.missense.vcf
	
			# Index
			echo "Indexing missense variants..."
			tabix -p vcf annotation/${STRAIN_ACC}.missense.vcf.gz
		done

		touch "${CHECKPOINT}"
	else
		echo "Compressing and indexing missense variants already done. Skipping..."
	fi

		echo "COMPRESSING AND INDEXING MISSENSE VARIANTS COMPLETE✅"
}


compare_missense_variants () {

	local CHECKPOINT="checkpoints/$(IFS=_; echo "${STRAINS[*]}")_missense_variant_comparison.done"

	echo "Comparing $(IFS=', '; echo "${STRAINS[*]}") missesnse variants..."
	if [ ! -f ${CHECKPOINT} ]; then

		local OUT_DIR="$1/$(IFS=_; echo "${STRAINS[*]}")"
		mkdir -p "$OUT_DIR"

		vcfs=()
		for STRAIN_ACC in "${STRAINS[@]}"; do
			vcfs+=("annotation/${STRAIN_ACC}.missense.vcf.gz")
		done

		if (( ${#vcfs[@]} < 2 )); then
			echo "Error: need at least two VCF files to compare" >&2
			return 1
		fi

		bcftools isec -p "${OUT_DIR}" "${vcfs[@]}"

		touch "${CHECKPOINT}"
	else
		echo "$(IFS=', '; echo "${STRAINS[*]}") missense variants already compared. Skipping..."
	fi
	echo "MISSENSE VARIANT COMPARISON COMPLETE✅"
}

#===============MAIN VARIANT CALLING EXECUTION==================
# Download reference genome
download_reference

# Download fastq files
download_fastq_files

# Perform quality checks
qc

# Index reference genome
index_ref_genome

# Align reads to obtain a Sequence Alignment Map (SAM)
align_reads

# Convert SAM to BAM, merge, sort and index the BAM
sam_to_indexed_bam

# Call variants
variant_call

# Index variants
index_variants

# Build reference variant annotation database
build_variant_annotation_database

# Annotate variants
annotate_variants

# Extract missense variants
extract_missense_variants

# Compress and index missense variants
compress_and_index_missense_variant_files

# Compare missense snps between strains
compare_missense_variants comparisons

echo "READ ALIGNMENT AND, VARIANT CALLING, ANNOTATION & EFFECT PREDICTION PIPELINE COMPLETE✅✅✅"
