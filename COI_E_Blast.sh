#!/bin/bash
#SBATCH --job-name=COI_BLAST
#SBATCH --account=__ACCOUNT_NAME__
#SBATCH --time=__TIME__
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --mem=__MEM__
#SBATCH --output=blast_job_%j.out
#SBATCH --error=blast_job_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__

# --- 1. Set Variables (dynamically set by R script) ---
QUERY_FILE="__QUERY_FILE__"
DB_NAME="__DB_NAME__"
OUTPUT_FILE="__OUTPUT_FILE__"
THREADS="__CPUS__"
EVALUE="__EVALUE__"

# --- 2. Load Required Modules ---
echo "Loading BLAST module..."
module load blast+

# --- 3. Run BLAST ---
echo "Starting BLAST search..."
blastn -db "${DB_NAME}" \
       -query "${QUERY_FILE}" \
       -out "${OUTPUT_FILE}" \
       -num_threads "${THREADS}" \
       -max_target_seqs 10 \
       -evalue "${EVALUE}" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

echo "BLAST search complete. Output saved to ${OUTPUT_FILE}"