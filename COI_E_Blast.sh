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
# These variables are set inside the script for clarity when viewing the job script file.
# The placeholders (__QUERY_FILE__, etc.) are replaced by the R MainScript.
QUERY_FILE="__QUERY_FILE__"
DB_NAME="__DB_NAME__"
OUTPUT_FILE="__OUTPUT_FILE__"
THREADS="__CPUS__"
EVALUE="__EVALUE__"

# --- 2. Load Required Modules ---
# This command may need to be adjusted based on your HPC's module system.
echo "Loading BLAST+ module..."
module load blast+

# --- 3. Run BLAST ---
echo "Starting BLAST search..."
echo "Database: ${DB_NAME}"
echo "Query: ${QUERY_FILE}"

# The -outfmt 6 command produces a standardized 12-column tabular output.
# The 'qcovs' field (Query Coverage per Subject) is added for downstream filtering.
blastn -db "${DB_NAME}" \
       -query "${QUERY_FILE}" \
       -out "${OUTPUT_FILE}" \
       -num_threads "${THREADS}" \
       -max_target_seqs 10 \
       -evalue "${EVALUE}" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"

echo "BLAST search complete. Output saved to ${OUTPUT_FILE}"