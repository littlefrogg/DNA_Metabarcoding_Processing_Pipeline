#!/bin/bash
#SBATCH --account=__ACCOUNT_NAME__
#SBATCH --mail-user=__EMAIL__
#SBATCH --time=__TIME__
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --mem=__MEM__
#SBATCH --job-name=blastn_alt
#SBATCH --output=blastn_alt_%j.out
#SBATCH --error=blastn_alt_%j.err

module load blast+/2.17.0

echo "Running alternative BLASTn..."
blastn -task blastn \
  -db __DB_NAME__ \
  -query __QUERY_FILE__ \
  -max_target_seqs __MAX_TARGET_SEQS__ \
  -word_size 11 \
  -evalue __EVALUE__ \
  -outfmt "6 qseqid sseqid evalue bitscore length qcovs nident pident" \
  -out __OUTPUT_FILE__

echo "BLASTn alternative run complete."
