#!/bin/bash
#SBATCH --account=__ACCOUNT_NAME__
#SBATCH --mail-user=__EMAIL__
#SBATCH --time=__TIME__
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --mem=__MEM__
#SBATCH --job-name=blastn
#SBATCH --output=blastn_%j.out
#SBATCH --error=blastn_%j.err

module load blast+/2.17.0

echo "Running BLASTn..."
blastn -db __DB_NAME__ \
  -query __QUERY_FILE__ \
  -num_threads __CPUS__ \
  -max_target_seqs __MAX_TARGET_SEQS__ \
  -evalue __EVALUE__ \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -out __OUTPUT_FILE__

echo "BLASTn run complete."
