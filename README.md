# DNA_Metabarcoding_Processing_Pipeline
Paige's DNA metabarcoding processing work in progress code
By Paige Smallman based on code written by Matthieu Leray, Helio Quintero, Luisa Meister, and Saul Fernando Rodriguez

  (parts of this code have been written or edited using generative AI, including ChatGPT-4o and Perplexity) 

This pipeline is made to be run + adjusted only through the Main Script. Additional scripts are internal and called using the source function in the Main Script.


**WORKFLOW:**

use cutadapt before running this code to trim adapters

1. Set up R environment and user-defined variables depending on project (ROHR05_12S)
2. Generate quality plots for sequencing reads
3. Adjust truncation parameters based on quality plots
4. Filter + Trim and infer Amplicon Sequence Variants (ASVs)
5. Assign Taxonomy using DADA2 Midori database; create a table of Operational Taxonomic Units (OTUs)
6. Create a phyloseq table (a tool to import, store, analyze, and graphically display complex phylogenetic sequencing data)
7. Continue processing steps based on primers used (12S vs. COI)

**Directory Organization:**

main folder:
<img width="749" alt="Screenshot 2025-04-20 at 10 53 45 AM" src="https://github.com/user-attachments/assets/02f271ef-3522-440c-a4d1-12914ae22dd8" />

inputs:

<img width="666" alt="Screenshot 2025-04-20 at 10 54 15 AM" src="https://github.com/user-attachments/assets/9207c30c-c148-4bc0-a811-981d60558c78" />

trimmed:

<img width="525" alt="Screenshot 2025-04-20 at 10 54 23 AM" src="https://github.com/user-attachments/assets/aaadb64a-198d-4a47-a01e-7d4fe028d615" />

Project with multiple primers (ex. ROHR01):

<img width="538" alt="Screenshot 2025-04-20 at 10 54 43 AM" src="https://github.com/user-attachments/assets/acca4af6-c866-4f7d-94f4-e13d2562e538" />
