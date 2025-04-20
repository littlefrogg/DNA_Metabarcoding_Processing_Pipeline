# eDNA_WiP
Paige's eDNA processing work in progress code

This pipeline is made to be run + adjusted only through the Main Script. Additional scripts are internal and called using the source function in the Main Script.

WORKFLOW:
(0). use cutadapt before running this code to trim adapters
1. Set up R environment and user-defined variables depending on project (ROHR05_12S)
2. Generate quality plots for sequencing reads
3. Adjust truncation parameters based on quality plots
4. Filter + Trim and infer Amplicon Sequence Variants (ASVs)
5. Assign Taxonomy using DADA2 Midori database; create a table of Operational Taxonomic Units (OTUs)
6. Create a phyloseq table (a tool to import, store, analyze, and graphically display complex phylogenetic sequencing data)
7. Continue processing steps based on primers used (12S vs. COI)

Directory Organization:
<img width="749" alt="Screenshot 2025-04-20 at 10 53 45 AM" src="https://github.com/user-attachments/assets/02f271ef-3522-440c-a4d1-12914ae22dd8" />

