#!/bin/bash

# Specify path to fasta file 
input="/Users/alicehong/R/16S-Analysis/fasta_out_June24"

# Specify path to database
database="/Users/alicehong/GTDB_DADA2/GTDB_bac120_arc53_ssu_r214_fullTaxo_for_vsearch_updated.fa"

# Specify path/name of output file
output="/Users/alicehong/R/16S-Analysis/vsearch_output_June24"

# Specify path/name of output uc file
uc="/Users/alicehong/R/16S-Analysis/vsearch_uc_June24"

# Run usearch_global command
vsearch --usearch_global "$input" --db "$database" --userout "$output" --userfields query+target+id --uc "$uc" --id 0.993 --iddef 0 --log --uc_allhits --top_hits_only --strand both --gapopen '*'