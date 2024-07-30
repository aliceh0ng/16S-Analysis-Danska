## Read in FASTA sequences from 16S folder -------------------------------------

library(tibble)
library(Biostrings)
library(purrr)
library(tidyr)
library(dplyr)
library(tibble)

# Reading in fasta files from WGS

fasta <- function(folder_path) {
  # List all fasta files in the folder
  fasta_files <- list.files(path = folder_path, pattern = "\\.fasta$", full.names = TRUE)
  
  # Initialize an empty list to store data
  fasta_data <- list()
  
  # Loop through each file and read its contents
  for (file in fasta_files) {
    fasta_content <- readDNAStringSet(file)
    fasta_names <- names(fasta_content)
    fasta_sequences <- as.character(fasta_content)
    
    # Create a temporary dataframe for each file
    temp_df <- data.frame(Name = fasta_names, Sequence = fasta_sequences, stringsAsFactors = FALSE)
    
    # Append to the list
    fasta_data[[file]] <- temp_df
  }
  
  # Combine all dataframes into one
  combined_df <- do.call(rbind, fasta_data)
  
  return(combined_df)
}

# Perform function to read fasta sequences 

folder_path <- "/Users/alicehong/R/16S-Analysis/      data/16S_whole_from_four_primers"
fasta.df <- fasta(folder_path)
  
fasta.df$Number <- as.numeric(gsub("[^0-9]", "", fasta.df$Name))
fasta.df$ID <- sprintf("%03d", fasta.df$Number)

fasta.df <- fasta.df %>% select("ID", "Sequence")

## July 4th, 2024
# ONLY RUN THIS WITH DANIEL'S FULL LENGTH SEQUENCES (some are in wrong direction)
# Reverse complement for ID: 
#   NS1: 010, 016, 017, 018, 020, 115, 028
#   NS6: 029, 031, 032, 034, 117, 036, 037, 039, 041, 043, 044, 046, 047, 050, 051, 053, 060
#   S2: 069, 082, 118, 086, 087, 088
#   S5: 090, 101, 108, 120

reverse_ID <- c("010", "016", "017", "018", "020", "115", "028", "029", "031", "032", "034", "117", 
                "036", "037", "039", "041", "043", "044", "046", "047", "050", "051", "053", "060", 
                "069", "082", "118", "086", "087", "088", "090", "101", "108", "120")

# Convert dataframe to DNAStringSet
sequences <- DNAStringSet(fasta.df$Sequence)
names(sequences) <- fasta.df$ID

# Identify sequences that need to be reversed and replace them with their reverse complements
reversed_sequences <- reverseComplement(sequences[names(sequences) %in% reverse_ID])
sequences[names(sequences) %in% reverse_ID] <- reversed_sequences

# Put back into dataframe
fasta.df$Sequence <- as.character(sequences)

## Final Sequence Dataframe ----------------------------------------------------

all <- read.xlsx("/Users/alicehong/R/16S-Analysis/      data/meta/WGS_Samples_June28.xlsx", sheet = "All")
all <- all %>%   mutate(ID = str_extract(Sample_ID, "\\d{3}$"))

# Left join: Keep all rows from df1, and only matching rows from df2
WGS <- left_join(all, fasta.df, by = "ID")

WGS_community <- WGS %>%
  group_by(mouse_inocula) %>%
  group_split()

# reference sequences to align to

NS1_ref <- DNAStringSet(WGS_community[[1]]$Sequence)
names(NS1_ref) <- as.character(seq_along(NS1_ref))

NS6_ref <- DNAStringSet(WGS_community[[2]]$Sequence)
names(NS6_ref) <- as.character(seq_along(NS6_ref))

S2_ref <- DNAStringSet(WGS_community[[3]]$Sequence)
names(S2_ref) <- as.character(seq_along(S2_ref))

S5_ref <- DNAStringSet(WGS_community[[4]]$Sequence)
names(S5_ref) <- as.character(seq_along(S5_ref))

## ASV Table aka input seqs ----------------------------------------------------

asv = read.table("/Users/alicehong/R/16S-Analysis/      data/asv-tables/2024_06_24_VSEARCH_GTDB_RDP_Silva_taxonomy_no_DECIPHER_cleanedup.txt", head = T, row.names = NULL, sep = "")
asv.temp <- asv %>% select(-contains("dummy")) %>% select(-contains("OTU"))
asv.temp <- column_to_rownames(asv.temp, var = "sequences")

communities <- c("NS1", "NS6", "S2", "S5")
asv.list <- list()

for (community in communities) {
  temp <- asv.temp %>% select(contains(community))
  temp <- temp[rowSums(temp) != 0, ]
  asv.list[[community]] <- temp
}

NS1_asv <-DNAStringSet(rownames(asv.list[["NS1"]]))
names(NS1_asv) <- as.character(seq_along(NS1_asv))

NS6_asv <- DNAStringSet(rownames(asv.list[["NS6"]]))
names(NS6_asv) <- as.character(seq_along(NS6_asv))

S2_asv <- DNAStringSet(rownames(asv.list[["S2"]]))
names(S2_asv) <- as.character(seq_along(S2_asv))

S5_asv <- DNAStringSet(rownames(asv.list[["S5"]]))
names(S5_asv) <- as.character(seq_along(S5_asv))

## MSA -------------------------------------------------------------------------

# Initialize a list to store alignment results
# alignment_results <- expand.grid(input_index = seq_along(NS6_asv), reference_index = seq_along(NS6_ref)) %>%
#   pmap(function(input_index, reference_index) {
#     alignment <- pairwiseAlignment(NS6_asv[input_index], NS6_ref[reference_index], type = "local")
#     list(name = paste0("Input_", input_index, "_Ref_", reference_index), alignment = alignment)
#   }) %>%
#   set_names(map_chr(., "name"))
# 
# top_hits <- map(names(NS6_asv), function(i) {
#   relevant_alignments <- keep(alignment_results, ~ grepl(paste0("Input_", i, "_"), .x$name))
#   scores <- map_dbl(relevant_alignments, ~ score(.x$alignment))
#   best_index <- which.max(scores)
#   best_alignment <- relevant_alignments[[best_index]]
#   list(name = best_alignment$name, alignment = best_alignment$alignment)
# }) %>%
#   set_names(map_chr(., "name"))
# 
# top_hits_df <- map_dfr(top_hits, ~ tibble(name = .x$name, score = .x$alignment@score))

### Making this into a function! -----------------------------------------------

perform_alignment <- function(asv_set, ref_set) {
  alignment_results <- expand.grid(input_index = seq_along(asv_set), reference_index = seq_along(ref_set)) %>%
    pmap(function(input_index, reference_index) {
      alignment <- pairwiseAlignment(asv_set[input_index], ref_set[reference_index], type = "local")
      list(name = paste0("Input_", input_index, "_Ref_", reference_index), alignment = alignment)
    }) %>%
    set_names(map_chr(., "name"))
  
  top_hits <- map(seq_along(asv_set), function(i) {
    relevant_alignments <- keep(alignment_results, ~ grepl(paste0("Input_", i, "_"), .x$name))
    scores <- map_dbl(relevant_alignments, ~ score(.x$alignment))
    best_index <- which.max(scores)
    best_alignment <- relevant_alignments[[best_index]]
    list(name = best_alignment$name, alignment = best_alignment$alignment)
  }) %>%
    set_names(map_chr(., "name"))
  
  # Create a dataframe with the name and score
  top_hits_df <- map_dfr(top_hits, ~ tibble(name = .x$name, score = .x$alignment@score))
  
  #return(top_hits_df)
  return(list(alignment_results = alignment_results, top_hits = top_hits, top_hits_df = top_hits_df))
}

# Run the function for each pair of datasets
results_NS6 <- perform_alignment(NS6_asv, NS6_ref)
results_NS1 <- perform_alignment(NS1_asv, NS1_ref)
results_S2 <- perform_alignment(S2_asv, S2_ref)
results_S5 <- perform_alignment(S5_asv, S5_ref)

## July 4th: ran alignment again with the reverse strands corrected, saving old results
old_top_hits_NS6 <- top_hits_NS6
old_top_hits_NS1 <- top_hits_NS1
old_top_hits_S2 <- top_hits_S2
old_top_hits_S5 <- top_hits_S5
old_all_top_hits <- all_top_hits

# Save alignment_results and top_hits for each dataset
top_hits_NS6 <- results_NS6$top_hits_df %>% mutate(seq = rownames(asv.list[["NS6"]])) %>%
  separate(name, into = c('Input', 'asv_num', 'Ref', 'ref_num'), sep = "_") %>%
  select(-Input, -Ref)

top_hits_NS1 <- results_NS1$top_hits_df %>% mutate(seq = rownames(asv.list[["NS1"]])) %>% 
  separate(name, into = c('Input', 'asv_num', 'Ref', 'ref_num'), sep = "_") %>%
  select(-Input, -Ref)

top_hits_S2 <- results_S2$top_hits_df %>% mutate(seq = rownames(asv.list[["S2"]])) %>% 
  separate(name, into = c('Input', 'asv_num', 'Ref', 'ref_num'), sep = "_") %>%
  select(-Input, -Ref)

top_hits_S5 <- results_S5$top_hits_df %>% mutate(seq = rownames(asv.list[["S5"]])) %>% 
  separate(name, into = c('Input', 'asv_num', 'Ref', 'ref_num'), sep = "_") %>%
  select(-Input, -Ref)

# Combine all top_hits results into one dataframe (optional)
all_top_hits <- bind_rows(
  top_hits_NS6 %>% mutate(set = "NS6"),
  top_hits_NS1 %>% mutate(set = "NS1"),
  top_hits_S2 %>% mutate(set = "S2"),
  top_hits_S5 %>% mutate(set = "S5")
)

## Collapse ASV table (for each community) -------------------------------------

# # Prepare reference asv data
# ref.S5 <- WGS_community[[4]] %>% select(-ID) %>% rownames_to_column(var = "ref_num")
# 
# # Create new asv table with reference asv ID 
# new.asv.S5 <- asv.list[["S5"]] %>% rownames_to_column(var = "seq") %>% 
#   left_join(top_hits_S5, by = "seq") %>% select(-score, -asv_num) %>% 
#   column_to_rownames(var = "seq") %>%
#   
#   # Collapse rows with same reference asv ID
#   group_by(ref_num) %>%
#   summarise(across(everything(), sum, .names = "sum_{col}")) %>%
#   ungroup() %>%
# 
#   # Add in reference asv sequence
#   left_join(ref.NS1, by = "ref_num") %>%
#   select(ref_num, everything(), -starts_with("sum_"), starts_with("sum_"))

## Making into a function!!

newASV <- function(ref_data, old_asv, top_hits) {
  
  # Prepare reference asv data
  ref <- ref_data %>% select(-ID) %>% rownames_to_column(var = "ref_num")
  
  # Create new asv table with reference asv ID 
  new.asv <- old_asv %>% rownames_to_column(var = "seq") %>% 
    left_join(top_hits, by = "seq") %>% select(-score, -asv_num) %>% # Join with alignment results
    column_to_rownames(var = "seq") %>%
    
    # Collapse rows with same reference asv ID
    group_by(ref_num) %>%
    summarise(across(everything(), sum, .names = "sum_{col}")) %>%
    ungroup() %>%
    
    # Add in reference asv sequence
    left_join(ref, by = "ref_num") %>%
    select(ref_num, everything(), -starts_with("sum_"), starts_with("sum_"))
  
  return(new.asv)
}

new.asv.NS1 <- newASV(WGS_community[[1]], asv.list[["NS1"]], top_hits_NS1)
new.asv.NS6 <- newASV(WGS_community[[2]], asv.list[["NS6"]], top_hits_NS6)
new.asv.S2 <- newASV(WGS_community[[3]], asv.list[["S2"]], top_hits_S2)
new.asv.S5 <- newASV(WGS_community[[4]], asv.list[["S5"]], top_hits_S5)


# -------------
  
NS1_row_sums <- rowSums(new.asv.NS1[,8:127]) %>% cbind(new.asv.NS1[,1:7])
NS6_row_sums <- rowSums(new.asv.NS6[,8:42]) %>% cbind(new.asv.NS6[,1:7])

S2_row_sums <- rowSums(new.asv.S2[,8:144]) %>% cbind(new.asv.S2[,1:7])
S5_row_sums <- rowSums(new.asv.S5[,8:42]) %>% cbind(new.asv.S5[,1:7])

# View the row sums
print(row_sums)

# Add the row sums to the original matrix (optional)
counts_with_row_sums <- cbind(counts, row_sums)

