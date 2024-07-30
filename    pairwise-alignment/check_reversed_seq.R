# Check if FASTA files are reverse compliment 

library(dplyr)

# Define the directory containing the files
dir_blast_output <- "/Users/alicehong/Desktop/16S_whole_from_four_primers copy/blast_output_gtdb"

# List all .fasta.txt files in the directory
files <- list.files(path = dir_blast_output, pattern = "*.fasta.txt", full.names = TRUE) 

# Function to read the first non-commented line from a file
read_first_hit <- function(file) {
  lines <- readLines(file)
  non_commented_lines <- lines[!grepl("^#", lines)]
  if (length(non_commented_lines) > 0) {
    return(non_commented_lines[1])
  } else {
    return(NA)
  }
}

# Apply the function to all files and store results in a dataframe
blast_output <- data.frame(
  file_name = basename(files),
  first_line = sapply(files, read_first_hit, USE.NAMES = FALSE),
  stringsAsFactors = FALSE
) %>% separate(first_line, into = c("taxonomy", "score", "evalue", "identical", "percent_identity", "query_length", "q_start", "q_end", "s_start"), sep = "\t", extra = "merge", fill = "right") %>% 
  mutate(q_start = as.numeric(q_start), q_end = as.numeric(q_end), s_start = as.numeric(s_start))

reversed <- blast_output %>%
  filter(s_start > q_end)
