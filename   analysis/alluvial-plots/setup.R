library(readxl)

# Setup ------------------------------------------------------------------------
# Declare file prefixes to attach

errPoolName = "2024_06_28"

# Read in ASV and metadata tables ----------------------------------------------

# Load in ASV Table
asv = read.table("/Users/alicehong/R/16S-Analysis/data/asv-tables/all5_2024_04_29_RDP_SILVA_taxonomy.txt", head = T, row.names = NULL, sep = "")

# Load in metadata table
sample_list <- read_excel("/Users/alicehong/R/16S-Analysis/data/meta/DIABIMMUNE_all_samples_AG7_231_240_262_AG8_299.xlsx", col_types = c('Danska Mice::EndpointCode' = 'text')) # The metadata table (DIABIMMUNE + Immunophenotyping) contains a column with the character '#'which must be read as a string

## Fix mislabeled samples and remove unwanted samples --------------------------

# Fix 036 cage mouse IDs, missing week9 suffix (see AW2-76 TCAG submission form, plate 4, missing 036_NS1_week9 suffix)
# Alessandra's initial sample name has it as 2446R_R_F_S2_Pellet_week9. Not sure why there is discrepancy, but corrected to use FMPro name
colnames(asv) = gsub("NS1$", "NS1_week9", colnames(asv)) # Fixes 4 sample names missing week suffix
colnames(asv) = gsub("2446R_R_F_S2_Pellet_week9", "2446R_R_F_S2_week9", colnames(asv)) # FMPro has this as 2446R_R_F_S2_week9;

# Drop irrelevant samples 
asv = asv[,!grepl("[0-9]+FCd|_Fecal|cleanup|ABX|diabetes", names(asv), perl = TRUE)] # In this case dropping chris/tiffany's, and ABX samples

# Gets rid of extra space sample IDs
sample_list$Sample_ID = gsub(" ", "", sample_list$Sample_ID)
sample_list <- sample_list %>% distinct(Sample_ID, .keep_all = TRUE)

## Organizing ASV Table --------------------------------------------------------

# Reorders the columns in the asv table so the samples within each treatment group are together 
asv = asv %>% 
  relocate(grep("(NS6_)", names(asv), value = T, perl = T), .after = (contains("BLASTmatch"))) %>%
  relocate(grep("(NS1_)", names(asv), value = T, perl = T), .after = (contains("BLASTmatch"))) %>%
  relocate(grep("(S5_)", names(asv), value = T, perl = T), .after = (contains("BLASTmatch"))) %>%
  relocate(grep("(S2_)", names(asv), value = T, perl = T), .after = (contains("BLASTmatch")))

# Renames ASV_Num to ASV_Number (consistency)
asv <- asv %>%
  dplyr::rename(ASV_Number = OTU_Num)

# Creates extra ASV_Number_Original column for ease of tracking
asv <- asv %>%
  dplyr::mutate(ASV_Number_Original = ASV_Number, .after = "ASV_Number") 

# Sort columns by number of reads
asv = asv [order(rowSums(asv[15:ncol(asv)]),decreasing = TRUE),] 
asv = asv [order(asv $dummy_BLASTphylum),] 
asv$ASV_Number = 1:nrow(asv)
asv$ASV_Number = sprintf("%04d", asv$ASV_Number) # Makes sure the asv number is 4 digits (e.g. 0001)

## Generate the relative abundance table ---------------------------------------
asv_norm = asv
asv_norm[, 15:ncol(asv_norm)] <- sapply(asv_norm[, 15:ncol(asv_norm)],function(x) x/sum(x) * 100) 
