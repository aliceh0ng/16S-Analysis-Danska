
# Phyloseq ---------------------------------------------------------------------

## Create/organize meta data table ---------------------------------------------
meta.asv <- as.data.table(grep("ASV|BLAST|sequences|per_id|bit|evalue", names(asv), value = TRUE, invert = "TRUE"))
meta.asv = meta.asv %>% 
  dplyr::rename(Sample = V1) %>%
  mutate(Sample_ID = Sample)

meta.asv$Sample_ID = gsub("^X", "", meta.asv$Sample_ID)

meta.asv = meta.asv %>%
  mutate(Sample_ID2 = Sample_ID)

meta.asv = merge(meta.asv, sample_list, by.x = "Sample_ID2", by.y = "Sample_ID", all.x = T)

meta.asv$Sample_ID2 = NULL

meta.asv$Plate = NA

plate1 = row.names(seqtab.nochim.plate1)
plate2 = row.names(seqtab.nochim.plate2)
plate3 = row.names(seqtab.nochim.plate3)
plate4 = row.names(seqtab.nochim.plate4)
plate5 = row.names(seqtab.nochim.plate5)

meta.asv$Plate[meta.asv$Sample_ID %in% plate5] = "5"
meta.asv$Plate[meta.asv$Sample_ID %in% plate3] = "4"
meta.asv$Plate[meta.asv$Sample_ID %in% plate2] = "3"
meta.asv$Plate[meta.asv$Sample_ID %in% plate4] = "2"
meta.asv$Plate[meta.asv$Sample_ID %in% plate1] = "1"

# Fixing controls and appending plate number
meta.asv$Plate[grep("plate1$", meta.asv$Sample_ID)] = "1"
meta.asv$Plate[grep("plate2$", meta.asv$Sample_ID)] = "2"
meta.asv$Plate[grep("plate3$", meta.asv$Sample_ID)] = "3"
meta.asv$Plate[grep("plate4$", meta.asv$Sample_ID)] = "4"
meta.asv$Plate[grep("plate5$", meta.asv$Sample_ID)] = "5"

## This is necessary for sample names that were changed on the asv table (seqtab stil has original name)
meta.asv$Plate[grep("036R", meta.asv$Sample_ID)] = "4"
meta.asv$Plate[grep("2446R_R_F_S2_Pellet_week9", meta.asv$Sample_ID)] = "3"
meta.asv$Plate[grep("2446R_R_F_S2_week9 ", meta.asv$Sample_ID)] = "3" # meta.asv[149,51] <- "3"
meta.asv$Plate[grep("T1D_", meta.asv$Sample_ID)] = "1"

## TODO: Write a loop for this later
meta.asv[ , 'Week'] = NA
meta.asv$Week[grep("week5", meta.asv$Sample_ID)] = "5"
meta.asv$Week[grep("week6", meta.asv$Sample_ID)] = "6"
meta.asv$Week[grep("week7", meta.asv$Sample_ID)] = "7"
meta.asv$Week[grep("week9", meta.asv$Sample_ID)] = "9"
meta.asv$Week[grep("week10", meta.asv$Sample_ID)] = "10"
meta.asv$Week[grep("week14", meta.asv$Sample_ID)] = "14"
meta.asv$Week[grep("week18", meta.asv$Sample_ID)] = "18"
meta.asv$Week[grep("week21", meta.asv$Sample_ID)] = "21"
meta.asv$Week[grep("week22", meta.asv$Sample_ID)] = "22"
meta.asv$Week[grep("week24", meta.asv$Sample_ID)] = "24"
meta.asv$Week[grep("week28", meta.asv$Sample_ID)] = "28"
meta.asv$Week[grep("week30", meta.asv$Sample_ID)] = "30"

# Controls Sex as Inocula/Controls
meta.asv$`Danska Mice::Sex`[grep("Negative|Positive|defined|inoculum|Pooled|control|ctrl", meta.asv$Sample_ID)] = "Inocula/Controls"
meta.asv$`Danska Mice::Sex` = tidyr::replace_na(meta.asv$`Danska Mice::Sex`, "Controls")

meta.asv$`Danska Mice::Age At Retirement` = as.numeric(meta.asv$`Danska Mice::Age At Retirement`)

meta.asv$`Danska Mice::Fate Primary`[grep("NS1", meta.asv$Sample_ID)] = "NS1"
meta.asv$`Danska Mice::Fate Primary`[grep("S2", meta.asv$Sample_ID)] = "S2"
meta.asv$`Danska Mice::Fate Primary`[grep("S5", meta.asv$Sample_ID)] = "S5"
meta.asv$`Danska Mice::Fate Primary`[grep("NS6", meta.asv$Sample_ID)] = "NS6"

# Changes Fate Primary of pos/neg controls to controls (did th mice recieve NS1/S2 community or control)
meta.asv$`Danska Mice::Fate Primary` = tidyr::replace_na(meta.asv$`Danska Mice::Fate Primary`, "Controls")

meta.asv$Sample_ID = gsub("postive", "positive", meta.asv$Sample_ID)

meta.asv$`Danska Mice::Fate Primary` = factor(meta.asv$`Danska Mice::Fate Primary`, levels = c("NS1", "NS6", "S2", "S5","Controls"))

meta.asv$`Danska Mice::Fate Primary` = factor(meta.asv$`Danska Mice::Diabetic?`, levels = c("yes", "no"))

meta.asv$Sample_ID = factor(meta.asv$Sample_ID, levels = c(grep("defined|_inocu|ctrl", meta.asv$Sample_ID, value = T),
                                                           #grep("_inocu", meta.asv$Sample_ID, value = T),
                                                           #grep("ctrl", meta.asv$Sample_ID, value = T),
                                                           grep("(^\\d.*)+(NS1)", meta.asv$Sample_ID, value = T),
                                                           grep("(^\\d.*)+(NS6)", meta.asv$Sample_ID, value = T),
                                                           grep("(^\\d.*)+(S2)", meta.asv$Sample_ID, value = T),
                                                           grep("(^\\d.*)+(S5)", meta.asv$Sample_ID, value = T),
                                                           grep("control|extraction|Pooled", meta.asv$Sample_ID, value = T)))

## Setup OTU table from absolute reads -----------------------------------------
otu = asv[,!grepl("ASV|BLAST|sequences|per_id|bit|evalue", names(asv))]

## Setup taxa ID table ---------------------------------------------------------
tax = asv[,grepl("BLAST", names(asv))]

## Set up metadata table -------------------------------------------------------
meta.asv <- meta.asv %>%
  tibble::column_to_rownames("Sample")

meta.asv$Plate = as.factor(meta.asv$Plate)

# Set up all OTU and taxonomic tables as matrices, metadata table can be left as is
otu_mat <- as.matrix(otu)
tax_mat <- as.matrix(tax)

# Create PS Object -------------------------------------------------------------
# Import OTU, taxonomic matrices and metadata table into Phyloseq object file 
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
meta = sample_data(meta.asv)
ps <- phyloseq(OTU, TAX, meta)

# Normalization
# This will normalize reads based on relative abundance across all samples 
normalizer = function(x) (100 * (x / sum(x)))
ps.norm = transform_sample_counts(ps, normalizer)

