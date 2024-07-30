library(phyloseq)

# PS Object for DA analysis ----------------------------------------------------

## Set up ----------------------------------------------------------------------

### Create/organize meta data table --------------------------------------------
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

### Organize meta.asv for ancom ps object --------------------------------------

meta.asv <- meta.asv %>% 
  select("Sample_ID", "CohortCode_MouseSample", "CollectionDateActual", "CollectionTime", "MouseID", "MouseID_Count", "SampleType", "Danska Mice::Age At Retirement", "Danska Mice::Cohort Code", "Danska Mice::Diabetic?", "Danska Mice::Experiment_ID", "Danska Mice::Fate Primary", "Danska Mice::Fate Secondary", "Danska Mice::insulitis score?", "Danska Mice::Sex", "Danska Mice::Spontaneous Insulitis Overall Score", "Danska Mice::Strain") %>%
  filter(!grepl('ABX', `MouseID`)) %>% # remove ABX mice
  filter(!grepl('Control', `Danska Mice::Sex`)) %>% # remove controls (cage pellets, defined inocula) 
  mutate(Experiment_Length = case_when(grepl("T1D", `Danska Mice::Cohort Code`) ~ "30 weeks",
                                       grepl("10week", `Danska Mice::Cohort Code`) ~ "10 weeks",
                                       TRUE ~ `Danska Mice::Cohort Code`)) %>%
  mutate(Diabetes = case_when(grepl("yes", `Danska Mice::Diabetic?`) ~ "yes",
                              grepl("no", `Danska Mice::Diabetic?`) ~ "no",
                              TRUE ~ `Danska Mice::Diabetic?`)) 

normalizer = function(x) (100 * (x / sum(x)))


## Comparison 1: Yes Diabetes, S2 only -----------------------------------------

meta2 <- meta.asv

# To do 1 timepoint only for each mouse?
#meta2 <- meta2 %>% filter(grepl('week9', CollectionTime))

meta2[ , 'DiabetesType'] = NA

meta2 <- meta2 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('yes', `Danska Mice::Diabetic?`)) %>%
  filter(grepl('S2', `Danska Mice::Cohort Code`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 141, "early", "late"))

# should i do just week 9 or all weeks
# meta2 <- meta2 %>% filter(grepl('week9', `Sample_ID`)) 

meta2$DiabetesType = factor(meta2$DiabetesType, levels = c("early", "late"))
#meta2$`Danska Mice::Cohort Code` = factor(meta2$`Danska Mice::Cohort Code`, levels = c("GermFree_DIABIMMUNE_S2_T1Dendpoint", "GermFree_DIABIMMUNE_S5_T1Dendpoint"))

meta_2 = sample_data(meta2)

ps2 <- phyloseq(OTU, TAX, meta_2)

ps2.norm = transform_sample_counts(ps2, normalizer)

test <- rowSums(ps2@otu_table)>20 #1? 20?
OTU.test <- OTU[test, ]
TAX.test <- TAX[test, ]
ps2.test <- phyloseq(OTU.test, TAX.test, meta_2)
ps2.test.norm <- transform_sample_counts(ps2.test, normalizer)

### testing

cols_to_keep <- rownames(meta_2)
otu_filtered <- otu %>% select(all_of(cols_to_keep))
otu_mat_filtered <- as.matrix(otu_filtered)

OTU_filtered = otu_table(otu_mat_filtered, taxa_are_rows = TRUE)
ps2.filtered <- phyloseq(OTU_filtered, TAX, meta_2)

ps2.filtered.norm = transform_sample_counts(ps2.filtered, normalizer)

g1_colsum <- colSums(otu_filtered)

####

# # central log transformation -> normalizes ther read to a standard normal dist. based on each bug (mean abundance is 0, devided by)
# library(microViz)
# ps2.cl <- tax_transform(ps2, "clr")


## Comparison 1b Include 'no' from the S Community --------------------

meta4 <- meta.asv 

# To do 1 timepoint only for each mouse?
#meta4 <- meta4 %>% distinct(MouseID, .keep_all = TRUE)

meta4[ , 'DiabetesType'] = NA

meta4 <- meta4 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  #filter(grepl('yes', `Danska Mice::Diabetic?`)) %>%
  filter(grepl('S2', `Danska Mice::Cohort Code`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 156, "early", "late"))

# should i do just week 9 or all weeks
# meta4 <- meta4 %>% filter(grepl('week9', `Sample_ID`)) 

meta4$DiabetesType = factor(meta4$DiabetesType, levels = c("early", "late"))
#meta2$`Danska Mice::Cohort Code` = factor(meta2$`Danska Mice::Cohort Code`, levels = c("GermFree_DIABIMMUNE_S2_T1Dendpoint", "GermFree_DIABIMMUNE_S5_T1Dendpoint"))

meta_4 = sample_data(meta4)

ps4 <- phyloseq(OTU, TAX, meta_4)

ps4.norm = transform_sample_counts(ps4, normalizer)


## Comparison 2: NS1 yes vs no -------------------------------------------------

meta3 <- meta.asv 

# To do 1 timepoint only for each mouse?
#meta3 <- meta3 %>% distinct(MouseID, .keep_all = TRUE)

meta3 <- meta3 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('NS1', `MouseID`)) %>%
  filter(!str_detect(MouseID, "3971R_RL_F_NS1"))

# should i do just week 9 or all weeks
# meta3 <- meta3 %>% filter(grepl('week9', `Sample_ID`)) 

meta3$Diabetes = factor(meta3$Diabetes, levels = c("yes", "no"))
#meta2$`Danska Mice::Cohort Code` = factor(meta2$`Danska Mice::Cohort Code`, levels = c("GermFree_DIABIMMUNE_S2_T1Dendpoint", "GermFree_DIABIMMUNE_S5_T1Dendpoint"))

meta_3 = sample_data(meta3)

ps3 <- phyloseq(OTU, TAX, meta_3)

ps3.norm = transform_sample_counts(ps3, normalizer)

test <- rowSums(ps3@otu_table)>20 #1? 20?
OTU.test <- OTU[test, ]
TAX.test <- TAX[test, ]
ps3.test <- phyloseq(OTU.test, TAX.test, meta_3)
ps3.test.norm <- transform_sample_counts(ps3.test, normalizer)


### testing

cols_to_keep3 <- rownames(meta_3)
otu_filtered3 <- otu %>% select(all_of(cols_to_keep3))
otu_mat_filtered3 <- as.matrix(otu_filtered3)

OTU_filtered3 = otu_table(otu_mat_filtered3, taxa_are_rows = TRUE)
ps3.filtered <- phyloseq(OTU_filtered3, TAX, meta_3)

ps3.filtered.norm = transform_sample_counts(ps3.filtered, normalizer)

# colsums
g2_colsum <- colSums(otu_filtered3)

## -------------------

cols_to_keep4 <- rownames(meta_4)
otu_filtered4 <- otu %>% select(all_of(cols_to_keep4))
otu_mat_filtered4 <- as.matrix(otu_filtered4)

OTU_filtered4 = otu_table(otu_mat_filtered4, taxa_are_rows = TRUE)
ps4.filtered <- phyloseq(OTU_filtered4, TAX, meta_4)

ps4.filtered.norm = transform_sample_counts(ps3.filtered, normalizer)
