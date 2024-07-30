## Script to create ps object from new (combined) asv table

library(readxl)
library(phyloseq)
library(dplyr)

# PS Object for DA analysis ----------------------------------------------------

## Set up ----------------------------------------------------------------------

# Need to run first:
# - load seqtab.nochim.plateX
# - run 'set-up' 
# - run 'create-ps'
# - run "MSA"

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

 # meta.asv <- meta.asv %>%
 #   tibble::column_to_rownames("Sample_ID")





## Create OTU and TAX tables from new asv table --------------------------------

# OTU tables - first create arbitrary otu numbers 
otu.NS1 <- new.asv.NS1[,grepl("sum", names(new.asv.NS1))] %>%
  rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))
otu.NS6 <- new.asv.NS6[,grepl("sum", names(new.asv.NS6))] %>%
  rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))
otu.S2 <- new.asv.S2[,grepl("sum", names(new.asv.S2))] %>%
  rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))
otu.S5 <- new.asv.S5[,grepl("sum", names(new.asv.S5))] %>%
  rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))

tax.NS1 = new.asv.NS1[,grepl("Sample_ID|isolate_ID|BLAST|16S", names(new.asv.NS1))]
tax.NS6 = new.asv.NS6[,grepl("Sample_ID|isolate_ID|BLAST|16S", names(new.asv.NS6))]
tax.S2 = new.asv.S2[,grepl("Sample_ID|isolate_ID|BLAST|16S", names(new.asv.S2))]
tax.S5 = new.asv.S5[,grepl("Sample_ID|isolate_ID|BLAST|16S", names(new.asv.S5))]

# Turn into Phyloseq object file 
OTU.NS1 = otu_table(as.matrix(otu.NS1), taxa_are_rows = TRUE)
OTU.NS6 = otu_table(as.matrix(otu.NS6), taxa_are_rows = TRUE)
OTU.S2 = otu_table(as.matrix(otu.S2), taxa_are_rows = TRUE)
OTU.S5 = otu_table(as.matrix(otu.S5), taxa_are_rows = TRUE)
TAX.NS1 = tax_table(as.matrix(tax.NS1))
TAX.NS6 = tax_table(as.matrix(tax.NS6))
TAX.S2 = tax_table(as.matrix(tax.S2))
TAX.S5 = tax_table(as.matrix(tax.S5))


## Set up metadata table -------------------------------------------------------

### Comparison 1b Include 'no' from the S Community ----------------------------

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

samples <- rownames(meta_4)
filtered.OTU.S2 <- otu_table(as.matrix(otu.S2 %>% select(all_of(samples))), taxa_are_rows = TRUE)

ps4 <- phyloseq(filtered.OTU.S2, TAX.S2, meta_4)

# Comparison 2 NS1 yes or no ---------------------------------------------------

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

samples <- rownames(meta_3)
filtered.OTU.NS1 <- otu_table(as.matrix(otu.NS1 %>% select(all_of(samples)) %>% select(contains("week9") | contains("week5"))), taxa_are_rows = TRUE) 

ps3 <- phyloseq(filtered.OTU.NS1, tax.NS1, meta_3)
