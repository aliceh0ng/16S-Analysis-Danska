# 
library(DESeq2)
library(dplyr)
library(tidyverse)
library(ALDEx2)
library(edgeR)


# Preparing data for comparison 2 NS1 yes or no --------------------------------

# countData
samples <- rownames(meta3)
countData <- as.matrix(otu.NS1 %>% select(all_of(samples)))
#countData <- as.matrix(otu.NS1 %>% select(all_of(samples)) %>% select(contains("week9") | contains("week5")))

# colData
meta3 <- meta.asv 

# To do 1 timepoint only for each mouse?
# meta3 <- meta3 %>% distinct(MouseID, .keep_all = TRUE)

meta3 <- meta3 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('NS1', `MouseID`)) %>%
  filter(!str_detect(MouseID, "3971R_RL_F_NS1"))

# should i do just week 9 or all weeks
# meta3 <- meta3 %>% filter(grepl('week9', `Sample_ID`)) 

#meta3$Diabetes = factor(meta3$Diabetes, levels = c("yes", "no"))
#meta2$`Danska Mice::Cohort Code` = factor(meta2$`Danska Mice::Cohort Code`, levels = c("GermFree_DIABIMMUNE_S2_T1Dendpoint", "GermFree_DIABIMMUNE_S5_T1Dendpoint"))

colData <- meta3 %>% select('Diabetes') %>% rownames_to_column(var = "Sample_ID")

meta_3 = sample_data(meta3)

samples <- rownames(meta_3)
filtered.OTU.NS1 <- otu_table(as.matrix(otu.NS1 %>% select(all_of(samples)) %>% select(contains("week9") | contains("week5"))), taxa_are_rows = TRUE) 

ps3 <- phyloseq(filtered.OTU.NS1, tax.NS1, meta_3)


# aldex (Comparison 2 test) -------
x.aldex <- aldex(countData, conds, mc.samples = 16)


# deseq (Comparison 2 test) ------

# Assuming you have a count matrix 'counts' and a metadata dataframe 'colData'
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Diabetes)
dds <- DESeq(dds)
results <- results(dds)


# edgeR NS1: Comparison 2 (NS1 yes or no) --------------------------------------

# Sample metadata
conds <- colData$Diabetes

# Create a DGEList object
d <- DGEList(counts=countData, group=conds)

# Normalize the data
d <- calcNormFactors(d)

# Estimate dispersion
d <- estimateDisp(d)

# Fit the model
fit <- glmFit(d, design=model.matrix(~ conds))

# Conduct the test
lrt <- glmLRT(fit, coef=2)
View(lrt)

topTags(lrt)


# edgeR S2: Comparison 1b Include 'no' from the S Community --------------------

# colData
meta4 <- meta.asv 

meta4[ , 'DiabetesType'] = NA

meta4 <- meta4 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('S2', `Danska Mice::Cohort Code`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 156, "early", "late"))

# should i do just week 9 or all weeks
# meta4 <- meta4 %>% filter(grepl('week9', `Sample_ID`)) 

#meta4$DiabetesType = factor(meta4$DiabetesType, levels = c("early", "late"))

colData2 <- meta4 %>% select('DiabetesType') %>% rownames_to_column(var = "Sample_ID")

# countData
samples2 <- rownames(meta4)
countData2 <- as.matrix(otu.S2 %>% select(all_of(samples2)))

conds2 <- colData2$DiabetesType

# Create a DGEList object
d2 <- DGEList(counts=countData2, group=conds2)

# Normalize the data
d2 <- calcNormFactors(d2)

# Estimate dispersion
d2 <- estimateDisp(d2)

# Fit the model
fit2 <- glmFit(d2, design=model.matrix(~ conds2))

# Conduct the test
lrt2 <- glmLRT(fit2, coef=2)
View(lrt2)

topTags(lrt2)



# Sample count matrix
counts <- matrix(sample(0:100, 2000, replace = TRUE), nrow = 20, ncol = 100)
rownames(counts) <- paste0("Taxon", 1:20)
colnames(counts) <- paste0("Sample", 1:100)

# Load necessary package
install.packages("ggplot2")
library(ggplot2)

# Convert to a long format dataframe for ggplot2
countData_df <- as.data.frame(countData)
countData_long <- reshape2::melt(countData_df)

# Histogram of all counts
ggplot(countData_long, aes(value)) +
  geom_histogram(binwidth = 1000, fill = "blue", color = "black") +
  labs(title = "Histogram of Count Data", x = "Counts", y = "Frequency") +
  theme_minimal()


# EdgeR S2 community exclusive -------------------------------------------------

# meta2, 141
meta5 <- meta.asv 

meta5[ , 'DiabetesType'] = NA

meta5 <- meta5 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('S2', `Danska Mice::Cohort Code`)) %>%
  filter(grepl('yes', `Danska Mice::Diabetic?`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 141, "early", "late"))

# should i do just week 9 or all weeks
# meta4 <- meta4 %>% filter(grepl('week9', `Sample_ID`)) 

#meta4$DiabetesType = factor(meta4$DiabetesType, levels = c("early", "late"))

colData3 <- meta5 %>% select('DiabetesType') %>% rownames_to_column(var = "Sample_ID")

# countData
samples3 <- rownames(meta5)
countData3 <- as.matrix(otu.S2 %>% select(all_of(samples3)))

conds3 <- colData3$DiabetesType

# Create a DGEList object
d3 <- DGEList(counts=countData3, group=conds3)

# Normalize the data
d3 <- calcNormFactors(d3)

# Estimate dispersion
d3 <- estimateDisp(d3)

# Fit the model
fit3 <- glmFit(d3, design=model.matrix(~ conds3))

# Conduct the test
lrt3 <- glmLRT(fit3, coef=2)
View(lrt3)

topTags(lrt3)
