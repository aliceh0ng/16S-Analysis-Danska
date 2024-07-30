# 
library(DESeq2)
library(dplyr)
library(tidyverse)
library(ALDEx2)
library(edgeR)

# S2 late vs. early comparison with continuous variable instead of binary ------

## From deseq.R script ->
# samples3
# counData3
# contConds
colData.cont <- meta5 %>% select('Danska Mice::Age At Retirement') %>% rownames_to_column(var = "Sample_ID")
contConds <- colData.cont$`Danska Mice::Age At Retirement`

# Create a DGEList object
d3 <- DGEList(counts=countData3)

# Normalize the data
d3 <- calcNormFactors(d3)

# Define design matrix
design <- model.matrix(~ contConds)

# Estimate dispersion
d3 <- estimateDisp(d3, design)

# Fit the model
fit3 <- glmFit(d3, design)

# Perform likelihood ratio test
lrt3 <- glmLRT(fit3, coef=2)

results <- topTags(lrt3, n = Inf)

View(lrt3)
topTags(lrt3)

# Plot the results (optional)
plotMD(lrt3, column = 2)
abline(h = 0, col = "red", lty = 2, lwd = 2)

# S2 late vs. early, inclusive!, comparison with continuous variable instead of binary ------

## From deseq.R script ->
# meta4
# samples2
# counData2
# contConds2

# colData
meta4 <- meta.asv 

meta4[ , 'DiabetesType'] = NA

meta4 <- meta4 %>% filter(grepl('30 weeks', Experiment_Length)) %>% 
  filter(grepl('S2', `Danska Mice::Cohort Code`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 156, "early", "late"))

# should i do just week 9 or all weeks
meta4 <- meta4 %>% filter(grepl('week9', `Sample_ID`)) 

colData.cont2 <- meta4 %>% select('Danska Mice::Age At Retirement') %>% rownames_to_column(var = "Sample_ID")
contConds2 <- colData.cont2$`Danska Mice::Age At Retirement`

# countData
samples2 <- rownames(meta4)
countData2 <- as.matrix(otu.S2 %>% select(all_of(samples2)))

# Create a DGEList object
d2 <- DGEList(counts=countData2)

# Normalize the data
d2 <- calcNormFactors(d2)

# Define design matrix
design2 <- model.matrix(~ contConds2)

# Estimate dispersion
d2 <- estimateDisp(d2, design2)

# Fit the model
fit2 <- glmFit(d2, design2)

# Perform likelihood ratio test
lrt2 <- glmLRT(fit2, coef=2)

results <- topTags(lrt2, n = Inf)

View(lrt2)
topTags(lrt2)


-------
  
# NS1 inclusive!

# meta3
# samples3
# countData
  
samples <- rownames(meta3)
countData <- as.matrix(otu.NS1 %>% select(all_of(samples)))
  
colData.cont3 <- meta3 %>% select('Danska Mice::Age At Retirement') %>% rownames_to_column(var = "Sample_ID")
contConds3 <- colData.cont3$`Danska Mice::Age At Retirement`

# Create a DGEList object
d4 <- DGEList(counts=countData)

# Normalize the data
d4 <- calcNormFactors(d4)

# Define design matrix
design4 <- model.matrix(~ contConds3)

# Estimate dispersion
d4 <- estimateDisp(d4, design4)

# Fit the model
fit4 <- glmFit(d4, design4)

# Perform likelihood ratio test
lrt4 <- glmLRT(fit4, coef=2)

results <- topTags(lrt4, n = Inf)

View(lrt4)
topTags(lrt4)

# Plot the results (optional)
plotMD(lrt4, column = 2)
abline(h = 0, col = "red", lty = 2, lwd = 2)

