# Differential abundance analysis taking the common bugs from the NS1 and S2 communities

# All 30 week samples for yes/no ---------

# Comment out this step to exclude Bl (longum vs. longum subsp. suillum)
new.asv.S2.Bl <- new.asv.S2 %>% 
  mutate(`BLAST_species.(as.done.by.EAV)` = if_else(`BLAST_species.(as.done.by.EAV)` == "Bifidobacterium longum subsp. suillum", "Bifidobacterium longum", `BLAST_species.(as.done.by.EAV)`))

# common.NS1.S2.collapsed <- inner_join(new.asv.NS1, new.asv.S2, by = "BLAST_species.(as.done.by.EAV)")
common.NS1.S2.collapsed <- inner_join(new.asv.NS1, new.asv.S2.Bl, by = "BLAST_species.(as.done.by.EAV)") %>%
  select(ends_with(".x"), ends_with(".y"), everything()) %>%
  select(-c("ref_num.x", "ref_num.y", "mouse_inocula.x", "mouse_inocula.y"))

# New meta data with just NS1 and S2 samples

meta.NS1.S2 <- meta.asv %>% filter(grepl('30 weeks', Experiment_Length)) %>%
  filter(grepl('NS1|S2', `MouseID`)) %>%
  filter(!str_detect(MouseID, "3971R_RL_F_NS1")) #%>% filter(grepl('week9', `CollectionTime`)) # For only week 5

samples <- rownames(meta.NS1.S2)
colData <- meta.NS1.S2 %>% select('Diabetes') %>% rownames_to_column(var = "Sample_ID")
conds <- colData$Diabetes
colData.cont <- meta.NS1.S2 %>% select('Danska Mice::Age At Retirement') %>% rownames_to_column(var = "Sample_ID")
conds.cont <- colData.cont$`Danska Mice::Age At Retirement`

# Visualisation of age at retirement & yes/no diabetes
p <- meta.NS1.S2 %>%
  ggplot(aes(x=`Danska Mice::Age At Retirement`, fill=`Danska Mice::Diabetic?`)) + geom_histogram(bins=12) + 
  #facet_wrap(~`Danska Mice::Fate Primary`) +
  scale_y_continuous(n.breaks=7)+
  theme_bw()

# count data from new ASV table
otu.NS1.S2 <- common.NS1.S2.collapsed %>% select(contains("sum")) %>% rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))
countData <- as.matrix(otu.NS1.S2 %>% select(all_of(samples)))

# edgeR binary
d <- DGEList(counts=countData, group=conds) # Create a DGEList object
d <- calcNormFactors(d) # Normalize the data
d <- estimateDisp(d) # Estimate dispersion
fit <- glmFit(d, design=model.matrix(~ conds)) # Fit the model
lrt <- glmLRT(fit, coef=2) # Conduct the test
View(lrt)
topTags(lrt)
# Plot the results
plotMD(lrt, column = 2)
abline(h = 0, col = "red", lty = 2, lwd = 2)

# edgeR cont
d <- DGEList(counts=countData) # Create a DGEList object
d <- calcNormFactors(d) # Normalize the data
design <- model.matrix(~ conds.cont) # Define design matrix
d <- estimateDisp(d, design) # Estimate dispersion
fit <- glmFit(d, design) # Fit the model
lrt.cont <- glmLRT(fit, coef=2) # Perform likelihood ratio test
results.cont <- topTags(lrt.cont, n = Inf)
View(lrt.cont)
topTags(lrt.cont)

## Visualization -----

# create new full joined ASV table for NS1 and S2 samples to calculate actual rel. abundance
NS1.S2.full.asv <- full_join(new.asv.NS1, new.asv.S2.Bl, by = "BLAST_species.(as.done.by.EAV)") %>%
  select(ends_with(".x"), ends_with(".y"), everything()) %>%
  select(-c("ref_num.x", "ref_num.y", "mouse_inocula.x", "mouse_inocula.y", "16S.Seq.result.(JD).x", "16S.Seq.result.(JD).y")) %>%
  mutate_if(is.numeric,coalesce,0)

# relative abundance table for all samples
NS1.S2.rel <- NS1.S2.full.asv 
NS1.S2.rel[, 8:ncol(NS1.S2.rel)] <- sapply(NS1.S2.rel[, 8:ncol(NS1.S2.rel)],function(x) x/sum(x) * 100)

NS1.S2.rel <- NS1.S2.rel %>% column_to_rownames("BLAST_species.(as.done.by.EAV)") %>% select(starts_with("sum"))

diabetes <- meta.NS1.S2 %>% select('Diabetes') %>% rownames_to_column(var = "Sample")

# Collinsella aerofaciens
#c.a. <- data.frame(Sample = colnames(NS1.S2.rel %>% select(starts_with("sum"))), Relative_Abundance = as.numeric(NS1.S2.rel[19,] %>% select(starts_with("sum"))))
c.a.<- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Collinsella aerofaciens",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  #mutate(Week = gsub(".*week", "", Sample))
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

c.a. <- inner_join(c.a., diabetes, by = "Sample")

plot <- ggplot(c.a., aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  facet_wrap(~Community) +
  theme_bw() 

# Eggerthella lenta
e.l. <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Eggerthella lenta",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

e.l. <- inner_join(e.l., diabetes, by = "Sample")

plot <- ggplot(e.l., aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  facet_wrap(~Community) +
  theme_bw() 

# Parabacteroides distasonis
p.d. <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Parabacteroides distasonis",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

p.d. <- inner_join(p.d., diabetes, by = "Sample")

plot <- ggplot(p.d., aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  facet_wrap(~Community) +
  theme_bw() 

# Alistipes finegoldii
a.f. <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Alistipes finegoldii",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

a.f. <- inner_join(a.f., diabetes, by = "Sample")

plot <- ggplot(a.f., aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  facet_wrap(~Community) +
  theme_bw() 
  
## Rest of common 12 species boxplots ------------

# Akkermansia muciniphila
am <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Akkermansia muciniphila",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

am <- inner_join(am, diabetes, by = "Sample")

p1 <- ggplot(am, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Akkermansia muciniphila")

# Clostridium bolteae
cb <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Clostridium bolteae",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

cb <- inner_join(cb, diabetes, by = "Sample")

p2 <- ggplot(cb, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Clostridium bolteae")

# Eisenbergiella tayi
et <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Eisenbergiella tayi",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

et <- inner_join(et, diabetes, by = "Sample")

p3 <- ggplot(et, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Eisenbergiella tayi")

# Flavonifractor plautii
fp <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Flavonifractor plautii",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

fp <- inner_join(fp, diabetes, by = "Sample")

p4 <- ggplot(fp, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Flavonifractor plautii")

# Hungatella effluvii
he <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Hungatella effluvii",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

he <- inner_join(he, diabetes, by = "Sample")

p5 <- ggplot(he, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Hungatella effluvii")

# Staphylococcus hominis subsp. novobiosepticus
SHN <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Staphylococcus hominis subsp. novobiosepticus",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

SHN <- inner_join(SHN, diabetes, by = "Sample")

p6 <- ggplot(SHN, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Staphylococcus hominis subsp. novobiosepticus")

# Bacteroides vulgatus
bv <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Bacteroides vulgatus",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

bv <- inner_join(bv, diabetes, by = "Sample")

p7 <- ggplot(bv, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Staphylococcus hominis subsp. novobiosepticus")

# Bifidobacterium longum
bl <- data.frame(Sample = colnames(NS1.S2.rel), Relative_Abundance = as.numeric(NS1.S2.rel["Bifidobacterium longum",])) %>%
  mutate(Sample = gsub("^sum_", "", Sample)) %>%
  mutate(Community = ifelse(grepl("NS1", Sample), "NS1", ifelse(grepl("S2", Sample), "S2", NA))) %>%
  mutate(Week = as.numeric(gsub(".*week", "", Sample)))

bl <- inner_join(bl, diabetes, by = "Sample")

p8 <- ggplot(bl, aes(x = factor(Week), y = Relative_Abundance)) +
  geom_boxplot(aes(color = Diabetes, fill = Diabetes)) +
  scale_fill_manual(values = c("yes" = "pink", "no" = "lightblue")) + # Custom colors for Diabetes status
  scale_color_manual(values = c("yes" = "red", "no" = "navy")) + # Custom colors for Diabetes status
  #scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  #facet_wrap(~Community) +
  theme_bw() +
  ggtitle("Staphylococcus hominis subsp. novobiosepticus")








# Just yes diabetes 30 week samples for "early"/"late" -------------------------

# New meta data table with just yes 30 week samples from NS1 and S2 mice

# Median for age at retirement is 146 (149.5 if including all week samples, but 146 if only 1 sample from each mouse)

meta.NS1.S2.yes <- meta.asv %>% filter(grepl('30 weeks', Experiment_Length)) %>%
  filter(grepl('NS1|S2', `MouseID`)) %>%
  filter(grepl('yes', `Danska Mice::Diabetic?`)) %>%
  mutate(DiabetesType = ifelse(`Danska Mice::Age At Retirement` < 147, "early", "late")) %>%
  filter(!str_detect(MouseID, "3971R_RL_F_NS1")) # %>% filter(grepl('week5', `CollectionTime`)) # For only week 5

samples <- rownames(meta.NS1.S2.yes)
colData <- meta.NS1.S2.yes %>% select('DiabetesType') %>% rownames_to_column(var = "Sample_ID")
conds <- colData$DiabetesType
colData.cont <- meta.NS1.S2.yes %>% select('Danska Mice::Age At Retirement') %>% rownames_to_column(var = "Sample_ID")
conds.cont <- colData.cont$`Danska Mice::Age At Retirement`

# Visualisation of age at retirement
p2 <- meta.NS1.S2.yes %>%
  ggplot(aes(x=`Danska Mice::Age At Retirement`)) + geom_histogram(bins=9) + 
  #facet_wrap(~`Danska Mice::Fate Primary`) +
  scale_y_continuous(n.breaks=7)+
  geom_vline(xintercept = 147, linetype = "dashed", color = "red") +
  theme_bw()

# count data from new ASV table
otu.NS1.S2 <- common.NS1.S2.collapsed %>% select(contains("sum")) %>% rename_with(~ str_remove(., "^sum_"), starts_with("sum_"))
countData <- as.matrix(otu.NS1.S2 %>% select(all_of(samples)))

# edgeR binary
d <- DGEList(counts=countData, group=conds) # Create a DGEList object
d <- calcNormFactors(d) # Normalize the data
d <- estimateDisp(d) # Estimate dispersion
fit <- glmFit(d, design=model.matrix(~ conds)) # Fit the model
lrt <- glmLRT(fit, coef=2) # Conduct the test
topTags(lrt)
# Plot the results
plotMD(lrt, column = 2)
abline(h = 0, col = "red", lty = 2, lwd = 2)

# edgeR cont
d <- DGEList(counts=countData) # Create a DGEList object
d <- calcNormFactors(d) # Normalize the data
design <- model.matrix(~ conds.cont) # Define design matrix
d <- estimateDisp(d, design) # Estimate dispersion
fit <- glmFit(d, design) # Fit the model
lrt.cont <- glmLRT(fit, coef=2) # Perform likelihood ratio test
results.cont <- topTags(lrt.cont, n = Inf)
topTags(lrt.cont)

