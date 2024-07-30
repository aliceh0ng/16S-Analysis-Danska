# From 'MSA':
# new.asv.NS1 <- newASV(WGS_community[[1]], asv.list[["NS1"]], top_hits_NS1)
# new.asv.NS6 <- newASV(WGS_community[[2]], asv.list[["NS6"]], top_hits_NS6)
# new.asv.S2 <- newASV(WGS_community[[3]], asv.list[["S2"]], top_hits_S2)
# new.asv.S5 <- newASV(WGS_community[[4]], asv.list[["S5"]], top_hits_S5)

# meta 3 -> samples included in the NS1 yes/no diabetes comparison -------------

LLA_NS1_yes <- meta3 %>% select(Diabetes) %>% rownames_to_column('Sample') %>% filter(Diabetes == 'yes')
samples_yes <- LLA_NS1_yes$Sample
LLA_NS1_no <- meta3 %>% select(Diabetes) %>% rownames_to_column('Sample') %>% filter(Diabetes == 'no')
samples_no <- LLA_NS1_no$Sample

LLA_NS1_yes <- cbind(tax.NS1, otu.NS1 %>% select(all_of(samples_yes)))
LLA_NS1_no <- cbind(tax.NS1, otu.NS1 %>% select(all_of(samples_no)))

#LLA_NS1_temp <- LLA_NS1
LLA_NS1_yes_weeks <- tax.NS1
LLA_NS1_no_weeks <- tax.NS1

# Obtain the sum of all samples in the same treatment and week and assign to new asv table
  
# All weeks present across all plates
weeks <- c("5","6","7","8","9","10")
weeks <- c("5","6","7","9")

for(j in weeks){
  
  # Move samples in the same community and week together
  asv.temp <- LLA_NS1_yes %>%
    dplyr::relocate(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)), .after = "16S.Seq.result.(JD)")
  
  # Create a new combined column of all the samples within the community/week group
  if(length(grep(paste0(sep = ".*", sep="week", j), names(asv.temp))) > 0) {
    row_sum <- rowSums(asv.temp[,(first(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))):(last(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))), drop = F])
    
    LLA_NS1_yes_weeks[[paste0(sep= "week", j)]] <- row_sum
  }
}

# Relative abundance table
LLA_NS1_yes_weeks_norm = LLA_NS1_yes_weeks %>% cbind(otu.NS1 %>% select(contains("inocu"))) %>%
  mutate(NS1_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_NS1_yes_weeks_norm[, 5:9] <- sapply(LLA_NS1_yes_weeks_norm[, 5:9],function(x) x/sum(x) * 100)

LLA_NS1_no_weeks_norm = LLA_NS1_no_weeks %>% cbind(otu.NS1 %>% select(contains("inocu"))) %>%
  mutate(NS1_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_NS1_no_weeks_norm[, 5:9] <- sapply(LLA_NS1_no_weeks_norm[, 5:9],function(x) x/sum(x) * 100)

# Add ASV Number column
LLA_NS1_yes_weeks_norm <- LLA_NS1_yes_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_NS1_yes_weeks_norm$Sample_ID)) %>%
  relocate(NS1_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)
LLA_NS1_no_weeks_norm <- LLA_NS1_no_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_NS1_no_weeks_norm$Sample_ID)) %>%
  relocate(NS1_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)

## Plot Alluvials. Representivity % (rel. abundance) ---------------------------------

inocula.list <- list(NS1_yes = LLA_NS1_yes_weeks_norm, NS1_no = LLA_NS1_no_weeks_norm)

# Make melted dataframe for plotting
alluvial.plots = plyr::llply(inocula.list, function(x){
  y = x
  
  # Add 'score' column: if asv is present in the defined inocula, return TRUE, and if the asv is present in the sample, return TRUE (e.g. TRUE.TRUE)
  y$score = interaction(rowSums(x[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0,
                        rowSums(x[,(first(grep("week", names(x)))):(last(grep("week", names(x)))), drop = F])  > 0 
  )
  
  # Add 'in_any_inoc' column: if asv is present in any of the inocula then return true
  y$in_any_inoc = rowSums(y[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0 
  
  y = y[order(y$in_any_inoc),]  # Move ASVs not in inocula to the top
  # y = y[order(y$score, decreasing = FALSE),]  # Uncomment this to order by presence/absence
  
  print(dim(y))  # Check dimensions
  print(head(y))
  
  # Ordering has to be done to ASV because of ggalluvial
  y$ASV_Number = seq_along(y$ASV_Number)
  y$ASV_Number = sprintf("%04d", y$ASV_Number)
  
  #colnames(y) = gsub("^X", "", colnames(y))
  #colnames(y) = gsub("week", "", colnames(y))
  #return(y)
  
  # Melt    
  y = reshape2::melt(y, measure.vars = 6:ncol(x), variable.name = "Sample", value.name = "Relative_abundance", as.is = T)
  
  y = y %>%
    mutate(Housing = Sample)
  y$Housing = gsub("^\\d.*", "GF", y$Housing)
  y$Housing = gsub("^T1D.*", "Defined inocula", y$Housing)
  
  y = tidyr::extract (y, Sample, into = c("Community", "Endpoint_temp"), regex = "(NS1|NS6|S2|S5)+(.*)", remove = FALSE)
  y = tidyr::extract(y, Endpoint_temp, into = "Endpoint", regex = "(week.*|day.*)", remove = TRUE)
  y$Endpoint = gsub("day63", "week9", y$Endpoint)
  y$Endpoint = gsub("day70", "week10", y$Endpoint)
  
  y$Endpoint = replace_na(y$Endpoint, "Defined inocula")
  y = y[y$Relative_abundance!=0,]
  
  #Comment these lines out if don't want to group by plate
  # This line sorts the plate 1 and 2 samples so that they are 'mixed' rather than separated
  # y$Sample = factor(as.character(y$Sample), levels = sort(levels(y$Sample))) 
  
  # This relevels s.t. defined community is moved to front
  # y$Sample = relevel(y$Sample, grep("T1D", levels(y$Sample), value = T, perl = T))
  
  return(y)
}
)

# Name lists in alluvial.plots list if subsetted
# names(alluvial.plots) = c("NS1", "NS6", "S2", "S5") # commented because already exists

# Convert to df
alluvial.df = ldply(alluvial.plots)
alluvial.df.NS1.yes <- alluvial.plots[["NS1_yes"]]
alluvial.df.NS1.no <- alluvial.plots[["NS1_no"]]

# Setting custom colours for phyla plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = as.numeric(length(unique(alluvial.df$`BLAST_species.(as.done.by.EAV)`)))
cols = gg_color_hue(n)

names(cols) <- sort(unique(alluvial.df$`BLAST_species.(as.done.by.EAV)`))


p = ggplot() +
  geom_stratum(data = alluvial.df.NS1.yes, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df.NS1.yes, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  #facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 

p2 = ggplot() +
  geom_stratum(data = alluvial.df, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 


# meta 4 -> samples included in the S2 "late vs early" diabetes comparison -----

LLA_S2_early <- meta4 %>% select(DiabetesType) %>% rownames_to_column('Sample') %>% filter(DiabetesType == 'early')
samples_early <- LLA_S2_early$Sample

LLA_S2_late <- meta4 %>% select(DiabetesType) %>% rownames_to_column('Sample') %>% filter(DiabetesType == 'late')
samples_late <- LLA_S2_late$Sample

LLA_S2_early <- cbind(tax.S2, otu.S2 %>% select(all_of(samples_early)))
LLA_S2_late <- cbind(tax.S2, otu.S2 %>% select(all_of(samples_late)))

LLA_S2_early_weeks <- tax.S2
LLA_S2_late_weeks <- tax.S2

# Obtain the sum of all samples in the same treatment and week and assign to new asv table

# All weeks present across all plates
weeks <- c("5","6","7","8","9","10")
weeks <- c("5","6","7","9")

for(j in weeks){
  
  # Move samples in the same community and week together
  asv.temp <- LLA_S2_late %>%
    dplyr::relocate(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)), .after = "16S.Seq.result.(JD)")
  
  # Create a new combined column of all the samples within the community/week group
  if(length(grep(paste0(sep = ".*", sep="week", j), names(asv.temp))) > 0) {
    row_sum <- rowSums(asv.temp[,(first(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))):(last(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))), drop = F])
    
    LLA_S2_late_weeks[[paste0(sep= "week", j)]] <- row_sum
  }
}

# Relative abundance table
LLA_S2_early_weeks_norm = LLA_S2_early_weeks %>% cbind(otu.S2 %>% select(contains("inocu"))) %>%
  mutate(S2_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_S2_early_weeks_norm[, 5:9] <- sapply(LLA_S2_early_weeks_norm[, 5:9],function(x) x/sum(x) * 100)

LLA_S2_late_weeks_norm = LLA_S2_late_weeks %>% cbind(otu.S2 %>% select(contains("inocu"))) %>%
  mutate(S2_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_S2_late_weeks_norm[, 5:9] <- sapply(LLA_S2_late_weeks_norm[, 5:9],function(x) x/sum(x) * 100)

# Add ASV Number column
LLA_S2_early_weeks_norm <- LLA_S2_early_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_S2_early_weeks_norm$Sample_ID)) %>%
  relocate(S2_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)
LLA_S2_late_weeks_norm <- LLA_S2_late_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_S2_late_weeks_norm$Sample_ID)) %>%
  relocate(S2_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)

## Plot Alluvials. Representivity % (rel. abundance) ---------------------------------

inocula.list.2 <- list(S2_early = LLA_S2_early_weeks_norm, S2_late = LLA_S2_late_weeks_norm)

# Make melted dataframe for plotting
alluvial.plots.2 = plyr::llply(inocula.list.2, function(x){
  y = x
  
  # Add 'score' column: if asv is present in the defined inocula, return TRUE, and if the asv is present in the sample, return TRUE (e.g. TRUE.TRUE)
  y$score = interaction(rowSums(x[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0,
                        rowSums(x[,(first(grep("week", names(x)))):(last(grep("week", names(x)))), drop = F])  > 0 
  )
  
  # Add 'in_any_inoc' column: if asv is present in any of the inocula then return true
  y$in_any_inoc = rowSums(y[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0 
  
  y = y[order(y$in_any_inoc),]  # Move ASVs not in inocula to the top
  # y = y[order(y$score, decreasing = FALSE),]  # Uncomment this to order by presence/absence
  
  print(dim(y))  # Check dimensions
  print(head(y))
  
  # Ordering has to be done to ASV because of ggalluvial
  y$ASV_Number = seq_along(y$ASV_Number)
  y$ASV_Number = sprintf("%04d", y$ASV_Number)
  
  #colnames(y) = gsub("^X", "", colnames(y))
  #colnames(y) = gsub("week", "", colnames(y))
  #return(y)
  
  # Melt    
  y = reshape2::melt(y, measure.vars = 6:ncol(x), variable.name = "Sample", value.name = "Relative_abundance", as.is = T)
  
  y = y %>%
    mutate(Housing = Sample)
  y$Housing = gsub("^\\d.*", "GF", y$Housing)
  y$Housing = gsub("^T1D.*", "Defined inocula", y$Housing)
  
  y = tidyr::extract (y, Sample, into = c("Community", "Endpoint_temp"), regex = "(NS1|NS6|S2|S5)+(.*)", remove = FALSE)
  y = tidyr::extract(y, Endpoint_temp, into = "Endpoint", regex = "(week.*|day.*)", remove = TRUE)
  y$Endpoint = gsub("day63", "week9", y$Endpoint)
  y$Endpoint = gsub("day70", "week10", y$Endpoint)
  
  y$Endpoint = replace_na(y$Endpoint, "Defined inocula")
  y = y[y$Relative_abundance!=0,]
  
  #Comment these lines out if don't want to group by plate
  # This line sorts the plate 1 and 2 samples so that they are 'mixed' rather than separated
  # y$Sample = factor(as.character(y$Sample), levels = sort(levels(y$Sample))) 
  
  # This relevels s.t. defined community is moved to front
  # y$Sample = relevel(y$Sample, grep("T1D", levels(y$Sample), value = T, perl = T))
  
  return(y)
}
)

# Name lists in alluvial.plots list if subsetted
# names(alluvial.plots) = c("NS1", "NS6", "S2", "S5") # commented because already exists

# Convert to df
alluvial.df.2 = ldply(alluvial.plots.2)
alluvial.df.S2.early <- alluvial.plots.2[["S2_early"]]
alluvial.df.S2.late <- alluvial.plots.2[["S2_late"]]

# Setting custom colours for phyla plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = as.numeric(length(unique(alluvial.df.2$`BLAST_species.(as.done.by.EAV)`)))
cols = gg_color_hue(n)

names(cols) <- sort(unique(alluvial.df.2$`BLAST_species.(as.done.by.EAV)`))


p = ggplot() +
  geom_stratum(data = alluvial.df.NS1.yes, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df.NS1.yes, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  #facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 

p4 = ggplot() +
  geom_stratum(data = alluvial.df.2, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df.2, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 




# Common among NS1 and S2 ----------------
# Run diffA-common.R first (to generate meta.NS1.S2)

## NS1 and S2 Yes and No Diabetes -----------------

LLA_NS1S2_yes <- meta.NS1.S2 %>% select(Diabetes) %>% rownames_to_column('Sample') %>% filter(Diabetes == 'yes')
samples_yes <- LLA_NS1S2_yes$Sample
LLA_NS1S2_no <- meta.NS1.S2 %>% select(Diabetes) %>% rownames_to_column('Sample') %>% filter(Diabetes == 'no')
samples_no <- LLA_NS1S2_no$Sample

tax.NS1S2 <- common.NS1.S2.collapsed %>% select(c("Sample_ID.x", "Sample_ID.y", "BLAST_species.(as.done.by.EAV)"))
  
LLA_NS1S2_yes <- cbind(tax.NS1S2, otu.NS1.S2 %>% select(all_of(samples_yes)))
LLA_NS1S2_no <- cbind(tax.NS1S2, otu.NS1.S2 %>% select(all_of(samples_no)))

#LLA_NS1S2_temp <- LLA_NS1S2
LLA_NS1S2_yes_weeks <- tax.NS1S2
LLA_NS1S2_no_weeks <- tax.NS1S2

# Obtain the sum of all samples in the same treatment and week and assign to new asv table

# All weeks present across all plates
weeks <- c("5","6","7","8","9","10")
#weeks <- c("5","6","7","9")

for(j in weeks){
  
  # Move samples in the same community and week together
  asv.temp <- LLA_NS1S2_yes %>%
    dplyr::relocate(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)), .after = "BLAST_species.(as.done.by.EAV)")
  
  # Create a new combined column of all the samples within the community/week group
  if(length(grep(paste0(sep = ".*", sep="week", j), names(asv.temp))) > 0) {
    row_sum <- rowSums(asv.temp[,(first(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))):(last(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))), drop = F])
    
    LLA_NS1S2_yes_weeks[[paste0(sep= "week", j)]] <- row_sum
  }
}

for(j in weeks){
  
  # Move samples in the same community and week together
  asv.temp <- LLA_NS1S2_no %>%
    dplyr::relocate(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)), .after = "BLAST_species.(as.done.by.EAV)")
  
  # Create a new combined column of all the samples within the community/week group
  if(length(grep(paste0(sep = ".*", sep="week", j), names(asv.temp))) > 0) {
    row_sum <- rowSums(asv.temp[,(first(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))):(last(grep(paste0(sep = ".*", sep="week", j), names(asv.temp)))), drop = F])
    
    LLA_NS1S2_no_weeks[[paste0(sep= "week", j)]] <- row_sum
  }
}

# Relative abundance table
LLA_NS1S2_yes_weeks_norm = LLA_NS1S2_yes_weeks %>% cbind(otu.NS1.S2 %>% select(contains("inocu"))) %>%
  mutate(NS1_S2_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_NS1S2_yes_weeks_norm[, 4:8] <- sapply(LLA_NS1S2_yes_weeks_norm[, 4:8],function(x) x/sum(x) * 100)

LLA_NS1S2_no_weeks_norm = LLA_NS1S2_no_weeks %>% cbind(otu.NS1.S2 %>% select(contains("inocu"))) %>%
  mutate(NS1_S2_defined = rowSums(select(., contains("inocu")), na.rm = TRUE)) %>% select(-contains("inocu"))
LLA_NS1S2_no_weeks_norm[, 4:8] <- sapply(LLA_NS1S2_no_weeks_norm[, 4:8],function(x) x/sum(x) * 100)

# Add ASV Number column
LLA_NS1S2_yes_weeks_norm <- LLA_NS1S2_yes_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_NS1S2_yes_weeks_norm$Sample_ID.x)) %>%
  relocate(NS1_S2_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)
LLA_NS1S2_no_weeks_norm <- LLA_NS1S2_no_weeks_norm %>% mutate(ASV_Number = seq_along(LLA_NS1S2_no_weeks_norm$Sample_ID.x)) %>%
  relocate(NS1_S2_defined, .before = week5) %>% relocate(ASV_Number, .before = 1)

### Plot Alluvials. Representivity % (rel. abundance) --------------------------

inocula.list.3 <- list(NS1S2_yes = LLA_NS1S2_yes_weeks_norm, NS1S2_no = LLA_NS1S2_no_weeks_norm)

# Make melted dataframe for plotting
alluvial.plots.3 = plyr::llply(inocula.list.3, function(x){
  y = x
  
  # Add 'score' column: if asv is present in the defined inocula, return TRUE, and if the asv is present in the sample, return TRUE (e.g. TRUE.TRUE)
  y$score = interaction(rowSums(x[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0,
                        rowSums(x[,(first(grep("week", names(x)))):(last(grep("week", names(x)))), drop = F])  > 0 
  )
  
  # Add 'in_any_inoc' column: if asv is present in any of the inocula then return true
  y$in_any_inoc = rowSums(y[,(first(grep("defined", names(x)))):(last(grep("defined", names(x)))), drop = F]) > 0 
  
  y = y[order(y$in_any_inoc),]  # Move ASVs not in inocula to the top
  # y = y[order(y$score, decreasing = FALSE),]  # Uncomment this to order by presence/absence
  
  print(dim(y))  # Check dimensions
  print(head(y))
  
  # Ordering has to be done to ASV because of ggalluvial
  y$ASV_Number = seq_along(y$ASV_Number)
  y$ASV_Number = sprintf("%04d", y$ASV_Number)
  
  #colnames(y) = gsub("^X", "", colnames(y))
  #colnames(y) = gsub("week", "", colnames(y))
  #return(y)
  
  # Melt    
  y = reshape2::melt(y, measure.vars = 5:ncol(x), variable.name = "Sample", value.name = "Relative_abundance", as.is = T)
  
  y = y %>%
    mutate(Housing = Sample)
  y$Housing = gsub("^\\d.*", "GF", y$Housing)
  y$Housing = gsub("^T1D.*", "Defined inocula", y$Housing)
  
  y = tidyr::extract (y, Sample, into = c("Community", "Endpoint_temp"), regex = "(NS1|NS6|S2|S5)+(.*)", remove = FALSE)
  y = tidyr::extract(y, Endpoint_temp, into = "Endpoint", regex = "(week.*|day.*)", remove = TRUE)
  y$Endpoint = gsub("day63", "week9", y$Endpoint)
  y$Endpoint = gsub("day70", "week10", y$Endpoint)
  
  y$Endpoint = replace_na(y$Endpoint, "Defined inocula")
  y = y[y$Relative_abundance!=0,]
  
  #Comment these lines out if don't want to group by plate
  # This line sorts the plate 1 and 2 samples so that they are 'mixed' rather than separated
  # y$Sample = factor(as.character(y$Sample), levels = sort(levels(y$Sample))) 
  
  # This relevels s.t. defined community is moved to front
  # y$Sample = relevel(y$Sample, grep("T1D", levels(y$Sample), value = T, perl = T))
  
  return(y)
}
)

# Name lists in alluvial.plots list if subsetted
# names(alluvial.plots) = c("NS1", "NS6", "S2", "S5") # commented because already exists

# Convert to df
alluvial.df.3 = ldply(alluvial.plots.3)
alluvial.df.3.yes <- alluvial.plots.2[["NS1S2_yes"]]
alluvial.df.3.no <- alluvial.plots.2[["NS1S2_no"]]

# Setting custom colours for phyla plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = as.numeric(length(unique(alluvial.df.3$`BLAST_species.(as.done.by.EAV)`)))
cols = gg_color_hue(n)

names(cols) <- sort(unique(alluvial.df.3$`BLAST_species.(as.done.by.EAV)`))


p = ggplot() +
  geom_stratum(data = alluvial.df.3.yes, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df.NS1.yes, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  #facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 

p4 = ggplot() +
  geom_stratum(data = alluvial.df.3, linetype = "solid", aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance , group = interaction(`BLAST_species.(as.done.by.EAV)`, in_any_inoc), fill = `BLAST_species.(as.done.by.EAV)`, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df.3, color = "black",aes(x = Sample, stratum = `BLAST_species.(as.done.by.EAV)`, alluvium = `BLAST_species.(as.done.by.EAV)`, y = Relative_abundance, fill = `BLAST_species.(as.done.by.EAV)`)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~.id, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 
