
# Taxonomic classification via SilvaDB -----------------------------------------
# In this method we're assigning taxonomy with the full (incl. species trainign set), then adding the species again
# Another option would be to assign taxonomy with the other (no speices) training set

set.seed(100) 
minBoot = 80 # Sets minimum bootstrap confidence

#merged_seqtab_silva = merged_seqtab #Specify input file here, in correct format

merged_seqtab_silva = seqtab.merge

colnames(merged_seqtab_silva) = newdf.sort$sequence #This adds back in sequences to column headers

## Initial round of taxonomic classification -----------------------------------
taxa.silva <- assignTaxonomy(merged_seqtab_silva, "/Users/alicehong/SilvaDB_138_1/silva_nr99_v138.1_wSpecies_train_set.fa", multithread = T,  minBoot = minBoot, verbose = T)

View(taxa.silva)

taxa.silva.df = as.data.frame(taxa.silva)
taxa.silva.df = taxa.silva.df %>% 
  dplyr::rename(
    Species1_temp = Species
  )
na.plot.silva = taxa.silva.df %>%
  summarise_all(~(sum(is.na(.))))
na.plot.silva = reshape2::melt(na.plot.silva)


# Show that there are still unassigned ASVs, particularly at species level 
na.plot.rdp.silva <- ggplot() + geom_bar(data = na.plot.silva, aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank (before additional species matching)") +
  xlab("") + ylab("") + 
  annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_silva))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_silva)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 

## Second round of classification, improving on species assignment -------------
taxa2.silva <- addSpecies(taxa.silva.df, "/Users/alicehong/SilvaDB_138_1/silva_species_assignment_v138.1.fa", tryRC = T, verbose = T)
taxa2.silva = taxa2.silva %>% 
  dplyr::rename(
    Species2_temp = Species
  )

taxa2.silva$Species = coalesce(taxa2.silva$Species1_temp, taxa2.silva$Species2_temp)

na.plot2.silva = taxa2.silva %>%
  summarise_all(~(sum(is.na(.))))

na.plot2.silva = reshape2::melt(na.plot2.silva)

# Show that there are still unassigned ASVs, even with additional matching. But slightly improved
na.plot2.rdp.silva <- ggplot() + geom_bar(data = subset(na.plot2.silva, !grepl("temp", na.plot2.silva$variable)), aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank (after additional species matching)") +
  xlab("") + ylab("") + 
  annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_silva))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_silva)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 

na.plots.rdp.silva <- grid.arrange(na.plot.rdp.silva, na.plot2.rdp.silva, nrow=1)
