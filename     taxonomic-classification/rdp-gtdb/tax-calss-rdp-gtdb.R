# Taxonomic classification via GTDB -----------------------------------------
# In this version, we are assigning Taxonomy with just the genus, then adding the species. This results in 1153 unassgined species
# Another option is to assign Taxonomy with the full taxonomy, then assign species again. This results in ____ unassgined species

set.seed(100)
minBoot = 80 # Sets minimum bootstrap confidence

merged_seqtab_GTDB = merged_seqtab #Specify input file here, in correct format
colnames(merged_seqtab_GTDB) = newdf.sort$sequence #This adds back in sequences to column headers

## Initial round of taxonomic classification -----------------------------------
taxa.GTDB <- assignTaxonomy(merged_seqtab_GTDB, "/Users/alicehong/R/16S-Analysis/      data/GTDB-database/GTDB_bac120_arc53_ssu_r214_fullTaxo.fa", multithread = T, minBoot = minBoot, verbose = T) # with full taxonomy

View(taxa.GTDB)

taxa.GTDB.df = as.data.frame(taxa.GTDB)
taxa.GTDB.df = taxa.GTDB.df %>% # This part not needed if assigning with genus dataset first
  dplyr::rename(
    Species1_temp = Species
  )

na.plot.GTDB = taxa.GTDB.df %>%
  summarise_all(~(sum(is.na(.))))
na.plot.GTDB = reshape2::melt(na.plot.GTDB)

# Show that there are still unassigned ASVs, particularly at species level
na.plot.rdp.GTDB <- ggplot() + geom_bar(data = na.plot.GTDB, aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank (before additional species matching)") +
  xlab("") + ylab("") +
  annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 

## Second round of classification, improving on species assignment -------------
taxa2.GTDB <- addSpecies(taxa.GTDB.df, "/Users/alicehong/GTDB_DADA2/GTDB_bac120_arc53_ssu_r214_species.fa", tryRC = T, verbose = T)

taxa2.GTDB = taxa2.GTDB %>% # commented, not needed if species doesnt exist in first dataset
  dplyr::rename(
    Species2_temp = Species
  )

taxa2.GTDB$Species = coalesce(taxa2.GTDB$Species1_temp, taxa2.GTDB$Species2_temp)

na.plot2.GTDB = taxa2.GTDB %>%
  summarise_all(~(sum(is.na(.))))
na.plot2.GTDB = reshape2::melt(na.plot2.GTDB)

# Show that there are still unassigned ASVs, even with additional matching. But slightly improved
na.plot2.rdp.GTDB <- ggplot() + geom_bar(data = subset(na.plot2.GTDB, !grepl("temp", na.plot2.GTDB$variable)), aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank (after additional species matching)") +
  xlab("") + ylab("") +
  annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 

grid.arrange(na.plot.rdp.GTDB, na.plot2.rdp.GTDB, nrow=1)
