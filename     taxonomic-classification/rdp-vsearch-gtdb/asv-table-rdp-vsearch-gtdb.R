# Create final ASV/OTU Tables --------------------------------------------------

# Combine ASV table with GTDB taxonomy, VSEARCH --------------------------------
merged_seqtab_GTDB = merged_seqtab
colnames(merged_seqtab_GTDB) = newdf.sort$sequence #This adds back in sequences to column headers

merged_seqtab_GTDB = ceiling(merged_seqtab_GTDB)
ASV_Table_rdp_GTDB = data.frame(OTU_Num = seq_along(colnames(merged_seqtab_GTDB)), sequences = colnames(merged_seqtab_GTDB),t(merged_seqtab_GTDB), stringsAsFactors = FALSE)
row.names(ASV_Table_rdp_GTDB) = NULL
ASV_Table_rdp_GTDB$OTU_Num = sprintf("%04d", ASV_Table_rdp_GTDB$OTU_Num)
ASV_Table_rdp_GTDB$OTU_Num = paste0("OTU_", ASV_Table_rdp_GTDB$OTU_Num)

final_asv_vsearch_GTDB = left_join(ASV_Table_rdp_GTDB, taxa.vsearch, by="sequences" )

na.plot.vsearch = final_asv_vsearch_GTDB %>% select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession")) %>%
  summarise_all(~(sum(is.na(.))))

na.plot.vsearch = reshape2::melt(na.plot.vsearch)

# Show that there are still unassigned ASVs, even with additional matching. But slightly improved
na.plot.vsearch.GTDB <- ggplot() + geom_bar(data = subset(na.plot.vsearch, !grepl("temp", na.plot.vsearch$variable)), aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank") +
  xlab("") + ylab("") + 
  annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 

final_asv_vsearch_GTDB = final_asv_vsearch_GTDB %>%
  relocate(Kingdom:Accession, .after = sequences) %>%
  mutate(dummy_bitscore = NA, .after = sequences) %>%
  mutate(dummy_evalue = NA, .after = sequences)%>%
  mutate(dummy_per_id = NA, .after = sequences)%>%
  mutate(dummy_BLASTmatch = NA, .after = Species) %>%
  dplyr::rename(dummy_BLASTsuperkingdom = Kingdom)%>%
  dplyr::rename(dummy_BLASTphylum = Phylum)%>%
  dplyr::rename(dummy_BLASTclass = Class)%>%
  dplyr::rename(dummy_BLASTorder = Order)%>%
  dplyr::rename(dummy_BLASTfamily = Family)%>%
  dplyr::rename(dummy_BLASTgenus = Genus)%>%
  dplyr::rename(dummy_BLASTspecies =  Species)

# Manually specify file name below and write out ASV table
# write.table(final_asv_vsearch_GTDB , file = (paste0(errPoolName, "_vsearch_GTDB_taxonomy.txt")) , sep = "\t", row.names = FALSE)

# Merging DADA2 and VSEARCH results, produce "parsed" asv table ----------------
# VSEARCH assignment to each ASV, if not it checks dada2 for assignment to any level
final_asv_combined <- final_asv_vsearch_GTDB

final_asv_rdp_silva_fill <- final_asv_rdp_silva %>% mutate(Accession = NA, .after = "dummy_BLASTmatch")

for(i in 1:nrow(final_asv_vsearch_GTDB)) {
  if(is.na(final_asv_vsearch_GTDB$dummy_BLASTKingdom[i]) == TRUE){
    final_asv_combined[i,6:14] <- final_asv_rdp_GTDB[i,6:14]
  }
}

na.plot.combined = final_asv_combined[,3:9] %>%
  summarise_all(~(sum(is.na(.))))

na.plot.combined = reshape2::melt(na.plot.combined)

# Show that there are still unassigned ASVs, even with additional matching. But slightly improved
na.plot.combined.GTDB <- ggplot() + geom_bar(data = subset(na.plot.combined, !grepl("temp", na.plot.combined$variable)), aes(x = variable, y = value), stat = "identity") +
  ggtitle("Number of unassigned ASVs at given taxonomic rank") +
  xlab("") + ylab("") + 
  #annotate(geom = "text", x = 2, y = 1500, label = paste0("minBoot = ", eval(minBoot)), color = "blue", size = 8) +
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 8) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme_bw() 


final_asv_combined = final_asv_combined %>%
  relocate(Kingdom:Accession, .after = sequences)

# Manually specify file name below and write out ASV table
write.table(final_asv_combined , file = (paste0(errPoolName, "_RDP_Silva_VSEARCH_GTDB_taxonomy_no_DECIPHER.txt")) , sep = "\t", row.names = FALSE)

# Generating Final NA Plot -----------------------------------------------------

silva <- na.plot2.silva %>% select(-contains("temp")) %>% mutate(Method = "RDP + Silva")
GTDB <- na.plot2.GTDB %>% select(-contains("temp")) %>% mutate (Method = "RDP + GTDB")
#vsearch <- na.plot.vsearch %>% mutate(Method = "VSEARCH + GTDB")
combined <- na.plot.combined %>% mutate(Method = "Combined + GTDB")
#na.plot.all <- bind_rows(silva, GTDB, vsearch, combined)
na.plot.all <- bind_rows(silva, GTDB, combined)

na.plot.all.30 <- na.plot.all
na.plot.all.30$value <- 1995 - na.plot.all.30$value

na.plot <- na.plot.all %>%
  subset(!grepl("temp|Accession", na.plot.all$variable)) %>%
  mutate(Method = fct_relevel(factor(Method), "RDP + Silva", "RDP + GTDB", "Combined + GTDB")) %>%
  
  ggplot(aes(x= variable, y = value, fill = Method)) +
  stat_summary(geom = "bar", stat= "identity" , position = position_dodge(width = 0.9), color = "black")  +
  #facet_wrap(~Strain) +
  #geom_hline(yintercept = seq(), linetype = 2, color = "red") +
  #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, 2), expand = c(0,0)) +
  scale_fill_npg() +
  theme_light() +
  ggtitle("Number of unassigned ASVs at given taxonomic rank") +
  xlab("Rank") + ylab("No. Unassigned ASVs") + 
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 5) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme(panel.border = element_rect(color = "black"),
        panel.grid = element_line(colour = "white"),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size=16,colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        axis.title = element_text(size = 12)
  )

nna.plot <- na.plot.all.30 %>%
  subset(!grepl("temp|Accession", na.plot.all.30$variable)) %>%
  mutate(Method = fct_relevel(factor(Method), "RDP + Silva", "RDP + GTDB", "Combined + GTDB")) %>%
  
  ggplot(aes(x= variable, y = value, fill = Method)) +
  stat_summary(geom = "bar", stat= "identity" , position = position_dodge(width = 0.9), color = "black")  +
  #facet_wrap(~Strain) +
  #geom_hline(yintercept = seq(), linetype = 2, color = "red") +
  #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, 2), expand = c(0,0)) +
  scale_fill_npg() +
  theme_light() +
  ggtitle("Number of assigned ASVs at given taxonomic rank") +
  xlab("Rank") + ylab("No. Assigned ASVs") + 
  annotate(geom = "text", x = 2.65, y = 1700, label = paste0("Number of collapsed ASVs: ", eval(ncol(merged_seqtab_GTDB))), colour = "red", size = 5) + 
  ylim(0,2000) + 
  geom_hline(yintercept = eval(ncol(merged_seqtab_GTDB)), color = "firebrick1", lwd = 1.5) +
  theme(panel.border = element_rect(color = "black"),
        panel.grid = element_line(colour = "white"),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size=16,colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        axis.title = element_text(size = 12)
  )
