
# Create final ASV/OTU Tables ------------------------------------------------
# This section is kept as is for both RDP or vsearch (OTU numbering arbritary)

merged_seqtab_GTDB = ceiling(merged_seqtab_GTDB)
ASV_Table_rdp_GTDB = data.frame(OTU_Num = seq_along(colnames(merged_seqtab_GTDB)), sequences = colnames(merged_seqtab_GTDB),t(merged_seqtab_GTDB), stringsAsFactors = FALSE)
row.names(ASV_Table_rdp_GTDB) = NULL
ASV_Table_rdp_GTDB$OTU_Num = sprintf("%04d", ASV_Table_rdp_GTDB$OTU_Num)
ASV_Table_rdp_GTDB$OTU_Num = paste0("OTU_", ASV_Table_rdp_GTDB$OTU_Num)

# Combine ASV table with GTDB taxonomy -----------------------------------------

taxa2.GTDB <- taxa2.GTDB %>% select(-c(Species1_temp, Species2_temp)) %>%
  separate(Species, into = c("Species", "Accession"), sep = "\\(") %>%
  mutate(Accession = substr(Accession, 1, nchar(Accession) - 1)) %>%
  mutate_at(vars(Species), ~str_replace_all(., "([A-Z][a-z]*.*) ([a-z]+.*)", "\\2")) %>% #OR just remove genus name from species
  mutate_at(vars(Species), ~str_replace_all(., "(\\d*-\\d*) ([a-z]+.*)", "\\2"))

final_asv_rdp_GTDB = cbind(ASV_Table_rdp_GTDB[,1:2], taxa2.GTDB %>% select(-contains("temp")), ASV_Table_rdp_GTDB[,3:ncol(ASV_Table_rdp_GTDB)])
row.names(final_asv_rdp_GTDB ) <- 1:nrow(final_asv_rdp_GTDB )

final_asv_rdp_GTDB = final_asv_rdp_GTDB %>%
  relocate(Kingdom:Accession, .after = sequences)

final_asv_rdp_GTDB = final_asv_rdp_GTDB  %>%
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
write.table(final_asv_rdp_GTDB , file = (paste0(errPoolName, "_RDP_GTDB_taxonomy.txt")) , sep = "\t", row.names = FALSE)

# if already generated:
# final_asv_rdp_GTDB <- read.table("/Users/alicehong/R/16S-Analysis/data/asv-tables/all5_2024_04_24_RDP_GTDB_taxonomy.txt", head = T, row.names = NULL, sep = "")
