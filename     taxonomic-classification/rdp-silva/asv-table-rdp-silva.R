
# Create final ASV/OTU Tables --------------------------------------------------

merged_seqtab_silva = ceiling(merged_seqtab_silva)
ASV_Table_rdp_silva = data.frame(OTU_Num = seq_along(colnames(merged_seqtab_silva)), sequences = colnames(merged_seqtab_silva),t(merged_seqtab_silva), stringsAsFactors = FALSE)
row.names(ASV_Table_rdp_silva) = NULL
ASV_Table_rdp_silva$OTU_Num = sprintf("%04d", ASV_Table_rdp_silva$OTU_Num)
ASV_Table_rdp_silva$OTU_Num = paste0("OTU_", ASV_Table_rdp_silva$OTU_Num)

## Combine ASV table with SILVA taxonomy ---------------------------------------
final_asv_rdp_silva = cbind(ASV_Table_rdp_silva[,1:2], taxa2.silva %>% select(-contains("temp")), ASV_Table_rdp_silva[,3:ncol(ASV_Table_rdp_silva)])
row.names(final_asv_rdp_silva ) <- 1:nrow(final_asv_rdp_silva )
final_asv_rdp_silva = final_asv_rdp_silva  %>%
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

## Manually specify file name below and write out ASV table --------------------
write.table(final_asv_rdp_silva , file = (paste0(errPoolName, "_RDP_SILVA_taxonomy_no_DECIPHER.txt")) , sep = "\t", row.names = FALSE)
