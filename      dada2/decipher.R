
# DECIPHER::treeline - Collapse ASVs to 98.8 -----------------------------------

# This generates the cluster grouping for all ASVs
clusters = DECIPHER_clusterASVs (colnames(seqtab.merge),
                                 identityCutoff = 0.988,
                                 testBounds = FALSE, upperBound = 1, lowerBound = 1, testIncrement = 0.1,
                                 ncores = NULL)

asv_sequences <- colnames(seqtab.merge)
sample_names <- rownames(seqtab.merge)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
#aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
#d <- DECIPHER::DistanceMatrix(aln, processors = nproc)

# Generation of ASV ("OTU") table based on DECIPHER collapsing
# From the clusters vector, this collapses the seqtab file to only show collapsed ASVs at the specified percent identity

clusters <- clusters %>%
  add_column(sequence = asv_sequences)

merged_seqtab <- seqtab.merge %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t 

newdf<-as.data.frame(cbind.data.frame(rownames(clusters),clusters$cluster,clusters$sequence))
colnames(newdf)<-c("ASV","C.ASV","sequence") 
newdf.sort<-newdf %>% arrange(C.ASV) %>% distinct(C.ASV,.keep_all = TRUE)

seqdf = data.frame(OTU_Num = seq(1:nrow(newdf.sort)), sequences = newdf.sort$sequence, stringsAsFactors = F)
seqdf$OTU_Num = sprintf("%04d", seqdf$OTU_Num)
seqdf$OTU_Num = paste0("OTU_", seqdf$OTU_N)


otu = merged_seqtab %>%
  t
otu = as.data.frame(otu)
otu = tibble::rownames_to_column(otu, "OTU_Num")
otu$OTU_Num = as.numeric(otu$OTU_Num)
otu$OTU_Num = sprintf("%04d", otu$OTU_Num)
otu$OTU_Num = paste0("OTU_", otu$OTU_Num)


# Collapsed ASV table created thus far - no taxonomic names yet
otu = merge(otu, seqdf, by = "OTU_Num")
otu = otu %>% relocate(sequences, .after = OTU_Num)
