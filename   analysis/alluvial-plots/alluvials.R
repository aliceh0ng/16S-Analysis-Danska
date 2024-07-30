# Alluvials

# Alluvials. Representivity % (rel. abundance) ---------------------------------

# Set df to either asv (raw reads), asv_norm (relative abundance), asv.weeks, or asv.weeks_norm

# df = asv_norm # run this line to have each sample show
df = asv.weeks_norm #run this line for combined by weeks

# Move defined community columns to front
df = df %>%
  dplyr::relocate(grep("defined|inocula", names(df)), .after = contains("BLASTmatch"))

# If there are NA at phyla level, replace with Unclassified
df$dummy_BLASTphylum = replace_na(df$dummy_BLASTphylum, "Unclassified")

# If there are NA at family level, replace with Other
df$dummy_BLASTfamily = replace_na(df$dummy_BLASTphylum, "Unclassified")

# Grep data from relative abundance tables for subsetting
NS1.sub = df[,grepl("ASV|BLAST|sequence|NS1", names(df))]
NS6.sub = df[,grepl("ASV|BLAST|sequence|NS6", names(df))]
S2.sub = df[,grepl("ASV|BLAST|sequences|S2", names(df))]
S5.sub = df[,grepl("ASV|BLAST|sequences|S5", names(df))]
Controls.sub = df[,grepl("ASV|BLAST|sequences|control|extraction|Pooled|defined|inocula", names(df))]

# For non-sero and seroconverter
# nonsero = df[,grepl("ASV|BLAST|sequences|NS1|NS6",names(asv_norm))]
# sero = df[,grepl("ASV|BLAST|sequences|S2|S5",names(asv_norm))]

# Moves the defined communities to the head for easy comparison
NS1.sub = NS1.sub %>%
  dplyr::relocate(grep("defined|inocula", names(NS1.sub)), .after = contains("BLASTmatch"))
NS6.sub = NS6.sub %>% 
  dplyr::relocate(grep("defined|inocula", names(NS6.sub)), .after = contains("BLASTmatch"))
S2.sub = S2.sub %>%
  dplyr::relocate(grep("defined|inocula", names(S2.sub)), .after = contains("BLASTmatch"))
S5.sub = S5.sub %>%
  dplyr::relocate(grep("defined|inocula", names(S5.sub)), .after = contains("BLASTmatch"))

# nonsero = nonsero %>%
#  dplyr::relocate(grep("inocula|defined", names(nonsero)), .after = contains("BLASTmatch"))
#  dplyr::relocate(grep("inocula&NS6", names(nonsero)), .after = contains("BLASTmatch"))
# sero = sero %>%
#  dplyr::relocate(grep("inocula|defined", names(sero)), .after = contains("BLASTmatch"))
#  dplyr::relocate(T1D_S2_defined_mouse_inocula_plate1, .after = contains("BLASTmatch"))

# Create list of all subsetted communities
# inocula.list <- list(NS1 = NS1.sub, NS6 =  NS6.sub, S2 = S2.sub, S5 = S5.sub, Controls = Controls.sub)

inocula.list <- list(NS1 = NS1.sub, NS6 =  NS6.sub, S2 = S2.sub, S5 = S5.sub)
# inocula.list <- list(NS1 = NS1.sub) #WITHOUT NS6/S5
# inocula.list <- list(Controls = Controls.sub) 

# inocula.list = list(df)

## `alluvial.plots` Make melted dataframe for plotting -------------------------

alluvial.plots = plyr::llply(inocula.list, function(x){
  y = x
  
  # Add 'score' column: if asv is present in the defined inocula, return TRUE, and if the asv is present in the sample, return TRUE (e.g. TRUE.TRUE)
  y$score = interaction(rowSums(x[,(first(grep("^T1D|defined|inocula", names(x)))):(last(grep("^T1D|defined|inocula", names(x)))), drop = F]) > 0
                        #rowSums(x[,(first(grep("^X\\d", names(x)))):(last(grep("^X\\d", names(x)))), drop = F])  > 0 # starts with X then digit
                        # rowSums(x[,(first(grep("_week", names(x)))):(last(grep("_week", names(x)))), drop = F])  > 0 # starts with X then digit
  )
  
  # Add 'in_any_inoc' column: if asv is present in ant of the inocula then return true
  y$in_any_inoc = rowSums(y[,(first(grep("^T1D|defined|inocula", names(x)))):(last(grep("^T1D|defined|inocula", names(x)))), drop = F]) > 0 
  
  y = y[order(y$in_any_inoc),]  # Move ASVs not in inocula to the top
  # y = y[order(y$score, decreasing = FALSE),]  # Uncomment this to order by presence/absence
  
  #y = y %>% mutate(dummy_BLASTphylum = fct_relevel(dummy_BLASTphylum, "Actinomycetota", "Bacteroidota", "Bdellovibrionota", "Chloroflexota", "Cyanobacteriota", "Deinococcota", "Patescibacteria", "Verrucomicrobiota", "Bacillota", "Bacillota_A", "Chlamydiota", "Unclassified"))
  
  # Ordering has to be done to ASV because of ggalluvial
  y$ASV_Number = seq_along(y$ASV_Number)
  y$ASV_Number = sprintf("%04d", y$ASV_Number)
  
  #colnames(y) = gsub("^X", "", colnames(y))
  #colnames(y) = gsub("week", "", colnames(y))
  #return(y)
  
  # Melt    
  y = reshape2::melt(y, measure.vars = 12:ncol(x), variable.name = "Sample", value.name = "Relative_abundance", as.is = T)
  
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
alluvial.df.NS1 <- alluvial.plots[["NS1"]] # To extract out just 1 community

# Setting custom colours for phyla plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = as.numeric(length(unique(alluvial.df$BLASTphylum)))
cols = gg_color_hue(n)

names(cols) <- sort(unique(alluvial.df$BLASTphylum))

# Comment:
# alluvial.plots adds the following columns: "score", "in_any_inoc", "Sample", "Community", "Endpoint", "Relative_abundance", "Housing"

# Plot Alluvials ---------------------------------------------------------------

# Using geom_stratum (this puts division between ASV, but busier)
# All communities on one plot
test <- alluvial.df %>%
  mutate(dummy_BLASTphylum = fct_relevel(dummy_BLASTphylum, "Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Unclassified"))

test2 <- alluvial.df %>%
  mutate(Phylum = fct_relevel(Phylum, "Actinomycetota", "Bacillota", "Bacillota_A", "Bacillota_C", "Bacteroidota","Pseudomonadota", "Verrucomicrobiota", "Unclassified"))

p_1 = ggplot() +
  geom_stratum(data = alluvial.df, linetype = "solid", aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(dummy_BLASTphylum, in_any_inoc), fill = dummy_BLASTphylum, alpha = in_any_inoc)) +
  geom_flow(data = alluvial.df, color = "black",aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, fill = dummy_BLASTphylum)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~Community, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 
scale_fill_manual(values = c("#D95F02", "#66A61E","#FC4E07","#00AFBB", "#A6761D", "#666666"))+
  ggtitle("Plates 1 to 5")

p_test = ggplot() +
  geom_stratum(data = test, linetype = "solid", aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(dummy_BLASTphylum, in_any_inoc), fill = dummy_BLASTphylum, alpha = in_any_inoc)) +
  geom_flow(data = test, color = "black",aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, fill = dummy_BLASTphylum)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~Community, scales = "free") +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() 
#scale_fill_manual(values = c("#D95F02", "#66A61E","#FC4E07","#00AFBB", "#A6761D", "#666666"))+
#ggtitle("Plates 1 to 5")

p_2 = ggplot() + 
  geom_stratum(data = alluvial.df, linetype = "solid", aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(Phylum, in_any_inoc), fill = Phylum, alpha = in_any_inoc)) + 
  #geom_flow(data = alluvial.df, color = "black",aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, fill = Phylum)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +   
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  #facet_wrap(~Community, scales = "free") + 
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() +
  #scale_fill_manual(values = c("#D95F02", "#1B9E77","#FFDB6D", "#7570B3", "#E7298A", "#66A61E","#C4961A", "#E6AB02", "#A6761D", "#666666")) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#7570B3", "#A6761D", "#666666", "#66A61E","#FFDB6D"))
#ggtitle("Plates 1 to 5")

p_test2 = ggplot() + 
  geom_stratum(data = test2, linetype = "solid", aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(Phylum, in_any_inoc), fill = Phylum, alpha = in_any_inoc)) + 
  geom_flow(data = test2, color = "black",aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, fill = Phylum)) +
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative abundance") +   
  scale_alpha_discrete(range = c(0.25, 1), guide = "none") + #Use (1,1) for showing only ASV contribution from inocula; if want to see all, use c(0.25, 1)
  facet_wrap(~Community, scales = "free") + 
  theme(legend.position="bottom", panel.spacing = unit(6, "lines")) +
  theme_bw() +
  scale_fill_manual(values = c("#D95F02",  "#1B9E77", "#7570B3",  "#E7298A","#66A61E", "#E6AB02", "#A6761D", "#666666")) +
  ggtitle("Plates 1 to 5")

# Using geom_bars, without stratifying ASVs (cleaner image)
defined1_nostratum = ggplot(data = alluvial.df,  aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, group = interaction(dummy_BLASTphylum, in_any_inoc),  fill = dummy_BLASTphylum)) +  
  geom_bar(position = "stack", stat = "identity") + 
  #geom_stratum(linetype = "solid", width = .7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "Relative abundance") +   
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#7570B3", "#A6761D", "#666666", "#66A61E", "#CC79A7", "#FFDB6D","#000000" )) +
  #scale_alpha_discrete(range = c(0.5, 1), guide = "none") +
  #facet_wrap(~Community, scales = "free") +  
  theme_bw() +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines"))

defined2_nostratum = ggplot(data = alluvial.df,  aes(x = Sample, stratum = ASV_Number, alluvium = ASV_Number, y = Relative_abundance, group = interaction(Phylum, in_any_inoc),  fill = Phylum)) +  
  geom_bar(position = "stack", stat = "identity") + 
  #geom_stratum(linetype = "solid", width = .7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "Relative abundance") +   
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#7570B3", "#A6761D", "#666666", "#66A61E", "#CC79A7", "#FFDB6D","#000000" )) +
  #scale_alpha_discrete(range = c(0.5, 1), guide = "none") +
  #facet_wrap(~Community, scales = "free") +  
  theme_bw() +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines"))

defined2 = ggplot() + 
  geom_stratum(data = alluvial.df, linetype = "solid", width = .7, aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(Phylum, in_any_inoc), fill = Phylum)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "Relative abundance") +   
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#7570B3", "#A6761D", "#666666", "#66A61E","#FFDB6D","#000000" )) +
  #scale_alpha_discrete(range = c(0.5, 1), guide = "none") +
  #facet_wrap(~Community, scales = "free") +  
  theme_bw() +
  theme(legend.position="bottom", panel.spacing = unit(6, "lines"))

defined1 = ggplot() +
  geom_stratum(data = alluvial.df, linetype = "solid", width = .7, aes(x = Sample, stratum = ASV_Number, y = Relative_abundance , group = interaction(dummy_BLASTphylum, in_any_inoc), fill = dummy_BLASTphylum)) +
  geom_bar(position = "stack", stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "Relative abundance") +   
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#7570B3", "#A6761D", "#666666", "#66A61E","#FFDB6D","#000000" )) +
  #scale_alpha_discrete(range = c(0.5, 1), guide = "none") +
  #facet_wrap(~Community, scales = "free") +  
  theme_bw() +
  theme(legend.position="bottom", panel.spacing = unit(1, "lines"))