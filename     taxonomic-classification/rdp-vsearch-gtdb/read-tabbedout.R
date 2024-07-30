
# Once usearch_global has been performed: Read tabbedout file ------------------

tabbedout <- read.table("/Users/alicehong/R/16S-Analysis/data/vsearch/vsearch_output_april24", sep="\t", header=FALSE)
#tabbedout_uc <- read.table("", setp="\t")

taxa.vsearch <- tabbedout %>% 
  `colnames<-` (c('sequences','tax','id')) %>%
  separate(tax, into = c("Accession","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ",|=") %>%
  mutate_at(vars(Accession), funs(gsub(";tax","",.))) %>%
  mutate_at(vars(Kingdom:Species), funs(gsub("[a-z]+:", "",.))) %>% 
  select(-c(id)) %>%
  # mutate_at(vars(Species), ~str_replace_all(., "([A-Z][a-z]*)_", "\\1 ")) #removes first underscore after word starting with capital
  mutate_at(vars(Species), ~str_replace_all(., "([A-Z][a-z]*.*)_([a-z]+.*)", "\\2")) %>% #OR just remove genus name from species
  mutate_at(vars(Species), ~str_replace_all(., "(\\d*-\\d*)_([a-z]+.*)", "\\2")) %>%
  relocate(Accession, .after = Species)

# IGNORE, previous 'sintax' code -----------------------------------------------
#mutate_at(vars(Kingdom:Species), funs(gsub("[A-Za-z0-9]+:|[A-Za-z0-9]+=|[A-Za-z0-9]+;|[0-9]+.|[A-Z]+_", "", .)))

# taxa_sintax <- tabbedout %>% select(V1, V4) %>%
#   # `row.names<-`(., NULL) %>%
#   # column_to_rownames(var = "V1") %>% # These two lines convert column V1 to rownames
#   separate(V4, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ",") %>%
#   mutate_at(vars(Kingdom:Species), funs(gsub("[A-Za-z0-9]+:", "", .)))

# Note: species column now contains the genus + species. May need to remove the Genus name from this column later