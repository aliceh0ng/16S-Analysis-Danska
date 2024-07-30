if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version='devel')
BiocManager::install("sangeranalyseR")

library(sangeranalyseR)

# Basic function, in a 3 line wrapper
# alignment = SangerAlignment(ABIF_Directory     = "./",
#                 REGEX_SuffixForward = ".*F_user_added.ab1",
#                 REGEX_SuffixReverse = ".*R_user_added.ab1")

# Advanced config function

#parentDir <- "./"
parentDir <- "/Users/alicehong/Desktop/Missing sequences"

alignment <- SangerAlignment(inputSource          = "ABIF",
                                      processMethod        = "REGEX",
                                      ABIF_Directory       = parentDir,
                                      REGEX_SuffixForward  = "_338F_user_added.ab1",
                                      REGEX_SuffixReverse  = "_805R_user_added.ab1",
                                      TrimmingMethod       = "M1",
                                      M1TrimmingCutoff     = 0.0001,
                                      M2CutoffQualityScore = NULL,
                                      M2SlidingWindowSize  = NULL,
                                      baseNumPerRow        = 100,
                                      heightPerRow         = 200,
                                      signalRatioCutoff    = 0.33,
                                      showTrimmed          = TRUE,
                                      refAminoAcidSeq      = "",
                                      minReadsNum          = 2,
                                      minReadLength        = 20,
                                      minFractionCall      = 0.5,
                                      maxFractionLost      = 0.5,
                                      geneticCode          = GENETIC_CODE,
                                      acceptStopCodons     = TRUE,
                                      readingFrame         = 1,
                                      processorsNum        = 4)

names = names(alignment@contigList)



# Generate detailed information for each strain, F/R reads and consensus sequence
full.alignment.list = NULL

for (i in 1:length(names)){
full.alignment.list[[i]] = eval(parse (text = paste0("alignment@contigList$", "'", (names[[i]]), "'", "@alignment")))
}
names(full.alignment.list) = names

# If you want to view them in your browser
library("DECIPHER")
BrowseSeqs(full.alignment.list[[1]])


# Only get consensus sequences for each strain
library(Biostrings)

consensus.list <- DNAStringSet()
  
  #vector(mode = "list", length = length(names))

for (i in 1:length(names)){
  consensus.list[i] = DNAStringSet( eval(parse (text = paste0("alignment@contigList$", "'", (names[[i]]), "'", "@contigSeq"))) )
}
names(consensus.list) = names

# Write out consensus sequences for each strain
writeXStringSet(consensus.list, "/Users/alicehong/Desktop/Missing sequences/consensus.fa", append = F, compress = F, format = "fasta")


# Example code from package
rawDataDir <- system.file("extdata", package = "sangeranalyseR")
parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')


ACHLO_contigs <- SangerAlignment(ABIF_Directory     = parentDir,
                                 REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                 REGEX_SuffixReverse = "_[0-9]*_R.ab1$")

writeFasta(ACHLO_contigs)

