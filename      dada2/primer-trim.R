
# Begin Processing of Raw Files -------------------------------------------

# Primer Trimming ---------------------------------------------------------
# Load F/W reads
fnFs <- sort(list.files(path = demultiplexed_path, pattern = "L001_R1", full.names = TRUE))
fnFs <- grep(fnFs, pattern = "md5", invert = TRUE, value = TRUE)
fnRs <- sort(list.files(demultiplexed_path, pattern="L001_R2", full.names = TRUE))
fnRs <- grep(fnRs, pattern = "md5", invert = TRUE, value = TRUE)

# Primers for Zymo Quick 16S Plus NGS Library Prep
# Kit uses mix of 2 forward primers
# CCTACGGGDGGCWGCAG
# CCTAYGGGGYGCWGCAG
# Set most degenerate sequence as forward primer

primerF = "CCTAYGGGDBGCWGCAG"
primerR = "GACTACNVGGGTMTCTAATCC"

# This starts trimming files
# primerTrim(path, fnFs, fnRs, primerF, primerR, linkerFlen = 20, linkerRlen = 20)

# Set Primer trimmed folder path
primerTrimmed_path <- "/Users/alicehong/Desktop/SK/Experiment_Data_Results/WON26772.20240306/240229_M07519_0079_000000000-LFGPK/Primer Trimmed"

# Ensure that Primer trimmed files are read
fnFs <- sort(list.files(primerTrimmed_path, pattern="L001_R1", full.names = TRUE))
fnRs <- sort(list.files(primerTrimmed_path, pattern="L001_R2", full.names = TRUE))

# Set sample names
sample.names = sapply(strsplit(basename(fnFs), "(_S\\d+_L\\d+)", perl = T), `[` ,1)
sample.names


primerLenF = 0 # If primer is still present, set # of chars to remove
primerLenR = 0
trimLeft = c(0,0) # If linker is still present, set # of chars to remove (forward and reverse)
trimLeftSelect = c(primerLenF+trimLeft[1],primerLenR+trimLeft[2]) # How much to remove for F and R reads

priorsF = character(0)
priorsR = character(0)


# Check Quality Scores after Trimming -------------------------------------

# Legend:
# For Quality Profile 
# green = mean
# orange = median
# dashed orange = 25th and 75th quantiles

plotQualityProfile2(fnFs[1:10], aggregate = TRUE)
plotQualityProfile2(fnRs[1:10], aggregate = TRUE)

# Find shortest/longest reads for F/R
plotAggregateLengths(fnFs[1:10])
plotAggregateLengths(fnRs[1:10])

# For Plot Quality Profile 3
# first black line = start of filtered reads
# second black line = end of filtered reads
# red line = overlap with opposite reads

seqlen = 430 # Set expected 16S fragment sequence length
# ? seqlen = 280

truncLenSelect = c(260,180) # How much to keep after left aligned trimming is done.
# First, estimate via plotQualityProfile where quality drops significantly
# Then finalize these numbers after running plotQualityProfile3 for fnFs and fnRs
# lenLimF/R is set to the shortest reads as determined by plotAggregateLengths
quality_table = plotQualityProfile3(fnFs[1:10],fnRs[1:10],primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF = 270, lenLimR = 264) 

quality_table = plotQualityProfile3(fnFs[1:10],fnRs[1:10],primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF = as.numeric(names(plotAggregateLengths(fnFs[1:10])[1])), lenLimR = as.numeric(names(plotAggregateLengths(fnRs[1:10])[1]))) 

# Go to table, look for n-mer overlap (in this case, 20mer), use these two numbers and reassign into truncLenSelect
View(quality_table)

truncLenSelect = c(268, 182) 