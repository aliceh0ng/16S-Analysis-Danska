
# DADA2 -------------------------------------------------------------------

# Set poolList, must be character list, poolList must match samples
# Samples that are pooled should be of the same housing condition and same inocula (must make logical sense)
# Group together like samples so when they're assigned to an asv, they're biased towards one particular treatment type

poolList = list(
  "GF_NS1" = sample.names[grep("_NS1_|NS1_mouse", sample.names, fixed = F, perl = T)], # NS1
  "GF_S2" = sample.names[grep("_S2_|S2_mouse", sample.names, fixed = F, perl = T)], # S2
  
  "Chemo_NS0" = sample.names[grep("^NS0", sample.names, fixed = F, perl = T)], # NS0
  "Chemo_NS1" = sample.names[grep("^NS1FC|^NS1_Fecal", sample.names, fixed = F, perl = T)],# NS1
  "Chemo_NS4" = sample.names[grep("^NS4", sample.names, fixed = F, perl = T)], # NS4
  "Chemo_NS6" = sample.names[grep("^NS6", sample.names, fixed = F, perl = T)], # NS6
  "Chemo_S2" = sample.names[grep("^S2FC|^S2_Fecal", sample.names, fixed = F, perl = T)], # S2
  "Chemo_S3" = sample.names[grep("^S3", sample.names, fixed = F, perl = T)], # S3
  "Chemo_S5" = sample.names[grep("^S5", sample.names, fixed = F, perl = T)] # S5
  # Total 88 samples
)


# Set error pool names
errPoolName = paste0("WON","_", gsub("-", "_", Sys.Date()))


# This starts DADA2 pipeline
process(path,
        fnFs, fnRs, sample.names, errPoolName = errPoolName, orientFR.split = FALSE,
        trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = "pseudo", poolList = poolList,
        #priorsF = priorsF, priorsR = priorsR,
        priorsF = character(0), priorsR = character(0),
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        # true -> starts from where it left off, false -> starts from beginning
        preFilteredToggle = FALSE,
        preCalcErrToggle = FALSE,
        preDerepToggle = FALSE,
        preDenoiseToggle = FALSE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)
