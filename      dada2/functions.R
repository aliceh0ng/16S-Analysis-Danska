
# Functions - Select and run to load all ---------------------------------------

## Primer Trimming Function -----------------------------------------------------
library(ShortRead)

# # Forward:  CCTACGGGAGGCAGCAG
# # Reverse:  CTACHVGGGTWTCTAAT
primerTrim = function(path, fnFs, fnRs, primerF, primerR, linkerFlen = 20, linkerRlen = 20){
  primerFseq = DNAString(paste0(paste0(rep.int("N",linkerFlen), collapse = ""),primerF, collapse = ""))
  primerRseq = DNAString(paste0(paste0(rep.int("N",linkerRlen), collapse = ""),primerR, collapse = ""))
  
  dir.create(path = file.path(path,"Primer Trimmed"), showWarnings = TRUE)
  
  fnFs2 = file.path(path,"Primer Trimmed",basename(fnFs))
  fnRs2 = file.path(path,"Primer Trimmed",basename(fnRs))
  
  for(i in 1:length(fnRs2)){
    R1Stream = FastqStreamer(fnFs[i])
    R2Stream = FastqStreamer(fnRs[i])
    on.exit(close(R1Stream))
    on.exit(close(R2Stream))
    
    repeat{
      R1Chunk <- yield(R1Stream)
      R2Chunk <- yield(R2Stream)
      if (
        length(R1Chunk) == 0 |
        length(R2Chunk) == 0){
        break
      }
      
      trimF = trimLRPatterns(Lpattern = primerFseq, subject = R1Chunk, Lfixed = FALSE, ranges = FALSE)
      trimR = trimLRPatterns(Lpattern = primerRseq, subject = R2Chunk, Lfixed = FALSE, ranges = FALSE)
      
      #write to R1 file
      writeFastq(file = fnFs2[i],
                 object = trimF,
                 mode = "a"
      )
      
      #write to R2 file
      writeFastq(file = fnRs2[i],
                 object = trimR,
                 mode = "a"
      )
    }
  }
}


## Pseudolog Function -----------------------------------------------------------
pseudo_log_breaks = function(n = 5, sigma = 1, base = 10){
  force(n)
  force(sigma)
  force(base)
  n_default = n
  
  function(x, n = n_default) {
    raw_rng <- suppressWarnings(range(x, na.rm = TRUE))
    if (any(!is.finite(raw_rng))) {
      return(numeric())
    }
    pos_x = raw_rng[raw_rng>0]
    neg_x = raw_rng[raw_rng<0]
    
    if(length(pos_x)>0){
      pos_breaks = log_breaks(n, base)(c(sigma,pos_x))
    }else{
      pos_breaks = numeric()
    }
    
    if(length(neg_x)>0){
      neg_breaks = -(log_breaks(n, base)(c(sigma,abs(neg_x))))
    }else{
      neg_breaks = numeric()
    }
    
    all_breaks = c(rev(neg_breaks),pos_breaks)
    
    return(c(rev(neg_breaks),pos_breaks))
  }
}

pseudo_log10_trans = function (sigma = 1, base = 10){
  force(sigma)
  force(base)
  trans_new(name = "pseudo_log10",
            transform = function(x) asinh(x/(2 * sigma))/log(base), 
            inverse = function(x) 2 * sigma * sinh(x * log(base)),
            breaks = pseudo_log_breaks(sigma = sigma,base = base)
  )
}
# when using pseudo_log10_trans transform:
# the scale will generally have reasonable resolution from (-Inf,-sigma)|(sigma,Inf)
# with ~linear behaviour between -2*sigma and 2*sigma, (~ y = 0.1913878/sigma * x) (slope = asinh(1)/(log(10)*2))
# with ~0.4 scale units between -sigma and sigma



## filterAndTrim with Windows parallelization function --------------------------
# Parameter documentation
# see dada2::filterAndTrim()
filterAndTrimWinPara = function(fwd, filt, rev = NULL, filt.rev = NULL, compress = TRUE, 
                                truncQ = 2, truncLen = 0, trimLeft = 0, trimRight = 0, maxLen = Inf, 
                                minLen = 20, maxN = 0, minQ = 0, maxEE = Inf, rm.phix = TRUE, 
                                rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, 
                                id.sep = "\\s", id.field = NULL, multithread = FALSE, n = 1e+05, 
                                OMP = TRUE, qualityType = "Auto", verbose = FALSE) 
{
  PAIRED <- FALSE
  if (!(is.character(fwd) && is.character(filt))) 
    stop("File paths must be provided as character vectors.")
  if (length(fwd) == 1 && dir.exists(fwd)) 
    fwd <- parseFastqDirectory(fwd)
  if (!all(file.exists(fwd))) 
    stop("Some input files do not exist.")
  if (length(filt) == 1 && length(fwd) > 1) 
    filt <- file.path(filt, basename(fwd))
  if (length(fwd) != length(filt)) 
    stop("Every input file must have a corresponding output file.")
  odirs <- unique(dirname(filt))
  for (odir in odirs) {
    if (!dir.exists(odir)) {
      message("Creating output directory: ", odir)
      dir.create(odir, recursive = TRUE, mode = "0777")
    }
  }
  fwd <- normalizePath(fwd, mustWork = TRUE)
  filt <- suppressWarnings(normalizePath(filt, mustWork = FALSE))
  if (any(duplicated(filt))) 
    stop("All output files must be distinct.")
  if (any(filt %in% fwd)) 
    stop("Output files must be distinct from the input files.")
  if (!is.null(rev)) {
    PAIRED <- TRUE
    if (is.null(filt.rev)) 
      stop("Output files for the reverse reads are required.")
    if (!(is.character(rev) && is.character(filt.rev))) 
      stop("File paths (rev/filt.rev) must be provided as character vectors.")
    if (length(rev) == 1 && dir.exists(rev)) 
      rev <- parseFastqDirectory(rev)
    if (!all(file.exists(rev))) 
      stop("Some input files (rev) do not exist.")
    if (length(rev) != length(fwd)) 
      stop("Paired forward and reverse input files must correspond.")
    if (length(filt.rev) == 1 && length(rev) > 1) 
      filt.rev <- file.path(filt.rev, basename(rev))
    if (length(rev) != length(filt.rev)) 
      stop("Every input file (rev) must have a corresponding output file (filt.rev).")
    odirs <- unique(dirname(filt.rev))
    for (odir in odirs) {
      if (!dir.exists(odir)) {
        message("Creating output directory:", odir)
        dir.create(odir, recursive = TRUE, mode = "0777")
      }
    }
    rev <- suppressWarnings(normalizePath(rev, mustWork = TRUE))
    filt.rev <- suppressWarnings(normalizePath(filt.rev, 
                                               mustWork = FALSE))
    if (any(duplicated(c(filt, filt.rev)))) 
      stop("All output files must be distinct.")
    if (any(c(filt, filt.rev) %in% c(fwd, rev))) 
      stop("Output files must be distinct from the input files.")
  }
  
  # added multithreading compatibility for Windows
  if(.Platform$OS.type == "windows"){
    parallelization <- "socket"
  }else{
    parallelization <- "fork"
  }
  # if (multithread && .Platform$OS.type == "unix") {
  if(multithread){
    OMP <- FALSE
    ncores <- detectCores()
    if (is.numeric(multithread))
      ncores <- multithread
    if (is.na(ncores))
      ncores <- 1
    if (ncores > 1)
      verbose <- FALSE
  }else{
    ncores <- 1
    # if (multithread && .Platform$OS.type == "windows") {
    #   message("Multithreading has been DISABLED, as forking is not supported on .Platform$OS.type 'windows'")
  }
  if (PAIRED) {
    if(!multithread | parallelization == "fork"){
      rval <- mcmapply(fastqPairedFilter, mapply(c, fwd, rev, 
                                                 SIMPLIFY = FALSE), mapply(c, filt, filt.rev, SIMPLIFY = FALSE), 
                       MoreArgs = list(truncQ = truncQ, truncLen = truncLen, 
                                       trimLeft = trimLeft, trimRight = trimRight, 
                                       maxLen = maxLen, minLen = minLen, maxN = maxN, 
                                       minQ = minQ, maxEE = maxEE, rm.phix = rm.phix, 
                                       rm.lowcomplex = rm.lowcomplex, orient.fwd = orient.fwd, 
                                       matchIDs = matchIDs, id.sep = id.sep, id.field = id.field, 
                                       n = n, OMP = OMP, qualityType = qualityType, 
                                       compress = compress, verbose = verbose), mc.cores = ncores, 
                       mc.silent = TRUE)
    }
    if(parallelization == "socket"){
      fn <- mapply(c, fwd, rev, SIMPLIFY = FALSE)
      fout <- mapply(c, filt, filt.rev, SIMPLIFY = FALSE)
      mat <- cbind(fn = fn, fout = fout)
      
      # starting parallel backend
      cl <- makeCluster(ncores, type = "PSOCK")
      registerDoParallel(cl)
      
      rval <- maply(.data = mat, .expand = FALSE, .fun = fastqPairedFilter,
                    truncQ = truncQ, truncLen = truncLen, 
                    trimLeft = trimLeft, trimRight = trimRight, 
                    maxLen = maxLen, minLen = minLen, maxN = maxN, 
                    minQ = minQ, maxEE = maxEE, rm.phix = rm.phix, 
                    rm.lowcomplex = rm.lowcomplex, orient.fwd = orient.fwd, 
                    matchIDs = matchIDs, id.sep = id.sep, id.field = id.field, 
                    n = n, OMP = OMP, qualityType = qualityType, 
                    compress = compress, verbose = verbose,
                    .parallel = TRUE, .paropts = list(.packages = "dada2"))
      
      #stopping parallel backend
      stopCluster(cl)
      
      rval <- t(rval)
    }
  }
  else {
    if(!multithread | parallelization == "fork"){
      rval <- mcmapply(fastqFilter, fwd, filt, MoreArgs = list(truncQ = truncQ, 
                                                               truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
                                                               maxLen = maxLen, minLen = minLen, maxN = maxN, minQ = minQ, 
                                                               maxEE = maxEE, rm.phix = rm.phix, rm.lowcomplex = rm.lowcomplex, 
                                                               orient.fwd = orient.fwd, n = n, OMP = OMP, qualityType = qualityType, 
                                                               compress = compress, verbose = verbose), mc.cores = ncores, 
                       mc.silent = TRUE)
    }
    if(parallelization == "socket"){
      mat <- cbind(fn = fwd, fout = filt)
      
      # starting parallel backend
      cl <- makeCluster(ncores, type = "PSOCK")
      registerDoParallel(cl)
      
      rval <- maply(.data = mat, .expand = FALSE, .fun = fastqFilter,
                    truncQ = truncQ, 
                    truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
                    maxLen = maxLen, minLen = minLen, maxN = maxN, minQ = minQ, 
                    maxEE = maxEE, rm.phix = rm.phix, rm.lowcomplex = rm.lowcomplex, 
                    orient.fwd = orient.fwd, n = n, OMP = OMP, qualityType = qualityType, 
                    compress = compress, verbose = verbose,
                    .parallel = TRUE,.inform = TRUE, .paropts = list(.packages = "dada2"))
      
      #stopping parallel backend
      stopCluster(cl)
      
      rval <- t(rval)
      colnames(rval) <- fwd
    }
  }
  if (!is(rval, "matrix")) {
    if (is(rval, "list")) {
      rval <- unlist(rval[sapply(rval, is.character)])
    }
    if (length(rval) > 5) 
      rval <- rval[1:5]
    stop("These are the errors (up to 5) encountered in individual cores...\n", 
         rval)
  }
  if (ncol(rval) != length(fwd)) {
    stop("Some input files were not processed, perhaps due to memory issues. Consider lowering ncores.")
  }
  colnames(rval) <- basename(fwd)
  if (all(rval["reads.out", ] == 0)) {
    warning("No reads passed the filter. Please revisit your filtering parameters.")
  }
  else if (any(rval["reads.out", ] == 0)) {
    message("Some input samples had no reads pass the filter.")
  }
  return(invisible(t(rval)))
}

## DADA2 Processing Pipeline function --------------------------------------
# Parameter documentation 
# Path parameters
# path : filepath to folder into which new folders/outputs will be written

# Sample and pipeline parameters 
# fnFs, fnRs : vector of paths to FASTQ files (gzipped is OK)

# sample.names : vector of sample names (used to name and save sequence table rows)

# errPoolName : name of sequencing run (used to name and save error model files)

# orientFR.split = FALSE : Are samples reads split into F and R orientations R1 and R2 files be used separately for error models?
# Split samples will be recombined at the end in the seqtab.nochim.
# If TRUE, expect files to be named "sampleName.orientF_R1[...]", "sampleName.orientF_R2[...]", "sampleName.orientR_R1[...]", "sampleName.orientR_R2[...]";
# where orientF/R designates the orientation of the amplicon, and R1/R2 designates the read index.

# Filter and trimming parameters 
# trimLeftSelect, truncLenSelect : vectors of length 2, # of nucleotides to trim from 5' end and total length to truncate reads at for R1 and R2

# filter.matchIDs = FALSE : should matching IDs for fwd and rev reads be checked and used to filter reads?

# Error model parameters
# ErrModelMonotonicity : should error model monotonicity be enforced? Do this if sequencing run has collapsed Q score binning.

# DADA2 sensitivity, pooling and priors parameters
# OMEGA_A = getDadaOpt(option = "OMEGA_A") : DADA2 denoising sensitivity parameter.
# Sets threshold for the creation of a new partition with a significantly overabundant sequence as the center. See DADA2 documentation for details.

# pool = FALSE : (FALSE, TRUE or "pseudo") DADA2 denoising pooling option. See DADA2 documentation for details.

# poolList = NULL : list of character vectors containing sample names of samples in each pool.
# Each list element should be named with the corresponding dadaPoolName.
# If orientFR.split = TRUE, should contain list elements for each of "dadaPoolName.orientF" and "dadaPoolName.orientR".
# Samples which are not named in poolList, but are in FnFs, FnRs, and sample.names will denoised unpooled.

# priorsF = character(0) : list of character vectors containing priors for forward reads, must be exactly same length as expected reads. 
# Each list element should be named with the corresponding sample name.
# If orientFR.split = TRUE, should contain list elements for each of "sampleName.orientF" and "sampleName.orientR" in the R1 read index.
# If pooling is used (pool = TRUE or "pseudo"), should contain list elements for each pool.
# If both orientFR.split = TRUE and pooling is used, should contain list elements for each pool "dadaPoolName.orientF" and "dadaPoolName.orientR" in the R1 read index.

# priorsR = character(0) : list of character vectors containing priors for reverse reads, must be exactly same length as expected reads.
# Each list element should be named with the corresponding sample name.
# If orientFR.split = TRUE, should contain list elements for each of "sampleName.orientF" and "sampleName.orientR" in the R2 read index.
# If pooling is used (pool = TRUE or "pseudo"), should contain list elements for each pool.
# If both orientFR.split = TRUE and pooling is used, should contain list elements for each pool "dadaPoolName.orientF" and "dadaPoolName.orientR" in the R2 read index.

# Plot and verbosity parameters
# plotToggle = FALSE : should graphs summarizing the read quality and error models be printed?
# printPriorHead = FALSE : should some fwd and rev reads and priors for each sample be printed to help troubleshoot trimming priors

# Preprocessing parameters
# preFilteredToggle = TRUE : has filtering reads already been done? Requires filtered FASTQ files in a "filtered" folder.
# preCalcErrToggle = TRUE :  has error model generation been done? Requires error model files, named by errPoolName, in a "Error Models" folder.
# preDerepToggle = TRUE :  has dereplicating reads already been done? Requires dereplicated reads files in a "dereplicated" folder.
# preDenoiseToggle = TRUE : has denoising reads already been done? Requires denoised reads files in a "dada denoised" folder.

# Data saving parameters
# saveFiltered = TRUE : should filtered FASTQ files be saved? If not, files will be removed after dereplication.
# saveDereplicated = TRUE : should dereplicated reads files be saved? If not, files will be removed after merging.
# saveDenoised = TRUE : should DADA denoised files be saved? If not, files will be removed after merging.
process = function(path,
                   fnFs, fnRs, sample.names, errPoolName, orientFR.split = FALSE,
                   trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
                   ErrModelMonotonicity = FALSE,
                   OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = character(0),
                   priorsF = character(0), priorsR = character(0),
                   
                   plotToggle = FALSE,
                   printPriorHead = FALSE,
                   
                   preFilteredToggle = TRUE,
                   preCalcErrToggle = TRUE,
                   preDerepToggle = TRUE,
                   preDenoiseToggle = TRUE,
                   
                   saveFiltered = TRUE,
                   saveDereplicated = TRUE,
                   saveDenoised = TRUE){
  # ---- filter and trimming ----
  cat("Filter and trimming...\n")
  dir.create(path = paste0(path,"/filtered/"), showWarnings = TRUE)
  
  # Set destination to filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  if(!preFilteredToggle){
    #filter and trim function
    #NOTE: consider relaxing maxEE if needed
    out <- filterAndTrimWinPara(fnFs, filtFs, fnRs, filtRs, truncLen = truncLenSelect, trimLeft = trimLeftSelect,
                                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                compress=TRUE, multithread=TRUE, matchIDs=filter.matchIDs)
    
    # creating sequence tracking table
    track <- cbind(out, NA, NA, NA, NA)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    
    print(track)
    
  }else{
    # creating sequence tracking table
    track <- cbind(rep(NA,length(sample.names)), rep(NA,length(sample.names)), rep(NA,length(sample.names)), rep(NA,length(sample.names)))
    colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
  }
  
  #plotting quality score profiles
  if(plotToggle){
    plotQualityProfile(filtFs,aggregate = TRUE)
    plotQualityProfile(filtRs,aggregate = TRUE)
  }
  
  # ---- generating error models ----
  cat("\n")
  cat("Generating error models...\n")
  dir.create(path = paste0(path,"/Error Models"), showWarnings = TRUE)
  
  if(!preCalcErrToggle){
    if(orientFR.split){
      #learning error models for split samples
      orientFerrF <- learnErrors(filtFs[grep(".orientF", filtFs, fixed = TRUE)], multithread=TRUE)
      orientFerrR <- learnErrors(filtRs[grep(".orientF", filtRs, fixed = TRUE)], multithread=TRUE)
      orientRerrF <- learnErrors(filtFs[grep(".orientR", filtFs, fixed = TRUE)], multithread=TRUE)
      orientRerrR <- learnErrors(filtRs[grep(".orientR", filtRs, fixed = TRUE)], multithread=TRUE)
      
      if(ErrModelMonotonicity){
        #set error values for Q-scores<40 to equal error values for Q-scores=40
        #orientFerrF
        orientFerrFOutMono = (getErrors(orientFerrF))
        orientFerrFOutMono = apply(orientFerrFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientFerrFOutMono[,40]] = orientFerrFOutMono[,40][y<orientFerrFOutMono[,40]]
          print(y)
          return(y)
        })
        orientFerrF$err_out = orientFerrFOutMono
        
        #orientFerrR
        orientFerrROutMono = (getErrors(orientFerrR))
        orientFerrROutMono = apply(orientFerrROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientFerrROutMono[,40]] = orientFerrROutMono[,40][y<orientFerrROutMono[,40]]
          print(y)
          return(y)
        })
        orientFerrR$err_out = orientFerrROutMono
        
        #orientRerrF
        orientRerrFOutMono = (getErrors(orientRerrF))
        orientRerrFOutMono = apply(orientRerrFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientRerrFOutMono[,40]] = orientRerrFOutMono[,40][y<orientRerrFOutMono[,40]]
          print(y)
          return(y)
        })
        orientRerrF$err_out = orientRerrFOutMono
        
        #orientRerrR
        orientRerrROutMono = (getErrors(orientRerrR))
        orientRerrROutMono = apply(orientRerrROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientRerrROutMono[,40]] = orientRerrROutMono[,40][y<orientRerrROutMono[,40]]
          print(y)
          return(y)
        })
        orientRerrR$err_out = orientRerrROutMono
      }
      
      #saving error model files as .rds
      saveRDS(orientFerrF,file = paste0(path,"/Error Models/",errPoolName,".orientF","_errF.rds"))
      saveRDS(orientFerrR,file = paste0(path,"/Error Models/",errPoolName,".orientF","_errR.rds"))
      saveRDS(orientRerrF,file = paste0(path,"/Error Models/",errPoolName,".orientR","_errF.rds"))
      saveRDS(orientRerrR,file = paste0(path,"/Error Models/",errPoolName,".orientR","_errR.rds"))
      
    }else{
      #learning error models
      errF <- learnErrors(filtFs, multithread=TRUE)
      errR <- learnErrors(filtRs, multithread=TRUE)
      
      if(ErrModelMonotonicity){
        #set error values for Q-scores<40 to equal error values for Q-scores=40
        #ErrF
        errFOutMono = (getErrors(errF))
        errFOutMono = apply(errFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<errFOutMono[,40]] = errFOutMono[,40][y<errFOutMono[,40]]
          print(y)
          return(y)
        })
        
        errF$err_out = errFOutMono
        
        #ErrR
        errROutMono = (getErrors(errR))
        errROutMono = apply(errROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<errROutMono[,40]] = errROutMono[,40][y<errROutMono[,40]]
          print(y)
          return(y)
        })
        
        errR$err_out = errROutMono
      }
      
      #saving error model files as .rds
      saveRDS(errF,file = paste0(path,"/Error Models/",errPoolName,"_errF.rds"))
      saveRDS(errR,file = paste0(path,"/Error Models/",errPoolName,"_errR.rds"))
    }    
  }else{
    if(orientFR.split){
      orientFerrF = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientF_errF.rds"))
      orientFerrR = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientF_errR.rds"))
      orientRerrF = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientR_errF.rds"))
      orientRerrR = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientR_errR.rds"))
    }else{
      errF = readRDS(file = paste0(path,"/Error Models/",errPoolName,"_errF.rds"))
      errR = readRDS(file = paste0(path,"/Error Models/",errPoolName,"_errR.rds"))
    }
  }
  
  #plotting error model summary graphs
  if(plotToggle){
    if(orientFR.split){
      print(plotErrors(orientFerrF, nominalQ=TRUE))
      print(plotErrors(orientFerrR, nominalQ=TRUE))
      print(plotErrors(orientRerrF, nominalQ=TRUE))
      print(plotErrors(orientRerrR, nominalQ=TRUE))
    }else{
      print(plotErrors(errF, nominalQ=TRUE))
      print(plotErrors(errR, nominalQ=TRUE))
    }
  }
  
  # ---- dereplicating filtered FASTQ files ----
  cat("\n")
  cat("Dereplicating sequences...\n")
  dir.create(path = paste0(path,"/dereplicated"), showWarnings = TRUE)
  
  if(!preDerepToggle){
    # dereplication is done multicore parallel, if possible
    # generating parallel backend
    nodes <- detectCores()
    cl <- makeCluster(nodes, type = "PSOCK")
    registerDoParallel(cl)
    
    # dereplicating fwd reads
    derepFs = aaply(filtFs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_F_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepF <- derepFastq(x, verbose=TRUE)
      saveRDS(derepF,file = paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
    }, filePath = path, .parallel = TRUE, .paropts = list(.packages = "dada2"))
    
    # dereplicating rev reads
    derepRs = aaply(filtRs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_R_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepR <- derepFastq(x, verbose=TRUE)
      saveRDS(derepR,file = paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
    }, filePath = path, .parallel = TRUE, .paropts = list(.packages = "dada2"))
    
    #stopping parallel backend
    stopCluster(cl)
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    # Remove filtered files if saveFiltered == FALSE and the process() generated filtered files (preFilteredToggle == FALSE)
    if(!saveFiltered & !preFilteredToggle){
      file.remove(filtFs,filtRs)
    }
  }else{
    derepFs = paste0(path,"/dereplicated/",sample.names,"_derepF.rds")
    derepRs = paste0(path,"/dereplicated/",sample.names,"_derepR.rds")
  }
  
  # ---- denoising dereplicated files ----
  cat("\n")
  cat("Denoising sequences...\n")
  dir.create(path = paste0(path,"/dada denoised"), showWarnings = TRUE)
  
  if(!preDenoiseToggle){
    dadaFs = data.frame(name = rep(NA,length(fnFs)),
                        count = rep(NA,length(fnFs)))
    
    dadaRs = data.frame(name = rep(NA,length(fnRs)),
                        count = rep(NA,length(fnRs)))
    
    # if pool != FALSE, denoise pooled samples according to poolList,
    # else set outPoolListNamesMatch (indices of non-pooled samples) to all samples
    if(pool != FALSE){
      
      # check if poolList is NULL
      if(is.null(poolList)){
        stop(("If pooling, poolList must not be NULL."))
      }
      
      # all samples in sample.names will be processed; 
      # samples in poolList will be pooled according to poolList,
      # samples not in poolList will default to non-pooled denoising
      
      # determining samples in/out poolList
      inPoolListNamesMatch = unlist(poolList, use.names = FALSE)
      
      outPoolListNamesMatch = which(!(sample.names %in% inPoolListNamesMatch))
      inPoolListNamesMatch = which(sample.names %in% inPoolListNamesMatch)
      
      a_ply(names(poolList),1, function(dadaPoolName){
        cat(paste0("denoising pool '",dadaPoolName,"'...\n"))
        
        poolSampleNames = poolList[[dadaPoolName]]
        poolSampleNamesMatch = match(poolSampleNames,sample.names)
        
        poolDerepFs = derepFs[poolSampleNamesMatch]
        poolDerepRs = derepRs[poolSampleNamesMatch]
        
        print(poolSampleNames)
        
        if(!isEmpty(priorsF)){
          priorMatchDf = adply(poolDerepFs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
            
            derepFFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
              cat("head of dereplicated Fwd sequences:\n")
              print(head(names(derepFFile$uniques)))
              cat("head of prior Fwd sequences:\n")
              print(head(priorsF[[dadaPoolName]]))
              
              # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
              cat("\n")
              cat(paste0(sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[dadaPoolName]]), " Fwd priors\n"))
              cat(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences\n"))
            }
            
            derepInPriorsCount = sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]])
            priorsCount = length(priorsF[[dadaPoolName]])
            derepInPriorsReadCount = sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]])
            readCount = sum(derepFFile$uniques)
            
            # same data as above in tabular format
            df = data.frame(sampleName = sampleName,
                            derepInPriorsCount = derepInPriorsCount,
                            priorsCount = priorsCount,
                            derepInPriorsReadCount = derepInPriorsReadCount,
                            readCount = readCount)
            
            return(df)
          },.progress = "text")
          print(priorMatchDf)
        }else{
          cat("no fwd priors\n")
        }
        
        # Loading Derep files into a list
        poolDerepFFileList = alply(poolDerepFs,1,function(x){
          derepFFile = readRDS(x)
        })
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", poolSampleNames, fixed = TRUE)
          if(!(all(is.orientR) | all(!is.orientR))){
            stop("If orientFR.split and pooling, pools should contain only samples with same orientation")
          }
          is.orientR = all(is.orientR)
          
          if(!is.orientR){
            dadaErrF = orientFerrF
          }else{
            dadaErrF = orientRerrF
          }
        }else{
          dadaErrF = errF
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsF)){
          dadaPriorsF = priorsF[[dadaPoolName]]
          if(is.null(dadaPriorsF)){
            dadaPriorsF = character(0)
          }
        }else{
          dadaPriorsF = character(0)
        }
        
        #DADA denoising
        dadaFList <- dada(poolDerepFFileList, err=dadaErrF, priors = dadaPriorsF, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaFList)){
          saveRDS(dadaFList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds"))
          dadaFs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds")
          dadaFs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaFList[[i]]))
        }
        
        track[,"denoisedF"] = as.numeric(dadaFs[,2])
        cat("\n")
        print(track[poolSampleNamesMatch,])
        
        if(!isEmpty(priorsR)){
          priorMatchDf = adply(poolDerepRs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
            
            derepRFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
              cat("\n")
              cat("head of dereplicated Rev sequences:\n")
              print(head(names(derepRFile$uniques)))
              cat("head of prior Rev sequences:\n")
              print(head(priorsR[[dadaPoolName]]))
              
              # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
              cat("\n")
              cat(paste0(sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]), " Rev dereplicated priors matches out of ", length(priorsR[[dadaPoolName]]), " Rev priors\n"))
              cat(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences\n"))
            }
            
            derepInPriorsCount = sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]])
            priorsCount = length(priorsR[[dadaPoolName]])
            derepInPriorsReadCount = sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]])
            readCount = sum(derepRFile$uniques)
            
            # same data as above in tabular format
            df = data.frame(sampleName = sampleName,
                            derepInPriorsCount = derepInPriorsCount,
                            priorsCount = priorsCount,
                            derepInPriorsReadCount = derepInPriorsReadCount,
                            readCount = readCount)
            
            return(df)
          },.progress = "text")
          print(priorMatchDf)
        }else{
          cat("no rev priors\n")
        }
        
        # Loading Derep files into a list
        poolDerepRFileList = alply(poolDerepRs,1,function(x){
          derepRFile = readRDS(x)
        })
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", poolSampleNames, fixed = TRUE)
          if(!(all(is.orientR) | all(!is.orientR))){
            stop("If orientFR.split and pooling, pools should contain only samples with same orientation")
          }
          is.orientR = all(is.orientR)
          
          if(!is.orientR){
            dadaErrR = orientFerrR
          }else{
            dadaErrR = orientRerrR
          }
        }else{
          dadaErrR = errR
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsR)){
          dadapriorsR = priorsR[[dadaPoolName]]
          if(is.null(dadapriorsR)){
            dadapriorsR = character(0)
          }
        }else{
          dadapriorsR = character(0)
        }
        
        #DADA denoising
        dadaRList <- dada(poolDerepRFileList, err=dadaErrR, priors = dadapriorsR, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaRList)){
          saveRDS(dadaRList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds"))
          dadaRs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds")
          dadaRs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaRList[[i]]))
        }
        
        track[,"denoisedR"] = as.numeric(dadaRs[,2])
        cat("\n")
        print(track[poolSampleNamesMatch,])
      })
      
      cat("\n")
      print(track[inPoolListNamesMatch,])
    }else{
      outPoolListNamesMatch = 1:length(sample.names)
    }
    
    # If any samples are non-pooled, denoise them
    if(length(outPoolListNamesMatch) > 0){
      if(!isEmpty(priorsF)){
        priorMatchDf = adply(derepFs[outPoolListNamesMatch],1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
          
          derepFFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
            cat("\n")
            cat("head of dereplicated Fwd sequences:")
            print(head(names(derepFFile$uniques)))
            cat("head of prior Fwd sequences:")
            print(head(priorsF[[sampleName]]))
            
            # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
            cat("\n")
            cat(paste0(sum(names(derepFFile$uniques) %in% priorsF[[sampleName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[sampleName]]), " Fwd priors\n"))
            cat(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[sampleName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepFFile$uniques) %in% priorsF[[sampleName]])
          priorsCount = length(priorsF[[sampleName]])
          derepInPriorsReadCount = sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[sampleName]]])
          readCount = sum(derepFFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
      }else{
        cat("no fwd priors\n")
      }
      
      #DADA denoising
      dadaFs[outPoolListNamesMatch,] = aaply(derepFs[outPoolListNamesMatch],1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
        
        derepFFile = readRDS(x)
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", sampleName, fixed = TRUE)
          if(!is.orientR){
            dadaErrF = orientFerrF
          }else{
            dadaErrF = orientRerrF
          }
        }else{
          dadaErrF = errF
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsF)){
          dadaPriorsF = priorsF[[sampleName]]
          if(is.null(dadaPriorsF)){
            dadaPriorsF = character(0)
          }
        }else{
          dadaPriorsF = character(0)
        }
        
        cat("\n")
        dadaF <- dada(derepFFile, err=dadaErrF, priors = dadaPriorsF, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaF,file = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds")
        count = sum(getUniques(dadaF))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedF"] = as.numeric(dadaFs[,2])
      print(track[outPoolListNamesMatch,])
      
      if(!isEmpty(priorsR)){
        priorMatchDf = adply(derepRs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
          
          derepRFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
            cat("\n")
            cat("head of dereplicated Rev sequences:\n")
            print(head(names(derepRFile$uniques)))
            cat("head of prior Rev sequences:\n")
            print(head(priorsR[[sampleName]]))
            
            # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
            cat(" ")
            cat(paste0(sum(names(derepRFile$uniques) %in% priorsR[[sampleName]]), " Rev dereplicated priors matches out of ", length(priorsR[[sampleName]]), " Rev priors\n"))
            cat(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[sampleName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepRFile$uniques) %in% priorsR[[sampleName]])
          priorsCount = length(priorsR[[sampleName]])
          derepInPriorsReadCount = sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[sampleName]]])
          readCount = sum(derepRFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
      }else{
        cat("no rev priors\n")
      }
      
      #DADA denoising
      dadaRs[outPoolListNamesMatch,] = aaply(derepRs[outPoolListNamesMatch],1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
        
        derepRFile = readRDS(x)
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", sampleName, fixed = TRUE)
          if(!is.orientR){
            dadaErrR = orientFerrR
          }else{
            dadaErrR = orientRerrR
          }
        }else{
          dadaErrR = errR
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsR)){
          dadaPriorsR = priorsR[[sampleName]]
          if(is.null(dadaPriorsR)){
            dadaPriorsR = character(0)
          }
        }else{
          dadaPriorsR = character(0)
        }
        
        cat("\n")
        dadaR <- dada(derepRFile, err=dadaErrR, priors = dadaPriorsR, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaR,file = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds")
        count = sum(getUniques(dadaR))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedR"] = as.numeric(dadaRs[,2])
      print(track[outPoolListNamesMatch,])
    }
    print(track)
    
  }else{
    dadaFs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaF.rds"))
    dadaRs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaR.rds"))
  }
  
  # ---- merging denoised F and R files ----
  cat("\n")
  cat("Merging F and R sequences...\n")
  dir.create(path = paste0(path,"/outputs"), showWarnings = TRUE)
  
  mergerArgs = data.frame(dadaFs = dadaFs[,1], derepFs, dadaRs = dadaRs[,1], derepRs, stringsAsFactors = FALSE)
  
  # merging fwd and rev sequences is done multicore parallel, if possible
  # generating parallel backend
  nodes <- detectCores()
  cl <- makeCluster(nodes, type = "PSOCK")
  registerDoParallel(cl)
  
  #merging fwd and rev sequences
  mergers = mlply(mergerArgs, function(dadaFs, derepFs, dadaRs, derepRs){
    dadaF = readRDS(dadaFs)
    derepF = readRDS(derepFs)
    dadaR = readRDS(dadaRs)
    derepR = readRDS(derepRs)
    
    merger <- mergePairs(dadaF, derepF, dadaR, derepR, maxMismatch = 0, verbose=TRUE)
    return(merger)
  }, .parallel = TRUE, .paropts = list(.packages = "dada2"))
  
  #stopping parallel backend
  stopCluster(cl)
  
  names(mergers) = sample.names
  
  # Remove dereplicated files if saveDereplicated == FALSE and the process() generated dereplicated files (preDerepToggle == FALSE)
  if(!saveDereplicated & !preDerepToggle){
    file.remove(derepFs,derepRs)
  }
  
  # Remove denoised files if saveDenoised == FALSE and the process() generated denoised files (preDenoiseToggle == FALSE)
  if(!saveDenoised & !preDenoiseToggle){
    file.remove(dadaFs[,1],dadaRs[,1])
  }
  
  #saving mergers
  for(i in 1:length(sample.names)){
    saveRDS(mergers[[i]],file = paste0(path,"/outputs/",sample.names[i],"_mergers.rds"))
  }
  
  track[,"merged"] = sapply(mergers, function(x) sum(getUniques(x)))
  print(track)
  
  # saving sequence tables
  for(i in 1:length(sample.names)){
    seqtabIndiv <- makeSequenceTable(mergers[[i]])
    saveRDS(seqtabIndiv,file = paste0(path,"/outputs/",sample.names[i],"_seqtab.rds"))
  }
  
  cat("Removing bimeras...\n")
  if(orientFR.split){
    seqtab.orientF <- makeSequenceTable(mergers[grep(".orientF", sample.names, fixed = TRUE)])
    saveRDS(seqtab.orientF,file = paste0(path,"/outputs/",errPoolName,".orientF_seqtab.rds"))
    
    seqtab.orientR <- makeSequenceTable(mergers[grep(".orientR", sample.names, fixed = TRUE)])
    saveRDS(seqtab.orientR,file = paste0(path,"/outputs/",errPoolName,".orientR_seqtab.rds"))
    
    #consensus-based chimera removal (sequences which appear to be composed of two parent sequences)
    seqtab.orientF.nochim <- removeBimeraDenovo(seqtab.orientF, method="consensus", multithread=TRUE, verbose=TRUE)
    seqtab.orientR.nochim <- removeBimeraDenovo(seqtab.orientR, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.orientF.nochim,file = paste0(path,"/outputs/",errPoolName,".orientF_seqtab.nochim.rds"))
    saveRDS(seqtab.orientR.nochim,file = paste0(path,"/outputs/",errPoolName,".orientR_seqtab.nochim.rds"))
    
    rownames(seqtab.orientF.nochim) = gsub(".orientF","",rownames(seqtab.orientF.nochim))
    rownames(seqtab.orientR.nochim) = gsub(".orientR","",rownames(seqtab.orientR.nochim))
    colnames(seqtab.orientR.nochim) = as.character(reverseComplement(DNAStringSet(colnames(seqtab.orientR.nochim))))
    
    seqtab.nochim <- mergeSequenceTables(tables = list(seqtab.orientF.nochim,seqtab.orientR.nochim), repeats = "sum")
    seqtab.nochim <- removeBimeraDenovo(seqtab.nochim, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.nochim,file = paste0(path,"/outputs/",errPoolName,"_seqtab.nochim.rds"))
    
    track[grep(".orientF", sample.names, fixed = TRUE),"nonchim"] = rowSums(seqtab.orientF.nochim)
    track[grep(".orientR", sample.names, fixed = TRUE),"nonchim"] = rowSums(seqtab.orientR.nochim)
  }else{
    seqtab <- makeSequenceTable(mergers)
    saveRDS(seqtab,file = paste0(path,"/outputs/",errPoolName,"_seqtab.rds"))
    
    #consensus-based chimera removal (sequences which appear to be composed of two parent sequences)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.nochim,file = paste0(path,"/outputs/",errPoolName,"_seqtab.nochim.rds"))
    
    track[,"nonchim"] = rowSums(seqtab.nochim)
  }
  
  cat("\n")
  cat("Done.\n")
  return(track)
}

## plotQualityProfile2 ----------------------------------------------------------
# edited by Chris
library(BiocParallel)

plotQualityProfile2 = function (fl, n = 5e+05, aggregate = FALSE) 
{
  ncores = detectCores()
  
  if(.Platform$OS.type == "windows"){
    cl <- makeCluster(ncores, type = "PSOCK")
  }else{
    cl <- makeCluster(ncores, type = "FORK")
  }
  
  clusterExport(cl, varlist = c("n"))
  
  test = parLapply(cl = cl, X = fl[!is.na(fl)], function(f){
    serialParam = BiocParallel::SerialParam()
    srqa <- ShortRead::qa(f, n = n, BPPARAM = serialParam)
    df <- cbind(srqa[["perCycle"]]$quality,file = basename(f))
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }
    else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), 
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    
    df
    
    
    statdf <- data.frame(Cycle = as.integer(rownames(means)),
                         Mean = means,
                         Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                         Q75 = as.vector(q75s), 
                         Cum = 10 * as.vector(cums)/min(rc,n), 
                         file = basename(f))
    
    anndf <- data.frame(minScore = min(df$Score),
                        label = basename(f), rclabel = rclabel, rc = rc,
                        file = basename(f))
    
    out = list(plotdf = df, statdf = statdf, anndf = anndf)
    out
  })
  
  stopCluster(cl)
  
  plotdf = lapply(test, function(x){
    x[[1]]
  })
  plotdf = do.call(rbind, plotdf)
  
  statdf = lapply(test, function(x){
    x[[2]]
  })
  statdf = do.call(rbind, statdf)
  
  anndf = lapply(test, function(x){
    x[[3]]
  })
  anndf = do.call(rbind, anndf)
  
  anndf$minScore <- min(anndf$minScore)
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, 
                                sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, 
                    plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                                 plotdf.summary$Cycle)
    
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                                 Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                 Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/sum(pmin(anndf$rc, 
                                                                                            n)))
    p <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) + 
      geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", 
                                                         high = "black") + geom_line(data = statdf.summary, 
                                                                                     aes(y = Mean), color = "#66C2A5") + geom_line(data = statdf.summary, 
                                                                                                                                   aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", 
                size = 0.25) + geom_line(data = statdf.summary, 
                                         aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      ylab("Quality Score") + xlab("Cycle") + annotate("text", 
                                                       x = 0, y = 0, label = sprintf("Total reads: %d", 
                                                                                     sum(anndf$rc)), color = "red", hjust = 0) + 
      theme_bw() + theme(panel.grid = element_blank()) + 
      guides(fill = FALSE) + facet_wrap(~label)
    if (length(unique(statdf$Cum)) > 1) {
      p <- p + geom_line(data = statdf.summary, aes(y = Cum), 
                         color = "red", size = 0.25, linetype = "solid") + 
        scale_y_continuous(limits = c(0, NA), sec.axis = sec_axis(~. * 
                                                                    10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), 
              axis.title.y.right = element_text(color = "red"))
    }
    else {
      p <- p + ylim(c(0, NA))
    }
  }
  else {
    p <- ggplot(data = plotdf, aes(x = Cycle, y = Score)) + 
      geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", 
                                                         high = "black") + geom_line(data = statdf, aes(y = Mean), 
                                                                                     color = "#66C2A5") + geom_line(data = statdf, aes(y = Q25), 
                                                                                                                    color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      geom_line(data = statdf, aes(y = Q50), color = "#FC8D62", 
                size = 0.25) + geom_line(data = statdf, aes(y = Q75), 
                                         color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      ylab("Quality Score") + xlab("Cycle") + theme_bw() + 
      theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
      geom_text(data = anndf, aes(x = 0, label = rclabel, 
                                  y = 0), color = "red", hjust = 0) + facet_wrap(~file)
    if (length(unique(statdf$Cum)) > 1) {
      p <- p + geom_line(data = statdf, aes(y = Cum), 
                         color = "red", size = 0.25, linetype = "solid") + 
        scale_y_continuous(limits = c(0, NA), sec.axis = sec_axis(~. * 
                                                                    10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), 
              axis.title.y.right = element_text(color = "red"))
    }
    else {
      p <- p + ylim(c(0, NA))
    }
  }
  p
}

## QscoreQuantiles function -----------------------------------------------------
QscoreQuantiles = function(scores, counts, prob){
  scoreProb = 1-10^(scores/-10)
  
  testintervals1 <- cumsum(counts)
  testintervals2 <- c(0,testintervals1[-length(testintervals1)]+1)
  testintervals <- as.vector(rbind(testintervals2,testintervals1))
  
  probCounts = prob*testintervals1[length(testintervals1)]
  
  probInterval = findInterval(probCounts, testintervals, rightmost.closed = TRUE)
  
  quant = numeric(length = length(prob))
  
  inRun = probInterval %% 2 == 1
  
  quant[inRun] = scores[(probInterval[inRun]+1)/2]
  
  outRun = !inRun
  
  scoreProb1 <- scoreProb[probInterval[outRun]/2]      ## x-value to the left of the jump
  scoreProb2 <- scoreProb[probInterval[outRun]/2 + 1]  ## x-value to the right of the jump
  count1 <- testintervals1[probInterval[outRun]/2]      ## percentile to the left of the jump
  p  <- probCounts[outRun]   ## probability on the jump
  ## evaluate the line `(pl, xl) -- (pr, xr)` at `p`
  # xq[on_jump] <- (xr - xl) / (pr - pl) * (p - pl) + xl
  
  
  quant[outRun] = log10(1-((scoreProb2 - scoreProb1) * (p - count1) + scoreProb1))*-10
  
  quant
}

##getAggregateQualityScores Function -------------------------------------------
getAggregateQualityScores = function (fl, n = 5e+05, quantiles = FALSE){
  plotdf = lapply(fl[!is.na(fl)], function(x){
    srqa <- ShortRead::qa(x, n = n)
    df <- srqa[["perCycle"]]$quality
  })
  plotdf = do.call(rbind, plotdf)
  
  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf,
                              sum)
  
  plotdf.summary$ScoreProb = 1-10^(plotdf.summary$Score/-10)
  
  means <- rowsum((1-10^(plotdf.summary$Score/-10)) * plotdf.summary$Count, group = plotdf.summary$Cycle)/
    rowsum(plotdf.summary$Count,plotdf.summary$Cycle)
  
  statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                               Mean = means)
  
  if(quantiles){
    plotdf.quantiles <- by(plotdf.summary, plotdf.summary$Cycle, function(x) QscoreQuantiles(x$Score, 
                                                                                             x$Count, c(0.25,0.5,0.75)), simplify = TRUE)
    
    plotdf.quantiles = data.frame(do.call(rbind,plotdf.quantiles))
    colnames(plotdf.quantiles) = c("q25","q50","q75")
    
    statdf.summary <- cbind(statdf.summary,plotdf.quantiles)
  }
  
  return(statdf.summary)
}

##plotAggregateLengths Function ------------------------------------------------
plotAggregateLengths = function (fl, n = 5e+05, table = TRUE){
  plotdf = lapply(fl[!is.na(fl)], function(x){
    df <- ShortRead::readFastq(x)
    df <- width(df)
  })
  plotdf = do.call(c, plotdf)
  
  if(table){
    print(table(plotdf))
  }else{
    plot = ggplot()+
      geom_histogram(bins = max(plotdf)-min(plotdf)+1,aes(x = plotdf))+
      scale_x_continuous(breaks = seq.int(min(plotdf), max(plotdf)))
    print(plot)
  }
}


##plotQualityProfile3 Function -------------------------------------------------
plotQualityProfile3 = function(fnFs,fnRs,primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF, lenLimR, overlapLen = 20){
  print("Getting aggregate quality scores...")
  qpF = getAggregateQualityScores(fnFs)
  qpR = getAggregateQualityScores(fnRs)
  # qpF = qpF1
  # qpR = qpR1
  
  
  qpF = qpF[1:lenLimF,]
  qpR = qpR[1:lenLimR,]
  
  qp = data.frame(CycleF = qpF$Cycle, MeanF = qpF$Mean) %>%
    full_join(data.frame(CycleF = seqlen-qpR$Cycle+1, CycleR = qpR$Cycle, MeanR = qpR$Mean))
  
  qp = qp[order(qp$CycleF),]
  
  FRoverlap = which(!is.na(qp$MeanF) & !is.na(qp$MeanR))
  
  # multiply log10(MeanF) and log10(MeanR) by 10^6, then truncate to integers
  # perform remaining calculations with transformed values to avoid floating point number errors
  #divide all values by 10^6 for display and output at the end
  qp = qp %>%
    mutate(MeanF_ = trunc(log10(MeanF)*10^6)) %>%
    mutate(MeanR_ = trunc(log10(MeanR)*10^6))
  
  # calculate cumulative sum of log10(mean error-free rate)
  qp$csumMeanF = NA
  qp$csumMeanF[which(!is.na(qp$MeanF))] = cumsum(qp$MeanF_[which(!is.na(qp$MeanF))])
  qp$csumMeanR = NA
  qp$csumMeanR[rev(which(!is.na(qp$MeanR)))] = cumsum(qp$MeanR_[rev(which(!is.na(qp$MeanR)))])
  
  print("Starting length n subarrays...")
  # find sum of log10(mean error-free rate)_fwd at position x and log10(mean error-free rate)_rev at position x-overlapLen+1
  # score represents log10(mean error-free rate) of trimmed fwd and rev reads with a given overlapLen
  subarray_csumMean = aaply(1:(length(FRoverlap)), 1, function(n){
    csumMean = rep_len(NA, nrow(qp))
    csumMean[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] = qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] + qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]]
    return(csumMean)
  })
  subarray_csumMean = t(subarray_csumMean)
  
  subarray_csumMeanMax = adply(1:ncol(subarray_csumMean), 1, function(n){
    data.frame(
      n = n,
      score = max(subarray_csumMean[,n],na.rm = TRUE),
      CycleF = qp$CycleF[which.max(subarray_csumMean[,n])],
      CycleR = qp$CycleR[which.max(subarray_csumMean[,n])-n+1]
    )
  })
  
  # find sum of log10(mean error-free rate)_fwd and log10(mean error-free rate)_rev over overlapLen
  # score represents log10(mean error-free rate) of overlapping region with a given overlapLen
  subarray_overlapMean = aaply(1:(length(FRoverlap)), 1, function(n){
    sumMean = rep_len(NA, nrow(qp))
    sumMean[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] = qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] - qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]-n] + 
      qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]] - qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]+n]
    return(sumMean)
  })
  subarray_overlapMean = t(subarray_overlapMean)
  
  subarray_overlapMeanMax = adply(1:ncol(subarray_overlapMean), 1, function(n){
    data.frame(
      n = n,
      score = max(subarray_overlapMean[,n],na.rm = TRUE),
      CycleF = qp$CycleF[which.max(subarray_overlapMean[,n])],
      CycleR = qp$CycleR[which.max(subarray_overlapMean[,n])-n+1]
    )
  })
  
  plotData = data.frame(CycleF = qp$CycleF,
                        errorF = 10^(qp$csumMeanF/10^6), errorR = 10^(qp$csumMeanR/10^6), errorSum = 10^(subarray_csumMean[,overlapLen]/10^6),
                        overlapError = 10^(subarray_overlapMean[,overlapLen]/10^6)
  )
  plotDataMelt = reshape2::melt(plotData, id.vars = 1)
  
  plot = ggplot()+
    geom_line(data = plotDataMelt, aes(x = CycleF, y = value, group = variable, color = variable))+
    # geom_line(data = qp, color = "red", aes(x = CycleF, y = 10^(csumMeanF/10^6)))+
    # geom_line(data = qp, color = "blue", aes(x = CycleF, y = 10^(csumMeanR/10^6)))+
    # geom_line(color = "black", aes(x = qp$CycleF, y = 10^(subarray_csumMean[,overlapLen])))+
    # geom_line(color = "grey", aes(x = qp$CycleF, y = 10^(subarray_overlapMean[,overlapLen])))+
    scale_color_manual(name = "Error-free Prob.",
                       labels = c(errorF = "Fwd read", errorR = "Rev read",
                                  errorSum = "Both reads summed",
                                  overlapError = "Overlap region, ending at CycleF"),
                       values=c(errorF = "red", errorR = "blue",
                                errorSum = "black",
                                overlapError = "grey")
    )+
    scale_y_continuous(name = "cumulative error-free prob.")+
    geom_vline(xintercept = trimLeftSelect[1], color = "black") + 
    geom_vline(xintercept = seqlen-trimLeftSelect[2], color = "black") + 
    geom_vline(xintercept = truncLenSelect[1], color = "pink") +
    geom_vline(xintercept = seqlen-truncLenSelect[2], color = "lightblue") +
    geom_vline(xintercept = subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen], color = "red") + 
    geom_vline(xintercept = seqlen-subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen], color = "blue") +
    geom_text(aes(x = c(trimLeftSelect[1],seqlen-trimLeftSelect[2],
                        truncLenSelect[1],seqlen-truncLenSelect[2],
                        subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen],seqlen-subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen]
    ), 
    y = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf), 
    label = c("ampl. start", "ampl. end",
              "user truncLenF","user truncLenR",
              "recom. truncLenF", "recom. truncLenR"),
    angle = 90,
    vjust = 1, hjust = 0))
  
  print(plot)
  print(paste0("Recommended truncLenSelect: c(",
               subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen],
               ",",
               subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen],
               ")"))
  invisible(subarray_csumMeanMax)
}



##seq2fasta Function -----------------------------------------------------------
seq2fasta = function(seqNames, sequences, fileName, savePath){
  fasta = paste0(">",seqNames,"\n",sequences,"\n")
  dir.create(path = savePath, showWarnings = TRUE)
  write(fasta,file = paste0(savePath,"/",fileName,".txt"), append = FALSE)
}


##collapseNoMismatch2 Function -------------------------------------------------
collapseNoMismatch2 = function (seqtab, minOverlap = 20, orderBy = "abundance", identicalOnly = F, 
                                vec = TRUE, band = -1, verbose = T) {
  dupes <- duplicated(colnames(seqtab))
  if (any(dupes)) {
    st <- seqtab[, !dupes, drop = FALSE]
    for (i in which(dupes)) {
      sq <- colnames(seqtab)[[i]]
      st[, sq] <- st[, sq] + seqtab[, i]
    }
    seqtab <- st
  }
  if (identicalOnly) {
    return(seqtab)
  }
  
  unqs.srt <- sort(getUniques(seqtab), decreasing = TRUE)
  seqs <- names(unqs.srt)
  seqs.prefix <- substr(seqs, 1, minOverlap)
  # print(length(unique(seqs.prefix)))
  seqs.unused <- rep.int(TRUE,length(seqs))
  seqs.out <- rep.int(FALSE,length(seqs))
  collapsed <- matrix(0L, nrow = nrow(seqtab), ncol = ncol(seqtab))
  colnames(collapsed) <- colnames(seqtab)
  rownames(collapsed) <- rownames(seqtab)
  
  while (any(seqs.unused)) {
    # first query is always added
    queryNum = which.max(seqs.unused)
    print(queryNum)
    query = seqs[queryNum]
    
    query.prefix <- seqs.prefix[queryNum]
    seqs.out[queryNum] <-TRUE
    seqs.unused[queryNum] <- FALSE
    collapsed[, query] <- seqtab[, query]
    if(!any(seqs.unused)){
      break
    }
    collapseCandidates = rep.int(FALSE,length(seqs))
    collapseCandidates[seqs.unused] = grepl(
      query.prefix, seqs[seqs.unused], fixed = TRUE) | 
      sapply(seqs.prefix[seqs.unused], function(x) grepl(x, query, fixed = TRUE)
      )
    seqs.used = FALSE
    if(any(collapseCandidates)){
      seqs.used = sapply(seqs[collapseCandidates], function(ref){
        if (nwhamming(query, ref, vec = vec, band = band) == 0) {
          collapsed[, query] <<- collapsed[, query] + seqtab[, ref]
          return(TRUE)
        }else{
          return(FALSE)
        }
      })
      
      seqs.unused[collapseCandidates] = !seqs.used
    }
  }
  
  collapsed <- collapsed[,seqs.out,drop = FALSE]
  if (!is.null(orderBy)) {
    if (orderBy == "abundance") {
      collapsed <- collapsed[, order(colSums(collapsed),
                                     decreasing = TRUE), drop = FALSE]
    }
    else if (orderBy == "nsamples") {
      collapsed <- collapsed[, order(colSums(collapsed >
                                               0), decreasing = TRUE), drop = FALSE]
    }
  }
  collapsed <- collapsed[, order(colSums(collapsed), decreasing = TRUE),
                         drop = FALSE]
  if (verbose) 
    message("Output ", ncol(collapsed), " collapsed sequences out of ", 
            ncol(seqtab), " input sequences.")
  collapsed
}



##kAlign -----------------------------------------------------------------------
# Parameter documentation 
# Query parameters 

# userEmail : required to temporarily identify REST API users
# fastafileName : filepath to Kalign query FASTA file
# Waiting parameters
# waitTime = 10 : Number of seconds to wait between EBI Kalign queries. Must be >10
# Output parameters
# outDir = dirname(fastafileName) : Filepath to destination directory of .pim file. File name will default to "fastafileName_Kalign_out.pim"
# outFile = NULL : Filepath to output file, overrides outDir.

# Kalign_pim function
#REST API access to Kalign
Kalign_pim = function(userEmail,
                      fastafileName,
                      
                      waitTime = 10,
                      
                      outDir = dirname(fastafileName),
                      outFile = NULL){
  
  if(is.null(outFile)){
    outFile = paste0(outDir,"/",basename(fastafileName),"_Kalign_out.txt")
  }
  params = list(
    email = userEmail,
    stype = "dna",
    format = "fasta",
    sequence = readChar(fastafileName, file.info(fastafileName)$size))
  
  kalignRunResponse <- POST("https://www.ebi.ac.uk/Tools/services/rest/kalign/run", body = params)
  kalignJobID = content(kalignRunResponse)
  
  # check job status
  kalignJobStatusResponse = "Checking"
  while (kalignJobStatusResponse != "FINISHED") {
    kalignJobStatuspb <- txtProgressBar(min = 0, max = waitTime, style = 3)
    for(i in 1:waitTime){
      Sys.sleep(1)
      # update progress bar
      setTxtProgressBar(kalignJobStatuspb, i)
    }
    close(kalignJobStatuspb)
    kalignJobStatus = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/status/",kalignJobID))
    kalignJobStatusResponse = content(kalignJobStatus)
  }
  
  # GET response
  kalignJob_pimResponse = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/result/",kalignJobID,"/pim"))
  kalignJob_pim = content(kalignJob_pimResponse)
  
  # format response
  kalignJob_pim = read.delim(header = FALSE, skip = 6, sep = "", text = kalignJob_pim, stringsAsFactors = FALSE)
  kalignJob_pim = kalignJob_pim[,-c(1:2)]
  seqNames = ShortRead::id(readFasta(fastafileName))
  rownames(kalignJob_pim) = seqNames
  colnames(kalignJob_pim) = seqNames
  
  write.table(kalignJob_pim, 
              file = outFile,
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = FALSE
  )
  return(kalignJob_pim)
}

##Rarefaction Function ---------------------------------------------------------
# Call up and change ggrare function/Edit plot aesthetics as needed
ggrare <- function (physeq_object, step = 10, label = NULL, color = NULL, 
                    plot = TRUE, title = "default", parallel = FALSE, se = TRUE) 
{
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) {
    x <- t(x)
  }
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), 
        sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, , drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    }
    else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  }
  else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), 
                       "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "Size", 
                                                        y = ".S", group = "Sample", color = color))
  p <- p + ggplot2::labs(x = "Number of Reads", 
                         y = "ASV Richness",
                         title = title)
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels, ggplot2::aes_string(x = "x", 
                                                                   y = "y", label = label, color = color), size = 4, show.legend = FALSE,
                                hjust = 0)
  }
  p <- p + ggplot2::geom_line(linewidth = 1)
  if (se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se", 
                                                      ymax = ".S + .se", color = NULL, fill = color), 
                                  alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
} 






##DECIPHER Treeline Function ---------------------------------------------------
DECIPHER_clusterASVs = function(seqs,
                                identityCutoff = 1,
                                testBounds = FALSE, upperBound = 1, lowerBound = 1, testIncrement = 0.1,
                                ncores = NULL){
  
  # detect cores; uses all cores if ncores = NULL
  if(is.null(ncores)){
    ncores <- detectCores()
  }
  
  dna <- Biostrings::DNAStringSet(seqs)
  
  # align sequences and calculate distance matrix
  cat("Aligning sequences...\n")
  aln <- DECIPHER::AlignSeqs(dna, processors = ncores)
  cat("Calculating distance matrix...\n")
  dmat <- DECIPHER::DistanceMatrix(aln, processors = ncores)
  
  if(!testBounds){
    cat(paste0("Clustering at ",identityCutoff * 100,"% identity...\n"))
    clusters <- DECIPHER::TreeLine(
      myDistMatrix = dmat,
      method = "complete",
      cutoff = 1-identityCutoff, # use cutoff = 0.03 for a 97% OTU
      type = "clusters",
      processors = ncores)
    
  }else{
    cat(paste0("Clustering between",lowerBound*100,"-",upperBound*100,"% identity, incrementing by ",testIncrement*100,"%...\n"))
    clusters <- DECIPHER::TreeLine(
      myDistMatrix=dmat,
      method = "complete",
      cutoff = 1-seq(lowerBound, upperBound, by = testIncrement),
      type = "clusters",
      processors = ncores)
  }
  return(clusters)
}
