
# Load in seqtab files from other plates ---------------------------------------

# Plate 1 RDS
# seqtab.nochim.plate1 = readRDS(paste0("/Users/alicehong/Desktop/SK/Experiment_Data_Results/02072024_220524_16S/220524_M07519_0023_000000000-KDDMC/outputs/KDDMC_2024_02_12_seqtab.nochim.rds"))
seqtab.nochim.plate1 = readRDS(paste0("/Users/alicehong/R/16S-Analysis/data/seqtabs/KDDMC_2024_02_12_seqtab.nochim.rds"))

# Plate 2 RDS
# seqtab.nochim.plate2 = readRDS(paste0("/Users/alicehong/Desktop/SK/Experiment_Data_Results/02072024_220524_16S/220830_M07519_0027_000000000-KKLCY/outputs/KKLCY_2022_09_09_seqtab.nochim.rds"))
seqtab.nochim.plate2 = readRDS(paste0("/Users/alicehong/R/16S-Analysis/data/seqtabs/KKLCY_2022_09_09_seqtab.nochim.rds"))

# Plate 3 RDS
# seqtab.nochim.plate3 = readRDS(paste0("/Users/alicehong/Desktop/SK/Experiment_Data_Results/02072024_220524_16S/230825_M07519_0053_000000000-L6PLN/outputs/L6PLN_2023_09_06_seqtab.nochim.rds"))
seqtab.nochim.plate3 = readRDS(paste0("/Users/alicehong/R/16S-Analysis/data/seqtabs/L6PLN_2023_09_06_seqtab.nochim.rds"))

# Plate 4 RDS
# seqtab.nochim.plate4 = readRDS(paste0("/Users/alicehong/Desktop/SK/Experiment_Data_Results/02072024_220524_16S/230828_M07519_0054_000000000-L6PM5/outputs/L6PM5_2023_09_06_seqtab.nochim.rds"))
seqtab.nochim.plate4 = readRDS(paste0("/Users/alicehong/R/16S-Analysis/data/seqtabs/L6PM5_2023_09_06_seqtab.nochim.rds"))

# Plate 5 RDS
# seqtab.nochim.plate5 = readRDS(paste0("/Users/alicehong/Desktop/SK/Experiment_Data_Results/WON26772.20240306/240229_M07519_0079_000000000-LFGPK/outputs/WON_2024_03_06_seqtab.nochim.rds"))
seqtab.nochim.plate5 = readRDS(paste0("/Users/alicehong/R/16S-Analysis/data/seqtabs/WON_2024_03_06_seqtab.nochim.rds"))

# Append plateX to control samples in each RDS file prior to merging to avoid errors
rds.list = list()

n = 5 # Set to number of plates you have
for (i in 1:n){
  rds.list[[(paste0("seqtab.nochim.plate",i))]] <- eval(parse(text = paste0("seqtab.nochim.plate",i)))
  rownames(rds.list[[i]]) <- gsub("control|Control", paste0("control_plate",i), rownames(rds.list[[i]]))
  rownames(rds.list[[i]]) <- gsub("extraction", paste0("extraction_plate",i), rownames(rds.list[[i]]))
  rownames(rds.list[[i]]) <- gsub("inocula", paste0("inocula_plate",i), rownames(rds.list[[i]]))
}

# If using more than two RDS files, need to perform merge and check again for bimeras
if (n != 1){
  seqtab.merge = mergeSequenceTables(tables = rds.list)
  seqtab.merge = removeBimeraDenovo(seqtab.merge, method="consensus", multithread=TRUE)
} else {
  seqtab.merge = rds.list[[1]]
}

# Diagnostics to check for shortest sequence
table(nchar(colnames(seqtab.merge))) # Check frequency of n-mer consensus sequences
hist((nchar(colnames(seqtab.merge))))