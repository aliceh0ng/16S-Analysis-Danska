
# Required Packages ------------------------------------------------------------
packages <- c("BiocManager", 
              "devtools", 
              "dada2", 
              "plyr", 
              "stringr", 
              "reshape2", 
              "ggplot2", 
              "doParallel", 
              "ggalluvial", 
              "R.utils", 
              "dplyr", 
              "vegan", 
              "taxize", 
              "ShortRead", 
              "Rfast", 
              "httr", 
              "scales", 
              "parallel",
              "phyloseq",
              "tidyr",
              "tibble",
              "remotes",
              "ranacapa",
              "data.table",
              "readxl",
              "xfun",
              "purrr",
              "ggpubr",
              "readr")


# Install regular packages if needed
# install.packages(setdiff(packages, rownames(installed.packages())))  

# # If needed, install BiocManager
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# #If needed, install Phyloseq
# BiocManager::install("phyloseq")
# 
# #If needed, install ggrare
# remotes::install_github("gauravsk/ranacapa")
#
# In case libraries are not being loaded properly (Windows), can manually specify library location
# and also repo to download from 
#
# BiocManager::repositories()
# getOption("repos")
# BiocManager::install("Biostrings", lib = "D:/Documents/R Libraries", repo = "https://bioconductor.org/packages/3.15/bioc")

# Feb 12 2024: dada2 and ShortRead "not available for this version of R"
# install.packages("devtools")
# library("devtools")
# devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

# or use:
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("dada2")

# output: "The downloaded binary packages are in /var/folders/b4/6c_5lqt14pb2cbp0vq45mdth0000gn/T//RtmpGyzXZE/downloaded_packages"

# Load packages ----------------------------------------------------------------
lapply(packages, require, character.only = TRUE)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "./16S_ribosomal_RNA/", sep = .Platform$path.sep))
library(tidyverse)
library(readxl) 
library(ggsci)
library(ggplot2)
library(gridExtra)
library(reshape2) 
library(reshape) 
library(tidyverse)
library(grid)
library(gridtext)
library("DECIPHER")
library(seqinr)
library(stringr)

# Set File Paths ---------------------------------------------------------------

# Input path for files
 path = "/Users/alicehong/Desktop/SK/Experiment_Data_Results/WON26772.20240306/240229_M07519_0079_000000000-LFGPK"

# Path that will house demultiplexed sequencing files
 demultiplexed_path <- "/Users/alicehong/Desktop/SK/Experiment_Data_Results/WON26772.20240306/240229_M07519_0079_000000000-LFGPK"
 