###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# MAIN R SCRIPT FOR: a case study using FIBOS, comparing packing densities 
# between experimental protein structures and AlphaFold predictions.
#
# AUTHORS: 
#
# Herson Soares, João Romanelli, Carlos Silveira
# Federal University of Itajubá (Unifei), Itabira, MG, Brazil.
#
# Patrick Fleming,
# Johns Hopkins University (JHU), Baltimore, MD, USA.
#
# CONTACT: carlos.silveira@unifei.edu.br
#
# FIBOS install on R:
# library(devtools)
# install_github("https://github.com/insilico-unifei/FIBOS-R.git") 

# LIBRARIES
library(FIBOS)
library(tidyverse)
library(fs)
library(httr)
library(bio3d)


# SOURCE OF AUXILIARY FUNCTIONS
source("FIBOS-study-case-fun.R")


# DATA INPUT: 
# If dataset is missing, download from the 1st sheet ("protein names") of the link below in 
# tsv format with the name: Fleming_base.tsv, and place it in a folder called "data".
# https://docs.google.com/spreadsheets/d/1TmBO15DQh9Sv9UXhXyMt3emQTztAk9YQt8kObdb5S3g/edit?usp=sharing
if(0){
  folder <- "data" # data input folder
  filename <- "Fleming_base.tsv" # dataset
  file <- path(folder,filename)
  db <- read_tsv(file) 
}

# DOWNLOAD EXPERIMENTAL PDB AND RETURN PDB PATHs
if(0){
  folder <- "data_exp_pdb" # experimental PDB output folder
  if (!dir.exists(folder)) dir_create(folder)
  pdb.ids = db$PDB_ID 
  pdb.path <- pdb.ids |> get.pdb(path = folder) |> as.list() 
}

# GET METADATA: UNIPROT ID AND AMINO ACID SEQUENCE OF EXPERIMENTAL MODELS
# AUXILIARY FUNCTIONS:
# get_pdb_metadata: get pdb metadata from source url and par
# get_pdb_uniprot_id: get uniprot id from pdb metadata
# get_pdb_aa_seq: get aa sequence from pdb metadata
if(0){
  source <- "https://data.rcsb.org/rest/v1/core/polymer_entity/"
  par <- "/1"
  #pdb.meta <- pdb.ids |> map(\(x) get_pdb_metadata(x, source, par))
  pdb.uniprot.id <- pdb.meta |> map(\(x) get_pdb_uniprot_id(x))
  pdb.aa.seq <- pdb.meta |> map(\(x) get_pdb_aa_seq(x))
}

# FORM CSM-AF NAMES FROM UNIPROT IDs
if(0){
  # for access to the AF website
  prefix <- "AF-"
  suffix <- "-F1-model_v4"
  csm.af.ids <- pdb.uniprot.id |>  map(\(x) paste0(prefix, x, suffix)) |> unlist()
  
  # for access to the RCSB website
  prefix <- "AF_AF"
  suffix <- "F1"
  csm.rcsb.ids <- pdb.uniprot.id |>  map(\(x) paste0(prefix, x, suffix)) |> unlist()
}

# DOWNLOAD PREDICTED AF AND RETURN AF PATHs
# AUXILIARY FUNCTIONS:
# get_csm: download AF in PDB format and return AF paths
if(0){
  folder <- "data_csm" # AF output folder
  if (!dir.exists(folder)) dir_create(folder)
  source <- "https://alphafold.ebi.ac.uk/files/"
  csm.path <- csm.af.ids |> map(\(x) get_csm(x, source = source, path = folder)) |> as.list()
  
  # NA here means a failure to get AF models
  csm.ok <- !(csm.path |> is.na()) 
  
  # get only success AF in PDB format with bio3d
  csm.bio3d <- csm.path[csm.ok] |> map(\(x) bio3d::read.pdb(x)) 
}

# GET AF METADATA AND AMINO ACID SEQUENCE OF AF MODELS
if(0){
  source <- "https://data.rcsb.org/rest/v1/core/polymer_entity/"
  par <- "/1"
  csm.meta <- csm.rcsb.ids[csm.ok] |> map(\(x) get_pdb_metadata(x, source, par)) 
  
  # NA here means a failure to get AF metadata
  csm.meta.ok <- !(csm.meta |> is.na()) 
  
  # get aa seq from only success AF metadata
  csm.aa.seq <- csm.meta[csm.meta.ok] |> map(\(x) get_pdb_aa_seq(x))
}

# CREATE THE TABLES FOR COMPARATIVE CALCULATIONS OF THE OSP
if(0){
  # main table with existing coordinates and metadata for AF models
  pdb.csm.tab <- tibble(PDB.ids = pdb.ids[csm.ok][csm.meta.ok], 
                        PDB.path = pdb.path[csm.ok][csm.meta.ok],
                        UniProt = pdb.uniprot.id[csm.ok][csm.meta.ok] ,
                        CSM.RCSB.ids = csm.rcsb.ids[csm.ok][csm.meta.ok],
                        CSM.path = csm.path[csm.ok][csm.meta.ok],
                        aaPDB = pdb.aa.seq[csm.ok][csm.meta.ok] |> nchar(),  
                        aaCSM = csm.aa.seq |> nchar()) 
  
  # complement with more information...
  pdb.csm.tab <- pdb.csm.tab |> mutate(aaDIFF = aaCSM-aaPDB) |> # get diff in aa seq sizes
                 mutate(pLDDT_global = csm.bio3d[csm.meta.ok] |> 
                 map(\(x) round(mean(x$atom$b),2)) |> unlist()) # calculate global pLDDT metric
  
  # filter main table in a work table with difference between AF and Experimental less than 5 aa in size
  pdb.csm.tab.work <- pdb.csm.tab |> filter(abs(aaDIFF) < 5)
  
  # create study case output table
  if(0) write_csv(pdb.csm.tab.work |> mutate(across(where(is.list), as.character)), 
                  file = "data/study-case-v1.csv")
}

#FOR CONVENIENCE, CREATE FOLDER FOR AF FILE MODELS WITH SAME NAME OF FILE EXPERIMENTAL MODELS
if(0){
  folder.from <- "data_exp_pdb" # folder of experimental models
  folder.to <- "data_csm_pdb" # folder of AF models
  if (!dir.exists(folder.to)) dir_create(folder.to)
  
  # adapt the paths in work table
  pdb.csm.tab.work <- pdb.csm.tab.work |>
    mutate(CSM.path.new = PDB.path) |> 
    mutate(CSM.path.new = gsub(folder.from, folder.to, CSM.path.new)) |> 
    relocate(CSM.path.new, .after = CSM.path)
  
  # copy AF downloaded models to new working folder
  file_copy(unlist(pdb.csm.tab.work$CSM.path), folder.to, overwrite = T)
  
  
  # rename AF file models with same name of file experimental models
  folder.from <- "data_csm_pdb"
  gsub(folder.from, folder.to, pdb.csm.tab.work$CSM.path) |> 
    map2(pdb.csm.tab.work$CSM.path.new, \(x,y) file.rename(x,y))
  
}

# PREPARE PATH FOR SRF FILES THAT WILL BE USED TO CALCULATE OSP
if(0){
  path <- "data_exp_pdb"
  pdb.csm.tab.work <- pdb.csm.tab.work |> mutate(SRF.path = gsub("/", "/prot_", PDB.path)) |> 
                           mutate(SRF.path = gsub("\\.pdb",".srf", SRF.path)) |>
                           mutate(SRF.path = gsub(path, "fibos_files", SRF.path)) |> 
                           relocate(SRF.path, .after = CSM.path.new)

}


if(0){
  
  pdb.exp.fibos = pdb.csm.tab.work$PDB.path %>% map(occluded_surface, method = "FIBOS")

  pdb.exp.osp <- pdb.csm.tab.work$SRF.path |> map(osp)

  file.rename("fibos_files","fibos_files_exp_pdb")
  
}


