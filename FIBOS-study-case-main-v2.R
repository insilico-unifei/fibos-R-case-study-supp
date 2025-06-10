###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# MAIN R SCRIPT FOR: a case study using FIBOS, comparing packing densities 
# between experimental protein structures (EXP) and AlphaFold predictions (CSM).
#
# AUTHORS: 
#
# Herson Soares, João Romanelli, Carlos Silveira
# Federal University of Itajubá (Unifei), Itabira, MG, Brazil.
#
# Patrick Fleming,
# Johns Hopkins University (JHU), Baltimore, MD, USA.
#
# CONTACTS: hersinsoares@gmail.com, carlos.silveira@unifei.edu.br
#
# FIBOS install on R:
# install.packages("fibos")
# See more here: https://github.com/insilico-unifei/fibos-R  

#LIBRARIES
library(fibos)
library(tidyverse)
library(fs)
library(httr)
library(bio3d)
library(furrr)
library(tictoc)
library(ggpubr)
library(tictoc)
library(effsize)
library(rPDBapi)
library(FactoMineR)
library(factoextra)
library(colorspace)
library(corrplot)
library(ggcorrplot)
library(Rpdb)
library(Biostrings)
library(msa)
library(scales)
library(moments)
library(rsample)
library(ggwordcloud)
library(Hmisc)
library(RColorBrewer)
library(ggseqlogo)
library(conflicted)
library(jsonlite)
library(grid)
library(ggrepel)
library(ggforce)
library(patchwork)
library(magick)

conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("desc", "dplyr")

# SOURCE OF AUXILIARY FUNCTIONS
source("FIBOS-case-study-expanded-more-fun.R")


# CHANGE CONDITIONALS TO 1 TO EXECUTE STEP BY STEP


# SETS UP VENV AND INSTALLS DEPENDENCIES FOR PYTHON SCRIPTS  
# AUXILIARY FUNCTIONS:  
# prepare_python_venv: creates Python venv in specified directory and installs  
# requirements
if(0){
  folder.venv <- ".venv"
  prepair_python_venv(folder.venv)
}

# CALLS GET_RAW_PDBS_PY WITH SPECIFIED CSV TO DOWNLOAD PDB IDS  
# AUXILIARY FUNCTIONS:  
# get_raw_pdbs_py: invokes Python script within venv to download raw PDB IDs
if(0){
  folder <- "data" # data reference folder
  if (!dir.exists(folder)) dir_create(folder)
  filename <- "pdb_raw_base.csv" # raw dataset pdb ids only
  file <- path(folder,filename)
  get_raw_pdbs_py(file)
}


# READS CSV OF RAW PDB IDS
if(0){ 
  folder <- "data"
  filename <- "pdb_raw_base.csv"  
  file <- path(folder, filename)
  db.raw <- read_csv(file, col_names = "PDB_ID") 
  db.raw <- db.raw |> separate(PDB_ID, c("PDB_ID", "X")) |> select(PDB_ID)
}


# DOWNLOADS EXPERIMENTAL PDBS FROM IDS IN DB.RAW AND RETURNS PDB PATHS
if(0){
  folder <- "data_exp_pdb" # experimental PDB output folder
  if (!dir.exists(folder)) dir_create(folder)
  # pdb.ids.raw <- db.raw$PDB_ID
  # pdb.ids.raw <- c("1G3P",tail(db.raw$PDB_ID, 20))
  pdb.ids.raw <- db.raw$PDB_ID
  pdb.return <- pdb.ids.raw |> get.pdb(path = folder) 
  pdb.path.raw <- pdb.return |> names()
  if(is.null(pdb.path.raw)) pdb.path.raw <- pdb.return
  
  # some PDBs fail to download
  pdb.get.ok <- !(pdb.return == "1")
  pdb.get.no.ok.size <- sum(!pdb.get.ok, na.rm = T)
  if (pdb.get.no.ok.size > 0){
    print(paste("PDBids with download problems:", pdb.get.no.ok.size))
    print(pdb.ids.raw[!pdb.get.ok])
  }
}

# FILTERS PDBS WITHOUT MISSING RESIDUES IN THE MIDDLE OF THE CHAIN  
# AUXILIARY FUNCTIONS:  
# verify_missing_residues: returns a table indicating possible missing residues in PDB  
# get_pdb_info: gets some experimental PDB parameters like title, resolution, r_factors
if(0){
  pdb.bio3d <- pdb.path.raw[pdb.get.ok] |> map(\(x) bio3d::read.pdb(x, verbose = FALSE))
  pdb.res.miss <- pdb.bio3d |> map2_dfr(pdb.ids.raw[pdb.get.ok], \(x,y) verify_missing_residues(x,y))
  pdb.res.ok <- pdb.res.miss$MISS_POINT == ""
  pdb.ids <- pdb.ids.raw[pdb.get.ok][pdb.res.ok]
  pdb.no.res.miss <- pdb.res.miss |> filter(pdb.res.ok) |> select(!MISS_POINT)
  db <- pdb.ids |> map_dfr(\(x) get_pdb_info(x)) |> left_join(pdb.no.res.miss, by = "PDB_ID")
  db$PDB.path <- pdb.path.raw[pdb.get.ok][pdb.res.ok]
  pdb.bio3d <- pdb.bio3d[pdb.res.ok]
}

# GETS METADATA: UNIPROT ID AND AMINO ACID SEQUENCE OF EXPERIMENTAL MODELS  
# AUXILIARY FUNCTIONS:  
# get_pdb_metadata: gets PDB metadata from source URL and PAR  
# get_pdb_uniprot_id: gets UniProt ID from PDB metadata
if(0){
  source <- "https://data.rcsb.org/rest/v1/core/polymer_entity/"
  par <- "/1"
  pdb.meta <- db$PDB_ID |> map(\(x) get_pdb_metadata(x, source, par))
  pdb.uniprot.id.raw <- pdb.meta |> map(\(x) get_pdb_uniprot_id(x))
  db$uniprot.id <- pdb.uniprot.id.raw |> unlist()
  uniprot.ok <- !(pdb.uniprot.id.raw |> is.na())
  db <- db |> filter(uniprot.ok)
  pdb.bio3d <- pdb.bio3d[uniprot.ok]
}

# FORMS CSM-AF NAMES FROM UNIPROT IDS
if(0){
  # for access to the AF website
  prefix <- "AF-"
  suffix <- "-F1-model_v4"
  db$csm.af.id <- db$uniprot.id |>  map(\(x) paste0(prefix, x, suffix)) |> unlist()
  # for access to the RCSB website
  prefix <- "AF_AF"
  suffix <- "F1"
  db$csm.rcsb.id <- db$uniprot.id |>  map(\(x) paste0(prefix, x, suffix)) |> unlist()
}

# DOWNLOADS PREDICTED AF AND RETURNS AF PATHS  
# AUXILIARY FUNCTIONS:  
# get_csm: downloads AF in PDB format and returns AF paths
if(0){
  folder <- "data_csm" # AF output folder
  if (!dir.exists(folder)) dir_create(folder)
  source <- "https://alphafold.ebi.ac.uk/files/"
  db$CSM.path <- db$csm.af.id |> map(\(x) get_csm(x, source = source, path = folder)) |> unlist()
  # NA here means a failure to get AF models
  csm.ok <- !(db$CSM.path |> is.na()) 
  db <- db |> filter(csm.ok)
  pdb.bio3d <- pdb.bio3d[csm.ok]
}

# GETS AF METADATA AND AMINO ACID SEQUENCE OF AF MODELS
if(0){
  source <- "https://data.rcsb.org/rest/v1/core/polymer_entity/"
  par <- "/1"
  csm.meta <- db$csm.rcsb.id |> map(\(x) get_pdb_metadata(x, source, par)) 
  # NA here means a failure to get AF metadata
  csm.meta.ok <- !(csm.meta |> is.na()) 
  csm.bio3d <- db$CSM.path |> map(\(x) bio3d::read.pdb(x))
  csm.res.miss <- csm.bio3d |> map2_dfr(db$PDB_ID, \(x,y) verify_missing_residues(x,y))
  #sum(abs(csm.res.miss$n_SEQSTR-csm.res.miss$n_SEQRES)>0)
  csm.res.ok <- csm.res.miss$MISS_POINT == ""
  if (sum(!csm.res.ok)>0) print("Warning! There is missing residues in AF")
  db <- db |> filter(csm.res.ok)
  pdb.bio3d <- pdb.bio3d[csm.res.ok]
  csm.bio3d <- csm.bio3d[csm.res.ok]
  db$n_SEQAF <- csm.res.miss$n_SEQSTR
  db <- db |> relocate(n_SEQAF, .after = "n_SEQRES")
}

# CALCULATES SEQDIFF, CREATES CLEANED EXP AND CSM PDB DIRECTORIES AND FILES  
# READS CLEANED PDB FILES WITH BIO3D AND ADDS A PDB_ID ATTRIBUTE  
# AUXILIARY FUNCTIONS  
# clean_pdb2: cleans PDB atom records by filtering out hydrogens, heteroatoms, and  
# unwanted residues.
if(0){
  db <- db |> mutate(n_SEQDIFF = abs(n_SEQAF-n_SEQSTR)) |> relocate(n_SEQDIFF, .after = "n_SEQAF")
  folder_cls <- "data_exp_pdb_cls"
  db.work <- db |> mutate(PDB.cls.path = path(folder_cls, PDB_ID, ext = "pdb")) |> 
             relocate(PDB.cls.path, .after = PDB.path)
  if (!dir.exists(folder_cls)) dir_create(folder_cls)
  db.work$PDB.path |> walk2(db.work$PDB.cls.path, \(x,y) clean_pdb2(x, y, i = 1))
  folder_cls <- "data_csm_pdb_cls"
  db.work <- db.work |> mutate(CSM.cls.path = path(folder_cls, PDB_ID, ext = "pdb")) |> 
             relocate(CSM.cls.path, .after = CSM.path)
  if (!dir.exists(folder_cls)) dir_create(folder_cls)
  db.work$CSM.path |> walk2(db.work$CSM.cls.path, \(x,y) clean_pdb2(x, y, i = 1))
  pdb.bio3d.cls.work.preali <- db.work$PDB.cls.path |> map(\(x) bio3d::read.pdb(x, verbose = FALSE))
  pdb.bio3d.cls.work.preali <- pdb.bio3d.cls.work.preali |> map2(db.work$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
  csm.bio3d.cls.work.preali <- db.work$CSM.cls.path |> map(\(x) bio3d::read.pdb(x, verbose = FALSE))
  csm.bio3d.cls.work.preali <- csm.bio3d.cls.work.preali |> map2(db.work$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
}

# CREATES EXP AND CSM ALIGNMENT FOLDERS AND ALIGNS PDB PAIRS WITH LOGGING  
# FILTERS SUCCESSFUL ALIGNMENTS, UPDATES FILE PATHS, AND READS ALIGNED PDBS  
# AUXILIARY FUNCTIONS:  
# align_and_save_pdb2: aligns two PDB sequences, identifies consensus region, and  
# trims atoms accordingly
if(0){
  folder_exp <- "data_exp_pdb_cls_ali"
  if (!dir.exists(folder_exp)) dir_create(folder_exp)
  folder_csm <- "data_csm_pdb_cls_ali"
  if (!dir.exists(folder_csm)) dir_create(folder_csm)
  aligned.ok <- pdb.bio3d.cls.work.preali |> 
                map2_lgl(csm.bio3d.cls.work.preali, 
                         \(x,y) align_and_save_pdb2(x, y, folder_exp, folder_csm, log = T, save = T))
  print(sum(aligned.ok))
  db.work <- db.work |> filter(aligned.ok)
  pdb.bio3d.cls.work <- pdb.bio3d.cls.work.preali[aligned.ok]
  csm.bio3d.cls.work <- csm.bio3d.cls.work.preali[aligned.ok]
  db.work <- db.work |> mutate(PDB.cls.ali.path = path(folder_exp, PDB_ID, ext = "pdb")) |> 
             relocate(PDB.cls.ali.path, .after = PDB.cls.path)
  db.work <- db.work |> mutate(CSM.cls.ali.path = path(folder_csm, PDB_ID, ext = "pdb")) |> 
             relocate(CSM.cls.ali.path, .after = CSM.cls.path)
  pdb.bio3d.cls.ali.work <- db.work$PDB.cls.ali.path |> map(\(x) bio3d::read.pdb(x, verbose = FALSE))
  pdb.bio3d.cls.ali.work <- pdb.bio3d.cls.ali.work |> map2(db.work$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
  csm.bio3d.cls.ali.work <- db.work$CSM.cls.ali.path |> map(\(x) bio3d::read.pdb(x, verbose = FALSE))
  csm.bio3d.cls.ali.work <- csm.bio3d.cls.ali.work |> map2(db.work$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
}

# COMPUTES SSE FOR SELECTED PDBS AND CALCULATES GLOBAL PLDDT
if(0){
  pdb.sse <- pdb.bio3d[db$PDB_ID %in% db.work$PDB_ID]|> map2_dfr(db.work$PDB.path, \(x,y) get_pdb_sse(x,y)) 
  db.work <- db.work |> left_join(pdb.sse, by = "PDB_ID") 
  db.work <- db.work |> relocate(HELIX, .after = "R_FACTOR_W") |> 
             relocate(STRAND, .after = "HELIX")
  db.work <- db.work |> mutate(pLDDT_global = csm.bio3d.cls.ali.work |> 
                        map(\(x) round(mean(x$atom$b),2)) |> unlist()) |> 
                        relocate(pLDDT_global, .after = "n_SEQDIFF")
}


# PREPARES PATH FOR SRF FILES THAT WILL BE USED TO CALCULATE OSP
if(0){
  folder_fibos <- "fibos_files"
  db.work <- db.work |> mutate(SRF.path = path(folder_fibos, paste0("prot_", PDB_ID), ext = "srf")) |> 
                        relocate(SRF.path, .after = CSM.cls.ali.path)
  folder_fibos <- "fibos_files_exp_pdb"
  db.work <- db.work |> mutate(PDB.srf.path = path(folder_fibos, paste0("prot_", PDB_ID), ext = "srf")) |> 
             relocate(PDB.srf.path, .after = PDB.cls.ali.path)
  folder_fibos <- "fibos_files_csm_pdb"
  db.work <- db.work |> mutate(CSM.srf.path = path(folder_fibos, paste0("prot_", PDB_ID), ext = "srf")) |> 
             relocate(CSM.srf.path, .after = CSM.cls.ali.path)
}

# SAVES DB.WORK TO CSV
if(0){
  folder <- "data"
  file <- "db-work-ali-v2.csv"
  db.work |> write_csv(path(folder,file))
}

# FILTERS OUT SOME PROBLEMATIC PDBS
if(0){
  exclude <- c("1FCX", "1BKR", "5EHA", "7EV5", "7NMQ")
  fibos.ok <- !(str_detect(db.work$SRF.path, paste(exclude, collapse = "|")))
  db.work <- db.work |> filter(fibos.ok)
  pdb.bio3d.cls.ali.work <- pdb.bio3d.cls.ali.work[fibos.ok]
  csm.bio3d.cls.ali.work <- csm.bio3d.cls.ali.work[fibos.ok]
}

# CALCULATES OCCLUDED SURFACE AT ATOM AND RESIDUE LEVEL FOR EXPERIMENTAL PDB  
# WARNING: THIS MAY TAKE A WHILE DEPENDING ON YOUR HARDWARE CONFIGURATION AND 
# THE NUMBER OF CORES AVAILABLE

if(0){

  fibos.exp.ok <- db.work$SRF.path |> file_exists()
  # calculate OS at atom level in parallel with max cores available according to PDB_ids size
  tic()
  default_mccores <- getOption("mc.cores")
  ideal_mccores <- min(parallel::detectCores(), sum(!fibos.exp.ok))
  if (ideal_mccores > 0) options(mc.cores = ideal_mccores)
  if (ideal_mccores > 1) plan(multisession, workers = ideal_mccores)
  aux <- db.work$PDB.cls.ali.path[!fibos.exp.ok] |> future_map(\(x) occluded_surface (x, method = "FIBOS"), 
                                        .options = furrr_options(seed = 123)) #antes no pdb.exp.fibos
  if (ideal_mccores > 0) options(mc.cores = default_mccores)
  toc()
  fibos.exp.ok <- db.work$SRF.path |> file_exists()
}

# VERIFIES EXP FIBOS RESULTS AND STOPS ON ERRORS  
# IF ALL ARE SUCCESSFUL, PRINTS MESSAGE AND RENAMES FIBOS_FILES FOLDER
if(0){
  if (exists("aux")) rm(aux)
  pdb.exp.probs <- db.work$PDB_ID[!fibos.exp.ok]
  if(length(pdb.exp.probs)>0){
    print("WARNING: something wrong in EXP FIBOS calculation for:")
    print(pdb.exp.probs)
    print("Eliminate them from db.work and run again")
    stop("EXP FIBOS calculation error", call. = FALSE)
  } else {
    print(paste("CONGRATULATIONS: all", length(fibos.exp.ok), "EXP FIBOS calculation done!"))
    file_in <- "fibos_files"
    file_out <- "fibos_files_exp_pdb"
    print(paste("Renaming", file_in, "to", file_out))
    file.rename(file_in, file_out)
  }
  # db.work <- db.work |> filter(fibos.exp.ok)
  # pdb.bio3d.cls.ali.work <- pdb.bio3d.cls.ali.work[fibos.exp.ok]
  # csm.bio3d.cls.ali.work <- csm.bio3d.cls.ali.work[fibos.exp.ok]
}

# CALCULATES OCCLUDED SURFACE AT ATOM AND RESIDUE LEVEL FOR AF MODELS  
# WARNING: THIS MAY TAKE A WHILE DEPENDING ON YOUR HARDWARE CONFIGURATION AND 
# THE NUMBER OF CORES AVAILABLE
if(0){
  
  fibos.csm.ok <- db.work$SRF.path |> file_exists()
  # calculate OS at atom level in parallel with max cores available according to PDB_ids size
  tictoc::tic()
  default_mccores <- getOption("mc.cores")
  ideal_mccores <- min(parallel::detectCores(), sum(!fibos.csm.ok))
  if (ideal_mccores > 0) options(mc.cores = ideal_mccores)
  if (ideal_mccores > 1) plan(multisession, workers = ideal_mccores)
  aux <- db.work$CSM.cls.ali.path[!fibos.csm.ok] |> future_map(\(x) occluded_surface(x, method = "FIBOS"), 
                                                               .options = furrr_options(seed = 123))
  if (ideal_mccores > 0) options(mc.cores = default_mccores)
  tictoc::toc()
  fibos.csm.ok <- db.work$SRF.path |> file_exists()
}

# VERIFIES CSM FIBOS RESULTS AND STOPS ON ERRORS  
# IF ALL ARE SUCCESSFUL, PRINTS MESSAGE AND RENAMES fibos_files FOLDER
if(0){
  if (exists("aux")) rm(aux)
  pdb.csm.probs <- db.work$PDB_ID[!fibos.csm.ok]
  if(length(pdb.csm.probs)>0){
    print("WARNING: something wrong in CSM FIBOS calculation for:")
    print(pdb.csm.probs)
    print("Eliminate them from db.work and run again")
    stop("CSM FIBOS calculation error", call. = FALSE)
  } else {
    print(paste("CONGRATULATIONS: all", length(fibos.csm.ok), "CSM FIBOS calculation done!"))
    file_in <- "fibos_files"
    file_out <- "fibos_files_csm_pdb"
    print(paste("Renaming", file_in, "to", file_out))
    file.rename(file_in, file_out)
  }
}

# LOADS PROTEIN FIBOS FOR EXP AND CSM STRUCTURES AND CALCULATES RESIDUE-LEVEL  
# OSP METRICS  
# AUXILIARY FUNCTIONS:  
# read_prot2: reads SRF files (OS by atom)
if(0){
  pdb.exp.fibos <- db.work$PDB.srf.path |> map(\(x) read_prot2(x))
  names(pdb.exp.fibos) <- paste0(db.work$PDB_ID,"_exp")
  #calculate the EXP OSP metric at the residue level
  pdb.exp.osp <- db.work$PDB.srf.path |> map(\(x) osp(x))
  names(pdb.exp.osp) <- paste0(db.work$PDB_ID,"_exp")
  pdb.csm.fibos <- db.work$CSM.srf.path |> map(\(x) read_prot2(x))
  names(pdb.csm.fibos) <- paste0(db.work$PDB_ID,"_csm")
  #calculate the CSM OSP metric at the residue level
  pdb.csm.osp <- db.work$CSM.srf.path |> map(\(x) osp(x))
  names(pdb.csm.osp) <- paste0(db.work$PDB_ID,"_csm")
}

# CHECKS INCONSISTENCIES INVOLVING PLDDT, EXP AND CSM SIZES, AND RESIDUE MATCHING
if(0){
  f1 <- db.work$pLDDT_global > 0
  x = pdb.exp.osp |> map_int(\(x) dim(x)[1])
  y = pdb.csm.osp |> map_int(\(x) dim(x)[1])
  f2 <- ((x - y) == 0)
  z = pdb.exp.osp |> map2_int(pdb.csm.osp, \(x,y) sum(x$Resname != y$Resname ))
  f3 <- z==0
  final.ok <- f1 & f2 & f3
  print(sum(!final.ok))
  db.eff <- db.work |> filter(final.ok) 
  exp.osp <- pdb.exp.osp[final.ok]
  csm.osp <- pdb.csm.osp[final.ok]
  exp.bio3d <- pdb.bio3d.cls.ali.work[final.ok]
  csm.bio3d <- csm.bio3d.cls.ali.work[final.ok]
}

# COMPUTES MEAN B-FACTORS FOR EXP AND CSM STRUCTURES INTO LISTS. IN CSM, B-FACTORS  
# REPRESENT PLDDT  
# AUXILIARY FUNCTIONS:  
# get_mean_b_factor: calculates mean B-factor per residue from atom table
if(0){
  exp.b <- exp.bio3d |> map(\(x) get_mean_b_factor(x$atom))
  names(exp.b) <- paste0(db.eff$PDB_ID)
  exp.b.uni <- exp.b |> bind_rows(.id = "PDB_ID") |> 
               rename(Resnum.exp = Resnum, Chain.exp = Chain)
  
  csm.b <- csm.bio3d |> map(\(x) get_mean_b_factor(x$atom))
  names(csm.b) <- paste0(db.eff$PDB_ID)
  csm.b.uni <- csm.b |> bind_rows(.id = "PDB_ID") |> 
    rename(Resnum.csm = Resnum, Chain.csm = Chain, pLDDT = b)
}

# PREPARES AND ADJUSTS EXP AND CSM OSP TABLES AND MERGES INTO UNIFIED OSP TABLE  
# AUXILIARY FUNCTIONS:  
# adjust_osp_table: recovers residue ID information from original PDB
if(0){
  exp.osp <- exp.osp |> map2(db.eff$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
  exp.osp <- exp.osp |> map2(exp.bio3d, \(x, y) adjust_osp_table(x, y, resnum.original=TRUE)) |> 
             map(~ mutate(.x, Resid = row_number()) |> relocate(Resid, .before = 1))
  csm.osp <- csm.osp |> map2(db.eff$PDB_ID, \(x,y) add_attr(x, y, "PDB"))
  csm.osp <- csm.osp |> map2(csm.bio3d, \(x, y) adjust_osp_table(x, y, resnum.original=TRUE)) |> 
             map(~ mutate(.x, Resid = row_number()) |> relocate(Resid, .before = 1))
  exp.osp.uni <- exp.osp |> bind_rows(.id = "ID") |> separate(ID, c("PDB_ID", "TYPE")) 
  csm.osp.uni <- csm.osp |> bind_rows(.id = "ID") |> separate(ID, c("PDB_ID", "TYPE")) 
  osp.uni <- exp.osp.uni |> left_join(csm.osp.uni, by = c("PDB_ID", "Resid", "Resname"),
                                      suffix = c(".exp", ".csm"))
  osp.uni <- osp.uni |> mutate(OSP.diff = OSP.csm - OSP.exp)
  osp.uni.order <- osp.uni |> arrange(desc(abs(OSP.diff))) |> filter(abs(OSP.diff) > 0.2)
  osp.uni <- osp.uni |> left_join(exp.b.uni, by = c("PDB_ID", "Resnum.exp", "Resname", "Chain.exp"))
  osp.uni <- osp.uni |> left_join(csm.b.uni, by = c("PDB_ID", "Resnum.csm", "Resname", "Chain.csm")) |> 
             relocate(pLDDT, .after = "b")
  if(0){
    folder <- "data"
    file <- "pdb-exp-csm-osp-diff-v4.csv"
    file <- "osp.exp-csm-uni-v5.csv"
    write_csv(osp.uni, path(folder, file))
  }
}

# COMPUTES TORSION ANGLES FOR EACH UNIQUE PDB IN EXP AND CSM  
# AUXILIARY FUNCTIONS:  
# get_torsions_py: gets side chain chi angles via Python script
if(0){
  folder_in <- "data_exp_pdb_cls_ali"
  folder_out <- "data_exp_cls_ali_tor"
  pdbs <- osp.uni$PDB_ID |> unique()
  pdbs |> walk(\(x) get_torsions_py(x, folder_in, folder_out))
  folder_in <- "data_csm_pdb_cls_ali"
  folder_out <- "data_csm_cls_ali_tor"
  pdbs |> walk(\(x) get_torsions_py(x, folder_in, folder_out))
}

# COMPUTES AND COMBINES CHI ANGLE DIFFS BETWEEN EXP AND CSM PDBS  
# AUXILIARY FUNCTIONS:  
# compare_chis_py: computes chi angle differences between two PDB torsion files and calculates angle diffs
if(0){
  folder_exp <- "data_exp_cls_ali_tor"
  folder_csm <- "data_csm_cls_ali_tor"
  chi_files <- paste0("chi",c(1,1,2,2,3:5)) |> paste0("_diffs.dat")
  chi_files[2] <- paste0("a",chi_files[2])
  chi_files[4] <- paste0("a",chi_files[4])
  tor.uni <- pdbs |> map_dfr(\(x) compare_chis_py(x, folder_exp, folder_csm, chi_files))
  tor.uni <- tor.uni |> rename(Resnum.exp = Resnum)
  file_delete(chi_files[file_exists(chi_files)])
  osp.uni <- osp.uni |> left_join(tor.uni, by = c("PDB_ID", "Resnum.exp", "Resname"))
  osp.uni <- osp.uni |> mutate(across(starts_with("OSP.dif"), ~ round(.x, 3)))
}

# COMPUTES ASA FOR EXP AND CSM STRUCTURES AND MERGES RESULTS  
# CALCULATES ASA DIFFERENCE AND CLASSIFIES RESIDUES AS SURFACE OR CORE  
# AUXILIARY FUNCTIONS:  
# frac_asa_py: computes relative SASA for each residue using Shrake-Rupley
if(0){
  pdbs <- osp.uni$PDB_ID |> unique()
  folder_in <- "data_exp_pdb_cls_ali"
  folder_out <- "data_exp_cls_ali_asa"
  exp.asa <- pdbs |> map(\(x) frac_asa_py(x, folder_in, folder_out))
  exp.asa <- exp.asa |> map(~ mutate(.x, Resid = row_number()) |> relocate(Resid, .after = 1)) |> 
    bind_rows()
  folder_in <- "data_csm_pdb_cls_ali"
  folder_out <- "data_csm_cls_ali_asa"
  csm.asa <- pdbs |> map(\(x) frac_asa_py(x, folder_in, folder_out))
  csm.asa <- csm.asa |> map(~ mutate(.x, Resid = row_number()) |> relocate(Resid, .after = 1)) |> 
    bind_rows()
  asa.uni <- exp.asa |> left_join(csm.asa, by=c("PDB_ID", "Resid", "Resname"), 
                                  suffix = c(".exp", ".csm"))
  osp.uni <- osp.uni |> left_join(asa.uni, 
                                  by=c("PDB_ID", "Resid", "Resnum.exp", "Resname", "Resnum.csm"))
  osp.uni <- osp.uni |> mutate(ASA.diff = round((ASA.csm - ASA.exp), 3))
  osp.uni <- osp.uni |> mutate(Position.exp = if_else(ASA.exp >= 0.25, "SURFACE", "CORE"))
  osp.uni <- osp.uni |> mutate(Position.csm = if_else(ASA.csm >= 0.25, "SURFACE", "CORE"))
}

# CALCULATES OSP MEAN, SD, AND MAD FOR EXP AND CSM STRUCTURES PER PDB_ID  
# MERGES STATISTICS WITH DB.EFF AND FILTERS PDBS WITH CUTOFF SEQUENCE DIFFERENCE
if(0){
  k <- 3
  res.stat <- tibble(PDB_ID = db.eff$PDB_ID,
                     osp.exp.mean = exp.osp |> map(\(x) mean(x$OSP)) |> unlist() |> round(k),
                     osp.csm.mean = csm.osp |> map(\(x) mean(x$OSP)) |> unlist() |> round(k),
                     osp.exp.sd = exp.osp |> map(\(x) sd(x$OSP)) |> unlist() |> round(k),
                     osp.csm.sd = csm.osp |> map(\(x) sd(x$OSP)) |> unlist() |> round(k),
                     osp.exp.mad = exp.osp |> map(\(x) mad(x$OSP)) |> unlist() |> round(k),
                     osp.csm.mad = csm.osp |> map(\(x) mad(x$OSP)) |> unlist() |> round(k))
  db.stat <- db.eff |> select(PDB_ID, !contains("path")) |> select(!contains(".id"))
  db.stat <- db.stat |> left_join(res.stat, by = "PDB_ID") 
  cutoff <- 0.05
  db.short <- db.stat |> filter(n_SEQDIFF <= cutoff*n_SEQSTR)
  osp.short <- osp.uni |> filter(PDB_ID %in% db.short$PDB_ID) 
  osp.short <- osp.short |> mutate(Position = paste(Position.exp,"|",Position.csm))
}

# FILTERS SELECTED PDBS AND RETRIEVES SSE VIA P-SEA FOR EXP AND CSM STRUCTURES  
# JOINS SSE DATA AND TRANSLATES CODES TO ALPHA/BETA/COIL LABELS
if(0){
  #x <- pdb.bio3d[db$PDB_ID %in% db.work$PDB_ID]|> map2_dfr(db.work$PDB.path, \(x,y) get_pdb_sse(x,y))
  aux <- db.eff |> filter(PDB_ID %in% db.short$PDB_ID) |> select(PDB.cls.ali.path, CSM.cls.ali.path)
  aux.exp <- db.short$PDB_ID |>  map2_dfr(aux$PDB.cls.ali.path, \(x,y) get_sse_psea(x,y)) |> 
    rename(SSE.exp = SSE)
  aux.csm <- db.short$PDB_ID |>  map2_dfr(aux$CSM.cls.ali.path, \(x,y) get_sse_psea(x,y)) |> 
    rename(SSE.csm = SSE)
  aux <- aux.exp |> left_join(aux.csm, by = c("PDB_ID", "Resid"))
  osp.short <- osp.short |> left_join(aux, by =  c("PDB_ID", "Resid"))
  #osp.short <- osp.short |> rename(SSE.csm = SSE_csm) |> rename(SSE.exp = SSE_exp)
  osp.short <- osp.short |> mutate(SSE.exp = case_when(
    SSE.exp == "a" ~ "ALPHA",
    SSE.exp == "b" ~ "BETA",
    SSE.exp == "c" ~ "COIL"
  ))
  osp.short <- osp.short |> mutate(SSE.csm = case_when(
    SSE.csm == "a" ~ "ALPHA",
    SSE.csm == "b" ~ "BETA",
    SSE.csm == "c" ~ "COIL"
  ))
  osp.short <- osp.short |> mutate(SSE = paste(SSE.exp,"|",SSE.csm))
}

# CALCULATES SQUARED DEVIATIONS AND INFLUENCE SCORES FOR EACH PDB_ID GROUP  
# FLAGS RESIDUES AS INFLUENTIAL USING VARIOUS CRITERIA (SQDIFF, INFL.DIFF, ...)  
# AUXILIARY FUNCTIONS:  
# add_stat_square_diff: calculates squared deviations of OSP.EXP and OSP.CSM from  
# values in V  
# add_stat_infl_sd: computes influence scores for OSP.EXP and OSP.CSM based on mean  
# and SD values  
# set_pack_type: defines label according to differences in metrics V in CSM and EXP  
# set_influential_res: identifies influential residues by interactively comparing  
# diff values using FUN
if(0){
  
  k = 3
  # aux <- db.short |> select(contains("mean")) |> t() |> as_tibble()
  db.short.stat <- db.short |> select(contains("osp")) |> select(1:4) |> t() |> as_tibble()
  osp.short <- osp.short |> group_by(PDB_ID) |> group_split() |>
               map2(db.short.stat |> slice(-(3:4)),
               \(x,y) add_stat_square_diff(x, y)) |> bind_rows()
  osp.short <- osp.short |> group_by(PDB_ID) |> group_split() |>
                map2(db.short.stat, \(x,y) add_stat_infl_sd(x, y)) |> bind_rows()
  osp.short$Influential = "NON-INFLUENTIAL"
  osp.short <- osp.short |> relocate(Influential, .after = Position.csm)
  osp.short$Influential_sc = "NON-INFLUENTIAL"
  osp.short <- osp.short |> relocate(Influential_sc, .after = Influential)
  osp.short$Influential_if = "NON-INFLUENTIAL"
  osp.short <- osp.short |> relocate(Influential_if, .after = Influential_sc)
  osp.short$Influential_th = "NON-INFLUENTIAL"
  osp.short <- osp.short |> relocate(Influential_th, .after = Influential_if)
  osp.short$Influential_vc = "NON-INFLUENTIAL"
  osp.short <- osp.short |> relocate(Influential_vc, .after = Influential_th)
  # db.short.sds <- db.short |> select(contains("sd")) |> t() |> as_tibble()
  # source("FIBOS-case-study-expanded-more-fun.R")
  osp.short.list <- osp.short |> group_by(PDB_ID) |> group_split()
  osp.short <- osp.short.list |> map2(db.short.stat |> slice(-(1:2)),
                                      \(x,y) set_pack_type(x, y)) |> bind_rows()
  
  osp.short.list <- osp.short |> group_by(PDB_ID) |> group_split()
  osp.short <- osp.short.list |> map2(db.short.stat |> slice(-(1:2)),
                                      \(x,y) set_influential_res(x, y, 
                                                                 cols = c("OSP.sqdiff", "Influential_sc"))) |>                 bind_rows() 
  osp.short.list <- osp.short |> group_by(PDB_ID) |> group_split()
  osp.short <- osp.short.list |> map2(db.short.stat |> slice(-(1:2)), 
                                      \(x,y) set_influential_res(x, y, 
                                                                 cols = c("OSP.infl.diff", "Influential_if"))) |>                 bind_rows() 
  infl.if <- (osp.short |> filter(Influential_if == "INFLUENTIAL"))$OSP.infl.diff
  # infl.if.med <- osp.short |> filter(Influential_if == "INFLUENTIAL") |> 
  #                select(OSP.infl.diff) |> unlist(use.names = FALSE) |> median()
  infl.if.med <- infl.if |> median()
  osp.short <- osp.short |> mutate(OSP.infl.vec = round(sqrt(OSP.infl.diff^2 + OSP.infl.sum^2), k),
                                   OSP.infl.quad = interaction(
                                     if_else(OSP.infl.diff >= 0, 1, -1),
                                     if_else(OSP.infl.sum >= 0, 1, -1)
                                   )) |> 
               relocate(OSP.infl.vec, .after = OSP.infl.sum) |> 
               relocate(OSP.infl.quad, .after = OSP.infl.vec)
  #infl.vc.cut <- median(osp.short$OSP.infl.vec)
  # infl.vc.cut <- 0.2 # pos modal
  infl.vc.cut <- osp.short |> filter(Influential_if == "INFLUENTIAL") |>
    select(OSP.infl.vec) |> unlist(use.names = FALSE) |> median()
  # infl.vc.cut <- osp.short$OSP.infl.vec |> quantile(probs = 0.9)
  osp.short <- osp.short |> 
    mutate(Influential_vc = if_else(abs(OSP.infl.vec) > infl.vc.cut, 
                                    "INFLUENTIAL", 
                                     Influential))
  osp.short <- osp.short |> 
    mutate(Influential = if_else(Influential_if == "INFLUENTIAL" |
                                   Influential_th == "INFLUENTIAL",
                                 "INFLUENTIAL",
                                 Influential))
}

# FILTERS INFLUENTIAL RESIDUES FROM THE ORIGINAL DISTRIBUTION
if(0){
  k = 3
  aux1 <- osp.short |> filter(Influential_if != "INFLUENTIAL", Packtype == "CSM > EXP") |> group_by(PDB_ID) |> 
    summarise(
      osp.csm.mean = round(mean(OSP.csm), k),
      osp.csm.sd = round(sd(OSP.csm), k)) |> 
    ungroup()
  x1 = db.short |> select(!contains("csm.mean") & !contains("csm.sd")) |> left_join(aux1, by = "PDB_ID") |> 
    select(PDB_ID, osp.exp.mean, osp.csm.mean, osp.exp.sd, osp.csm.sd) |> drop_na()
  aux2 <- osp.short |> filter(Influential_if != "INFLUENTIAL", Packtype == "CSM < EXP") |> group_by(PDB_ID) |> 
    summarise(
      osp.exp.mean = round(mean(OSP.exp), k),
      osp.exp.sd = round(sd(OSP.exp), k)) |> 
    ungroup()
  x2 = db.short |> select(!contains("exp.mean") & !contains("exp.sd")) |> left_join(aux2, by = "PDB_ID") |> 
    select(PDB_ID, osp.exp.mean, osp.csm.mean, osp.exp.sd, osp.csm.sd) |> drop_na()
  aux3 = osp.short |> filter(Packtype == "CSM == EXP") |> select(PDB_ID) |> 
    unlist(use.names = FALSE) |> unique()
  x3 <- db.short |> filter(PDB_ID %in% aux3) |> 
    select(PDB_ID, osp.exp.mean, osp.csm.mean, osp.exp.sd, osp.csm.sd)
  db.short.clean <- bind_rows(x1, x2, x3) |> slice(match(db.short$PDB_ID, PDB_ID))
  osp.short.infl.if <- osp.short |> filter(Influential_if == "INFLUENTIAL")
  osp.short.infl <- osp.short |> filter(Influential == "INFLUENTIAL")
}

# SAVES CSV OF DATASETS
if(0){
  
  # extended dataset
  if(0) {
    folder <- "data"
    file <- "FIBOS_extended-dataset-titles-v1.csv"
    db.eff |> select(PDB_ID, TITLE) |> write_csv(path(folder, file))
    file <- "FIBOS_extended-dataset-info-v1.csv"
    db.eff |> select(PDB_ID, !TITLE) |> select(!contains("path")) |> write_csv(path(folder, file))
    file <- "FIBOS_extended-dataset-all-v1.csv"
    db.eff |> select(!contains("path")) |> write_csv(path(folder, file))
  }
  
  # strict dataset
  if(0) {
    folder <- "data"
    file <- "FIBOS_strict-dataset-titles-v1.csv"
    db.short |> select(PDB_ID, TITLE) |> write_csv(path(folder, file))
    file <- "FIBOS_strict-dataset-info-v1.csv"
    db.short |> select(PDB_ID, !TITLE) |> select(!contains("path")) |> write_csv(path(folder, file))
    file <- "FIBOS_strict-dataset-all-v1.csv"
    db.short |> select(!contains("path")) |> select(!ends_with("mad")) |> write_csv(path(folder, file))
    file <- "FIBOS_strict-dataset-clean-v1.csv"
    db.short.clean |>  write_csv(path(folder, file))
  }
  
  # osp by residue dataset
  if(0) {
    folder <- "data"
    file <- "FIBOS_by_residue-dataset-all-v1.csv"
    osp.short |> select(!contains("path")) |> rename(inf = Influential_if) |> 
      select(!contains("Influential")) |> rename(Influential_if = inf) |> 
      write_csv(path(folder, file))
  }
}

