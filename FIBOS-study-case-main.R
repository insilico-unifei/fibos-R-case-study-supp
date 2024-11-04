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
# CONTACTS: hersinsoares@gmail.com, carlos.silveira@unifei.edu.br
#
# FIBOS install on R:
# See tutorial here: https://github.com/insilico-unifei/fibos-R 

# LIBRARIES
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


# SOURCE OF AUXILIARY FUNCTIONS
source("FIBOS-study-case-fun.R")

# CHANGE CONDITIONALS TO 1 TO EXECUTE EACH STEP

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
  pdb.meta <- pdb.ids |> map(\(x) get_pdb_metadata(x, source, par))
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

# CREATE THE TABLES FOR COMPARATIVE CALCULATIONS OF THE OS AND OSP
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
  
  # filter main table to produce work table with difference between AF and Experimental less than 5 aa in size
  pdb.csm.tab.work <- pdb.csm.tab |> filter(abs(aaDIFF) < 5)
  
  # create case study output table
  if(0) write_csv(pdb.csm.tab.work |> mutate(across(where(is.list), as.character)), 
                  file = "data/case_study_v1.csv")
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
  pdb.csm.tab.work<- pdb.csm.tab.work |> 
                     mutate(SRF.path = path("fibos_files", paste0("prot_", PDB.ids), ext = "srf")) |> 
                     relocate(SRF.path, .after = CSM.path.new)
}

# CALCULATE OCCLUDE SURFACE AT ATOM AND RESIDUE LEVEL FOR EXPERIMENTAL PDB
# WARNING: This may take a while depending on your hardware configuration and the number of cores available
if(0){

  # calculate OS at atom level in parallel with max cores available according to PDB_ids size
  tic()
  my_default_mccores = getOption("mc.cores")
  my_ideal_mccores = min(parallel::detectCores(), length(pdb.csm.tab.work$PDB.path))
  if (my_ideal_mccores > 0) options(mc.cores = my_ideal_mccores)
  if (my_ideal_mccores > 1) plan(multisession, workers = my_ideal_mccores)
  pdb.exp.fibos <- pdb.csm.tab.work$PDB.path |> future_map(\(x) occluded_surface (x, method = "FIBOS"), 
                                                .options = furrr_options(seed = 123))
  if (my_ideal_mccores > 0) options(mc.cores = my_default_mccores)
  toc()
  names(pdb.exp.fibos) <- paste0(pdb.csm.tab.work$PDB.ids,"_exp")
  
  #calculate the OSP metric at the residue level
  pdb.exp.osp <- pdb.csm.tab.work$SRF.path |> map(\(x) osp(x))
  names(pdb.exp.osp) <- paste0(pdb.csm.tab.work$PDB.ids,"_exp")
  
  # save the calculation files in the "fibos_files_pdb_exp" folder
  file.rename("fibos_files","fibos_files_exp_pdb")
  

}

# CALCULATE OCCLUDE SURFACE AT ATOM AND RESIDUE LEVEL FOR AF MODELS
if(0){
  
  # calculate OS at atom level in parallel with max cores available according to PDB_ids size
  tictoc::tic()
  my_default_mccores = getOption("mc.cores")
  my_ideal_mccores = min(parallel::detectCores(), length(pdb.csm.tab.work$CSM.path.new))
  if (my_ideal_mccores > 0) options(mc.cores = my_ideal_mccores)
  if (my_ideal_mccores > 1) plan(multisession, workers = my_ideal_mccores)
  pdb.csm.fibos <- pdb.csm.tab.work$CSM.path.new |> future_map(\(x) occluded_surface(x, method = "FIBOS"), 
                                                               .options = furrr_options(seed = 123))
  if (my_ideal_mccores > 0) options(mc.cores = my_default_mccores)
  tictoc::toc()
  names(pdb.csm.fibos) <- paste0(pdb.csm.tab.work$PDB.ids,"_af")
  
  #calculate the OSP metric at the residue level
  pdb.csm.osp <- pdb.csm.tab.work$SRF.path |> map(\(x) osp(x))
  names(pdb.csm.osp) <- paste0(pdb.csm.tab.work$PDB.ids,"_af")
  
  # save the calculation files in the "fibos_files_pdb_exp" folder
  file.rename("fibos_files","fibos_files_csm_pdb")
  
}

# MAKE OSP UNIFIED TABLES
if(0){
  
  folder <- "data"
  
  file <- "case-study-osp-pdb-unique-tab.csv"
  pdb.exp.osp.uni <- pdb.exp.osp |> bind_rows(.id = "ID")
  if(0) write_csv(pdb.exp.osp.uni, path(folder, file))
  
  file <- "case-study-osp-af-unique-tab.csv"
  pdb.csm.osp.uni <- pdb.csm.osp |> bind_rows(.id = "ID")
  if(0) write_csv(pdb.csm.osp.uni, path(folder, file))
  
}

# MAKE OSP STATISTICAL TABLE
if(0){
  
  res.stat <- tibble(ID = pdb.csm.tab.work$PDB.ids,
                     osp.exp.mean = pdb.exp.osp |> map(\(x) mean(x$OSP)) |> unlist(),
                     osp.af.mean = pdb.csm.osp |> map(\(x) mean(x$OSP)) |> unlist(),
                     osp.exp.sd = pdb.exp.osp |> map(\(x) sd(x$OSP)) |> unlist(),
                     osp.af.sd = pdb.csm.osp |> map(\(x) sd(x$OSP)) |> unlist()
  )
  
  # SET TO 1 TO SAVE IT
  if(0) res.stat |> mutate(across(where(is.numeric), ~round(.x, 3))) |> 
        write_csv("data/osp-mean-sd-exp-af.csv")
}

# MAKE PLOT OF MAIN FIGURE 02 OF THE ARTCILE
if(0){
  par.lines <- tibble(
    mean = round(c(mean(res.stat$osp.exp.mean),mean(res.stat$osp.af.mean)),3),
    sd =  round(c(mean(res.stat$osp.exp.sd),mean(res.stat$osp.af.sd)),3),
    type = c("solid", "81"),
    color = c("#F8766D","#00BFC4")
  )
  
  
  p1 <- res.stat |> rename(`Experimental PDB` = osp.exp.mean, `Predicted AF` = osp.af.mean) |>
    select(!ID) |> 
    pivot_longer(!ends_with("sd"), names_to = "TYPE:", 
                 values_to = "mean(OSP)") |> 
    ggplot(aes(x = `mean(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity",  linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) + 
    scale_x_continuous(breaks = seq(0.38, 0.48, by = 0.04)) +
    theme_classic() + 
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust = -1.5),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
  #print(p1)
  
  p2 <- res.stat |> rename(`Experimental PDB` = osp.exp.sd, `Predicted AF` = osp.af.sd) |>
    select(!ID) |>
    pivot_longer(!ends_with("mean"), names_to = "TYPE:",
                 values_to = "sd(OSP)") |>
    ggplot(aes(x = `sd(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust = -1.5),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
  #print(p2)
  #
  title <- "Protein Packing Density Distributions for Predicted and Experimental Models\n"
  p12 <- ggarrange(p1, p2, ncol=2, labels=c("A)","B)"), common.legend = T,
                   legend = "bottom")
  p12 <- annotate_figure(p12, top = text_grob(title,
                                              color = "black", face = "plain", size = 16))
  print(p12)

}

# SOME USEFUL STATISTICS
if(0){
  
  wilcox.test(res.stat$osp.exp.mean, res.stat$osp.af.mean, paired = T)
  ks.test(res.stat$osp.exp.mean, res.stat$osp.af.mean)
  cohen.d(res.stat$osp.exp.mean, res.stat$osp.af.mean)
  cohen.d(res.stat$osp.exp.mean, res.stat$osp.af.mean, hedges.correction=T)
  
  wilcox.test(res.stat$osp.exp.sd, res.stat$osp.af.sd, paired = T)
  ks.test(res.stat$osp.exp.sd, res.stat$osp.af.sd)
  cohen.d(res.stat$osp.exp.sd, res.stat$osp.af.sd)
  cohen.d(res.stat$osp.exp.sd, res.stat$osp.af.sd, hedges.correction=T)
  
}

