###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# R SCRIPT OF AUXILIARY FUNCTIONS USED BY FIBOS-study-case.main.R
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

# GET PDB METADATA FROM SOURCE URL AND PAR AND RETURN IT IN JSON FORMAT
get_pdb_metadata <- function(pdb,
                             source,
                             par = "",
                             verbose = T){
  
  # form url
  URL <- paste0(source, pdb, par)
  
  # test url
  con <- suppressWarnings(try(url(URL,"r"), silent = T))
  
  # if error, then warning 
  if (inherits(con, "try-error")) {
    warning(paste("FAILED to get metadata for", URL))
    return(NA)
  }
  
  if (verbose) paste("Getting metadata for", pdb) |> print()
  
  # read the content
  res <- con |> readLines(warn = F)
  
  # close connection
  close(con)
  
  # put the content in json format
  data <- res |> paste(collapse="") |> jsonlite::fromJSON()
  
  return(data)
  
}

# GET UNIPROT ID FROM PDB METADATA
get_pdb_uniprot_id <- function(data, verbose = T) {
  
  # get pdb id from metadata
  pdb.id <- data$rcsb_polymer_entity_container_identifiers$entry_id
  
  if (verbose) paste("Getting UniProt id for", pdb.id) |> print()
  
  # get vector of identifiers
  data.seq.ids <- data$rcsb_polymer_entity_container_identifiers$reference_sequence_identifiers |>
    list()
  
  # get uniprot id
  uniprot.id <- data.seq.ids |> map_vec(\(x) ifelse(x["database_name"] == "UniProt", 
                                                    x["database_accession"], ""))
  # If more than one uniprot id, warning
  if (length(uniprot.id) > 1) warning(paste0("WARNING: there is more than one UniProt id for "
                                             ,pdb.id)) |> print()

  return(unlist(uniprot.id))
}

# GET AA SEQUENCE FROM PDB METADATA
get_pdb_aa_seq <- function(data, verbose = T){
  
  # get pdb id from metadata
  pdb.id <- data$rcsb_polymer_entity_container_identifiers$entry_id
  
  if (verbose) paste("Getting aa sequence id for", pdb.id) |> print()
  
  # get aa sequence from metadata
  data.seq <- data$entity_poly$pdbx_seq_one_letter_code_can
  
  return(data.seq)
  
}

# DOWNLOAD AF IN PDB FORMAT AND RETURN AF PATHs
get_csm <- function(csm.id, 
                    source = "https://alphafold.ebi.ac.uk/files/",
                    path = "",
                    suffix = ".pdb",
                    verbose = T,
                    overwrite = F){

  # form url and file name
  csm.name <- paste0(csm.id, suffix)
  url <- paste0(source, csm.name)
  file <- path(path, csm.name)
  
  # if file doesn't exist or overwrite
  if (!file.exists(file) | overwrite){
    
    # get url
    res <- GET(url)
    
    # if get ok ... 
    if (res$status_code == 200) {
      
      if (verbose) paste("Getting CSM", url) |> print()
      
      # write content as text in file name
      writeBin(content(res, "raw"), file)
      
      # return path
      return(file)
    # if get not ok, warning
    } else {
      paste("FAILED to get CSM", url, "Status code:", res$status_code) |> 
        print()
      
      # return NA
      return(NA)
    }
  # if file exist, warning skipping download
  } else {
    warning(paste(file, " exists. Skipping download"))
    
    # return path
    return(file)
  }
}


