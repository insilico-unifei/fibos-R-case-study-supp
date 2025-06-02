###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# R SCRIPT OF AUXILIARY FUNCTIONS USED BY FIBOS-case-study-expanded-main.R
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
  
  #browser()
  # get pdb id from metadata
  pdb.id <- data$rcsb_polymer_entity_container_identifiers$entry_id
  
  if (verbose) paste("Getting UniProt id for", pdb.id) |> print()
  
  # get vector of identifiers
  data.seq.ids <- data$rcsb_polymer_entity_container_identifiers$reference_sequence_identifiers |>
    list()
  #browser()
  # get uniprot id
  
  if(is.null(data.seq.ids[[1]])){
    print(paste("WARNING: there is some problem with UniProt id for", pdb.id, "Returning NA"))
    return(NA)
  }
  
  if(dim(data.seq.ids[[1]])[1]>1) {
    print(paste("WARNING: there is more than one UniProt id for", pdb.id, "Returning NA"))
    return(NA)
  }
  uniprot.id <- data.seq.ids |> map_vec(\(x) ifelse(x["database_name"] == "UniProt", 
                                                    x["database_accession"], ""))
  # If more than one uniprot id, warning
  if (length(uniprot.id) > 1) {
    print(paste("WARNING: there is more than one UniProt id for", pdb.id, "Returning NA"))
    return(NA)
  }

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

# GET SOME EXPERIMENTAL PDB PARAMETERS LIKE TITLE, RESOLUTION, R_FACTORS
get_pdb_info = function(pdb, k = 2){
  
  aux = get_info(pdb)
  #browser()
  res = tibble(PDB_ID = pdb,
               TITLE = aux$struct$title,
               #RESIDUES = aux$rcsb_entry_info$deposited_polymer_monomer_count,
               RESOLUTION = aux$rcsb_entry_info$resolution_combined,
               R_FACTOR = aux$refine[1,"ls_rfactor_obs"],
               R_FACTOR_W = aux$refine[1,"ls_rfactor_rwork"])
  
  return(res)
  
}


# RETURN A TIBBLE TABLE INDICATING POSSIBLE PDB WITH MISSING RESIDUES
verify_missing_residues = function(pdb, pdbid){
  
  #pdbid <- read.pdb(pdb.path, verbose = FALSE, rm.alt = FALSE)

  
  str_tab <- pdb$atom |>  filter(type == "ATOM") |> select(resid, resno) |> 
              group_by(resid, resno) |> distinct()
  
  
  seqres <- pdb$seqres
  
  res_tab <- tibble(PDB_ID = pdbid,
                    n_SEQSTR = dim(str_tab)[1],
                    n_SEQRES = length(seqres)) #,
                    #n_SEQDIFF = n_SEQRES - n_SEQSTR)

  miss_ind <- c(str_tab$resno, last(str_tab$resno)+1) |> diff() |> {\(x) abs(x) > 1}()
  
  miss_tab <- str_tab[miss_ind,] |> mutate(residno = paste0(aa321(resid),resno))
  
  res_tab$MISS_POINT = paste(miss_tab$residno, collapse = ",")
  
  return(res_tab)
}

# GET SSE INFORMATION USING DSSP
get_pdb_sse <-  function(pdb, pdb.path, exefile = "/usr/bin/dssp", 
                       a = c("G","H","I"), b = c("B","E"), k = 0){
  
  #browser()
  #sse = dssp(pdb)
  sse = dssp_adapted(pdb, file.path(getwd(), pdb.path), exefile)
  types = table(sse$sse)
  all = sum(pdb$calpha)
  #browser()
  f = names(types) %in% a
  alpha = round(sum(types[f])/all*100,k)
  f = names(types) %in% b
  beta = round(sum(types[f])/all*100,k)
  res = tibble(PDB_ID = pdb.path |> path_file() |> path_ext_remove(),
               HELIX = alpha,
               STRAND = beta)
  return(res)
}

# GET SSE INFORMATION USING P-SEA
get_sse_psea <- function(pdb, pdb_path,
                         python_cmd = "python3",
                         script     = "sse_psea.py") {

  sse_json <- system2(python_cmd,
                      args   = c(script, pdb_path),
                      stdout = TRUE,
                      stderr = TRUE)
  sse_str <- paste(sse_json, collapse = "")
  sse_vec <- fromJSON(sse_str) |> as.character()
  res <- tibble(PDB_ID = pdb,
                Resid = 1:length(sse_vec),
                SSE = sse_vec)
  return(res)
}

# ADAPTED DSSP FUNCTION FROM BIO3D
# (Need dssp installed in SO)
dssp_adapted = function (pdb, pdb.path, exefile = "/usr/bin/dssp", par = "--output-format=dssp", 
                    resno = TRUE, full = FALSE, verbose = FALSE, ...) 
{
  cl <- match.call()
  os1 <- Sys.info()["sysname"]
  #exefile <- .get.exepath(exefile)
  #message(exefile)
  #success <- .test.exefile(exefile)
  success <- file_exists(exefile)
  if (!success) {
    stop(paste("Launching external program 'dssp' (or 'mkdssp') failed\n", 
               "  make sure '", exefile, "' is in your search path", 
               sep = ""))
  }
  checkatoms <- TRUE
  if (checkatoms) {
    inds <- atom.select(pdb, "backbone", verbose = verbose)
    tmp <- trim.pdb(pdb, inds)
    resid <- paste(tmp$atom$resno, tmp$atom$chain, sep = "-")
    musthave <- c("C", "CA", "N", "O")
    incomplete <- sapply(unique(resid), function(x) {
      inds <- which(resid == x)
      elety <- sort(tmp$atom$elety[inds])
      if (!all(musthave %in% elety)) 
        return(TRUE)
      else return(FALSE)
    })
    if (all(incomplete)) 
      stop("No residues found with a complete set of backbone atoms")
    if (any(incomplete)) 
      warning(paste("Residues with missing backbone atoms detected:", 
                    paste(unique(resid)[incomplete], collapse = ", "), 
                    collapse = " "))
    inds <- atom.select(pdb, "protein", verbose = verbose, 
                        inverse = TRUE)
    if ((length(inds$atom) > 0) && (verbose))
      warning(paste("Non-protein residues detected in input PDB:", 
                    paste(unique(pdb$atom$resid[inds$atom]), collapse = ", ")))
  }
  #infile <- tempfile()
  infile <- pdb.path
  outfile <- tempfile()
  #browser()
  #write.pdb(pdb, file = infile)
  cmd <- paste(exefile, par, infile, outfile)
  if (verbose) 
    cat(paste("Running command:\n ", cmd, "\n"))
  if (os1 == "Windows") {
    status <- shell(paste(shQuote(exefile), infile, outfile), 
                    ignore.stderr = !verbose, ignore.stdout = !verbose)
  }
  else {
    status <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)
  }
  if (!(status %in% c(0, 1))) {
    stop(paste("An error occurred while running command\n '", 
               cmd, "'", sep = ""))
  }
  trim <- function(s) {
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s == "")] <- NA
    s
  }
  split.line <- function(x, split = " ") {
    tmp <- unlist(strsplit(x, split = split))
    inds <- which(tmp != "")
    return(trim(tmp[inds]))
  }
  #browser()
  raw.lines <- base::readLines(outfile)
  #unlink(c(infile, outfile))
  unlink(c(outfile))
  type <- substring(raw.lines, 1, 3)
  raw.lines <- raw.lines[-(1:which(type == "  #"))]
  aa <- substring(raw.lines, 14, 14)
  if (any(aa == "!")) 
    raw.lines <- raw.lines[-which(aa == "!")]
  cha <- substring(raw.lines, 12, 12)
  sse <- substring(raw.lines, 17, 17)
  res.name <- substring(raw.lines, 14, 14)
  res.id <- as.numeric(substring(raw.lines, 1, 5))
  res.num <- as.numeric(substring(raw.lines, 6, 10))
  res.insert <- substring(raw.lines, 11, 11)
  res.ind <- 1:length(res.num)
  ins <- trim(res.insert)
  ins[ins == ""] = NA
  names(sse) <- paste(res.num, cha, ins, sep = "_")
  if (any(res.insert != " ")) {
    if ((resno)&&(verbose)) {
      warning("Insertions are found in PDB: Residue numbers may be incorrect.\n                Try again with resno=FALSE")
    }
    else {
      ii <- diff(res.num)
      ii[ii == 0] <- 1
      ii[ii < 0] <- 2
      res.num <- res.num[1] + c(0, cumsum(ii))
    }
  }
  if (full) {
    diff <- res.id - res.ind
    names(diff) <- res.id
    bp1 <- as.numeric(substring(raw.lines, 26, 29))
    bp2 <- as.numeric(substring(raw.lines, 30, 33))
    bp1[bp1 == 0] <- NA
    bp2[bp2 == 0] <- NA
    bp1[!is.na(bp1)] <- as.vector(bp1[!is.na(bp1)] - diff[as.character(bp1[!is.na(bp1)])])
    bp2[!is.na(bp2)] <- as.vector(bp2[!is.na(bp2)] - diff[as.character(bp2[!is.na(bp2)])])
    hbonds <- split.line(split.line(substring(raw.lines, 
                                              40, 83), split = ","), split = " ")
    hbonds <- matrix(as.numeric(hbonds), ncol = 8, byrow = TRUE)
    hbonds <- as.data.frame(hbonds)
    for (i in seq(1, 7, by = 2)) {
      hbonds[[i]][which(hbonds[[i]] == 0)] <- NA
      hmm <- res.id + hbonds[[i]]
      hbonds[[i]][!is.na(hmm)] <- as.vector(hmm[!is.na(hmm)] - 
                                              diff[as.character(hmm[!is.na(hmm)])])
    }
    hbonds <- cbind(bp1, bp2, hbonds)
    cnames <- c("BP1", "BP2", "NH-O.1", "E1", "O-HN.1", 
                "E2", "NH-O.2", "E3", "O-HN.2", "E4")
    colnames(hbonds) <- cnames
    if (resno) {
      tmp.map <- cbind(res.num, cha)
      row.names(tmp.map) <- res.ind
      hbonds <- cbind(hbonds, data.frame(matrix(NA, ncol = 6, 
                                                nrow = nrow(tmp.map)), stringsAsFactors = FALSE))
      colnames(hbonds) <- c(cnames, "ChainBP1", "ChainBP2", 
                            "Chain1", "Chain2", "Chain3", "Chain4")
      tmp.inds <- which(!is.na(hbonds[, "BP1"]))
      tmp.names <- as.character(hbonds[tmp.inds, "BP1"])
      hbonds[tmp.inds, "BP1"] <- as.numeric(tmp.map[tmp.names, 
                                                    "res.num"])
      hbonds[tmp.inds, "ChainBP1"] <- tmp.map[tmp.names, 
                                              "cha"]
      tmp.inds <- which(!is.na(hbonds[, "BP2"]))
      tmp.names <- as.character(hbonds[tmp.inds, "BP2"])
      hbonds[tmp.inds, "BP2"] <- as.numeric(tmp.map[tmp.names, 
                                                    "res.num"])
      hbonds[tmp.inds, "ChainBP2"] <- tmp.map[tmp.names, 
                                              "cha"]
      tmp.inds <- which(!is.na(hbonds[, "NH-O.1"]))
      tmp.names <- as.character(hbonds[tmp.inds, "NH-O.1"])
      hbonds[tmp.inds, "NH-O.1"] <- as.numeric(tmp.map[tmp.names, 
                                                       "res.num"])
      hbonds[tmp.inds, "Chain1"] <- tmp.map[tmp.names, 
                                            "cha"]
      tmp.inds <- which(!is.na(hbonds[, "O-HN.1"]))
      tmp.names <- as.character(hbonds[tmp.inds, "O-HN.1"])
      hbonds[tmp.inds, "O-HN.1"] <- as.numeric(tmp.map[tmp.names, 
                                                       "res.num"])
      hbonds[tmp.inds, "Chain2"] <- tmp.map[tmp.names, 
                                            "cha"]
      tmp.inds <- which(!is.na(hbonds[, "NH-O.2"]))
      tmp.names <- as.character(hbonds[tmp.inds, "NH-O.2"])
      hbonds[tmp.inds, "NH-O.2"] <- as.numeric(tmp.map[tmp.names, 
                                                       "res.num"])
      hbonds[tmp.inds, "Chain3"] <- tmp.map[tmp.names, 
                                            "cha"]
      tmp.inds <- which(!is.na(hbonds[, "O-HN.2"]))
      tmp.names <- as.character(hbonds[tmp.inds, "O-HN.2"])
      hbonds[tmp.inds, "O-HN.2"] <- as.numeric(tmp.map[tmp.names, 
                                                       "res.num"])
      hbonds[tmp.inds, "Chain4"] <- tmp.map[tmp.names, 
                                            "cha"]
      row.names(hbonds) <- apply(tmp.map, 1, paste, collapse = "-")
    }
  }
  else {
    hbonds <- NULL
  }
  phi <- as.numeric(substring(raw.lines, 104, 109))
  psi <- as.numeric(substring(raw.lines, 110, 115))
  acc <- as.numeric(substring(raw.lines, 35, 38))
  h.res <- bounds(res.num[which(sse == "H")], pre.sort = FALSE)
  g.res <- bounds(res.num[which(sse == "G")], pre.sort = FALSE)
  e.res <- bounds(res.num[which(sse == "E")], pre.sort = FALSE)
  t.res <- bounds(res.num[which(sse == "T")], pre.sort = FALSE)
  h.ind <- h.res
  g.ind <- g.res
  e.ind <- e.res
  t.ind <- t.res
  if (length(h.res) > 0) {
    res.ind <- which(sse == "H")
    h.ind[, "end"] <- res.ind[cumsum(h.res[, "length"])]
    h.ind[, "start"] <- h.ind[, "end"] - h.res[, "length"] + 
      1
  }
  if (length(g.res) > 0) {
    res.ind <- which(sse == "G")
    g.ind[, "end"] <- res.ind[cumsum(g.res[, "length"])]
    g.ind[, "start"] <- g.ind[, "end"] - g.res[, "length"] + 
      1
  }
  if (length(e.res) > 0) {
    res.ind <- which(sse == "E")
    e.ind[, "end"] <- res.ind[cumsum(e.res[, "length"])]
    e.ind[, "start"] <- e.ind[, "end"] - e.res[, "length"] + 
      1
  }
  if (length(t.res) > 0) {
    res.ind <- which(sse == "T")
    t.ind[, "end"] <- res.ind[cumsum(t.res[, "length"])]
    t.ind[, "start"] <- t.ind[, "end"] - t.res[, "length"] + 
      1
  }
  if (!resno) {
    h.res <- h.ind
    g.res <- g.ind
    e.res <- e.ind
    t.res <- t.ind
  }
  sheet = list(start = NULL, end = NULL, length = NULL, chain = NULL)
  helix = list(start = NULL, end = NULL, length = NULL, chain = NULL, 
               type = NULL)
  turn = sheet
  if (length(h.res) > 1) {
    helix$start = c(helix$start, h.res[, "start"])
    helix$end = c(helix$end, h.res[, "end"])
    helix$length = c(helix$length, h.res[, "length"])
    helix$chain = c(helix$chain, cha[h.ind[, "start"]])
    helix$type = c(helix$type, sse[h.ind[, "start"]])
  }
  if (length(g.res) > 1) {
    helix$start = c(helix$start, g.res[, "start"])
    helix$end = c(helix$end, g.res[, "end"])
    helix$length = c(helix$length, g.res[, "length"])
    helix$chain = c(helix$chain, cha[g.ind[, "start"]])
    helix$type = c(helix$type, sse[g.ind[, "start"]])
  }
  if (length(helix$start) > 0) 
    helix <- lapply(helix, function(x) {
      names(x) <- 1:length(helix$start)
      return(x)
    })
  if (length(e.res) > 1) {
    sheet$start = c(sheet$start, e.res[, "start"])
    sheet$end = c(sheet$end, e.res[, "end"])
    sheet$length = c(sheet$length, e.res[, "length"])
    sheet$chain = c(sheet$chain, cha[e.ind[, "start"]])
  }
  if (length(sheet$start) > 0) 
    sheet <- lapply(sheet, function(x) {
      names(x) <- 1:length(sheet$start)
      return(x)
    })
  if (length(t.res) > 1) {
    turn$start = c(turn$start, t.res[, "start"])
    turn$end = c(turn$end, t.res[, "end"])
    turn$length = c(turn$length, t.res[, "length"])
    turn$chain = c(turn$chain, cha[t.ind[, "start"]])
  }
  if (length(turn$start) > 0) 
    turn <- lapply(turn, function(x) {
      names(x) <- 1:length(turn$start)
      return(x)
    })
  out <- list(helix = helix, sheet = sheet, hbonds = hbonds, 
              turn = turn, phi = phi, psi = psi, acc = acc, sse = sse, 
              call = cl)
  class(out) <- "sse"
  return(out)
}

# ADD PDB-ID AS ATTR
add_attr = function(m, value, attr){
  #browser()
  attr(m, attr) <- value
  return(m)
}

# FILTER RES ID LESS THAN ZERO
filter_fibos_id_zero = function(m, colres, pdbid, debug_column ="ATOM"){
  
  #browser()
  aux <- m |> getElement(colres) |> as.integer()
  f <- aux <= 0
  res.zero <- sum(f)
  if (res.zero){
    print(paste("WARNING: in", pdbid, "found residue id zero or less in", res.zero,
                "contacts, column:", colres, ". Filtering them..."))
    print(m |> filter(f) |> select(all_of(debug_column)))
    m <- m |> filter(!f)
  }
  return(m)
}

# RECOVER RESIDUE ID INFORMATIONS FROM ORIGINAL PDB 
adjust_osp_table = function(m, pdb,
                            atom.type = "ATOM",
                            resname = c("Resnum"),
                            chainname = c("Chain"),
                            sufix = "o",
                            chain_ids = NULL,
                            resnum.original = FALSE,
                            verbose = FALSE){
  
  #browser()
  
  pdbid <- attr(m, "PDB")
  if(verbose) print(paste("adjusting OSP table for", pdbid))
  pdb <- pdb$atom |> dplyr::filter(type %in% atom.type) |> dplyr::select(chain, resno) |> dplyr::distinct(chain, resno)
  pdb.len <- dim(pdb)[1]
  resnameo <- paste0(resname,sufix)
  m <- m |> mutate(across(all_of(c(resname)), as.integer)) 
  m <- filter_fibos_id_zero(m, resname, pdbid, c(resname, chainname))
  fib.id1 <- m |> getElement(resname) |> as.integer()
  fib.len <- fib.id1 |> unique() |> length()
  if (pdb.len != fib.len){
    print(paste("WARNING: something wrong in", pdbid, "pdb.len != fib.len", pdb.len, fib.len))
    #return(m)
  }
  #browser()
  m <- m |> mutate({{resnameo}} := pdb$resno[fib.id1]) |> 
    mutate(!!{{chainname}} := pdb$chain[fib.id1]) 
  #browser()
  col_seq <- c(resname, resnameo)
  m <- m %>% relocate(all_of(col_seq))
  col_seq <- chainname
  #browser()
  m <- m |> relocate(all_of(col_seq), .after = 3)
  if(resnum.original) {
    #browser()
    m <- m |> dplyr::select(!all_of(resname))
    colnames(m)[1] = resname
  }
  attr(m, "PDB") <- pdbid
  return(m) 
}

# GET SIDE CHAIN CHI ANGLES VIA PYTHON SCRIPT
get_torsions_py <- function(pdb, folder_in, folder_out, 
                                  script_path = "get_torsions.py", chis = "-c 1 2 3 4 5",
                                  venv = ".venv",
                                  verbose = TRUE){
  
  pdb <- toupper(pdb)
  pdb_path <- path(folder_in, pdb, ext = "pdb")
  if (!dir.exists(folder_out)) dir_create(folder_out)
  tor_path <- path(folder_out, pdb, ext = "dat")
  args <- paste(pdb_path, "-o", tor_path, chis)
  #browser()
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    python_venv <- path(venv, "Scripts", "python3.exe")
  } else {
    python_venv <- path(venv, "bin", "python3")
  }
  if (!file.exists(python_venv)) {
    stop("Python not found in path: ", python_venv)
  }
  if (verbose) print(paste("Getting torsions for: ", pdb))
  system2(python_venv, args = c(script_path, args))
}

# COMPUTES CHI ANGLE DIFFERENCES BETWEEN TWO PDB TORSION FILES AND CALCULATES ANGLE DIFFS
# WRITES RESULTS FOR CHI1–CHI5 AND ALT CHI ANGLES TO SEPARATE OUTPUT DAT FILES
# REUNITE THEM IN A TIBBLE
compare_chis_py <- function(pdb, folder_exp, folder_csms, chi_files,
                            script_path = "compare_chis2.py",
                            venv = ".venv",
                            verbose = TRUE){
  
  pdb <- toupper(pdb)
  pdb_path <- path(folder_exp, pdb, ext = "dat")
  csm_path <- path(folder_csm, pdb, ext = "dat")
  args <- paste(pdb_path, csm_path)
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    python_venv <- path(venv, "Scripts", "python3.exe")
  } else {
    python_venv <- path(venv, "bin", "python3")
  }
  if (!file.exists(python_venv)) {
    stop("Python not found in path: ", python_venv)
  }
  if (verbose) print(paste("Comparing torsions for: ", pdb))
  system2(python_venv, args = c(script_path, args))
  tor_tab <- chi_files[1] |> read_csv(col_names = c("PDB_ID", "Resname", "Resnum","chi1"), 
                                      show_col_types = FALSE)
  #browser()
  aux <-  chi_files[-1] |> 
          map(\(x) read_csv(x, col_names = c("resid","chi"), show_col_types = FALSE)) |> 
          bind_cols() |> select(starts_with("chi"))
  tor_tab <- tor_tab |> bind_cols(aux)
  colnames(tor_tab)[-(1:3)] <- gsub("_diffs.dat", "", chi_files)
  #browser()
  tor_tab <- tor_tab |> mutate(PDB_ID = pdb) |> relocate(PDB_ID) |> 
             relocate(Resname, .after = Resnum)
  return(tor_tab)
}

# CREATES PYTHON VENV IN SPECIFIED DIR AND INSTALLS REQUIREMENTS.
# HANDLES OS DIFFERENCES TO LOCATE PIP EXECUTABLE.
# VALIDATES PRESENCE OF PIP AND REQUIREMENTS FILE.
prepair_python_venv <- function(folder.venv, 
                                python_exe = "python3",
                                pip_exe = "pip",
                                requirements = "requirements.txt"){
  
  #browser()
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    pip_exe <- path(pip_exe, ext = "exe")
    pip_venv <- path(folder.venv, "Scripts", pip_exe)
    python_exe <- path(python_exe, ext = "exe")
  } else {
    pip_venv <- path(folder.venv, "bin", pip_exe)
  }

  args <- paste("-m venv", folder.venv)
  system2(python_exe, args = args)
  
  if (!file.exists(pip_venv)) {
    stop("pip not found in path: ", pip_venv)
  }
  
  if (!file.exists(requirements)) {
    stop("requirements file not found: ", requirements)
  }
  args <- paste("install -r", requirements)
  system2(pip_venv, args = args)

}

# INVOKES PYTHON SCRIPT WITHIN VENV TO DOWNLOAD RAW PDB IDs TO 
# SPECIFIC PATH IN "file".
get_raw_pdbs_py <- function(file,
                            script_path = "get_raw_pdbs.py",
                            python_exe = "python3",
                            folder.venv = ".venv",
                            verbose = TRUE){
  
  #browser()
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    python_exe <- path(python_exe, ext = "exe")
    python_venv <- path(folder.venv, "Scripts", python_exe)
  } else {
    python_venv <- path(folder.venv, "bin", python_exe)
  }
  
  if (!file.exists(python_venv)) {
    stop("Python not found in path: ", python_venv)
  }
  
  print(paste("Putting PDB files in:", file))
  args <- paste(file)
  #browser()
  system2(python_venv, c(script_path, args = args))
  
}


# COMPUTES RELATIVE SASA FOR EACH RESIDUE USING SHRAKE RUPLEY.
# NORMALIZES BY MAX AMINO ACID SASA FROM TIEN ET AL. 2013.
frac_asa_py <- function(pdb, folder_in, folder_out, 
                        script_path = "Frac_ASA.py",
                        venv = ".venv",
                        file_temp = "REL_SASA.dat",
                        verbose = TRUE){
  
  pdb <- toupper(pdb)
  pdb_path <- path(folder_in, pdb, ext = "pdb")
  if (!dir.exists(folder_out)) dir_create(folder_out)
  asa_path <- path(folder_out, pdb, ext = "dat")
  args <- paste(pdb_path, asa_path)
  #browser()
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    python_venv <- path(venv, "Scripts", "python3.exe")
  } else {
    python_venv <- path(venv, "bin", "python3")
  }
  if (!file.exists(python_venv)) {
    stop("Python not found in path: ", python_venv)
  }
  if (verbose) print(paste("Getting ASA for: ", pdb))
  system2(python_venv, args = c(script_path, args))
  if (!file.exists(file_temp)) {
    stop("Some problem with path: ", file_temp)
  }
  #browser()
  file_out <- paste0(pdb,".dat")
  file.rename(file_temp, file_out)
  file_move(file_out, folder_out)
  res <- read_table(path(folder_out, file_out), col_names = c("Resnum", "Resname", "ASA"), 
                    show_col_types = FALSE) 
  res <- res |> mutate(PDB_ID = pdb) |> relocate(PDB_ID)
  #browser()
  return(res)  
}

# CLEANS PDB ATOM RECORDS BY FILTERING OUT HYDROGENS, HETEROATOMS AND UNWANTED RESIDUES.
# RENAMES SPECIFIC RESIDUE/TAG LABELS (E.G., HSD/HSE TO HIS, OT1 TO O, OT2 TO OXT).
# ADJUSTS CHAIN/RESIDUE COLUMNS (LIKE THAT INVOLVED IN OCCUPANCIES).
clean_pdb2 = function(name_file, name_new_file = "temp1.cln", i = 0){
  
  #browser()
  #i = 1 # i = 0 before
  #l1 = c("B", "2", "L")
  l1 = c("B", "C", "D", "2", "L")
  l2 = c("A","B","C","D","E","F","G","H","I")
  l3 = c("A","1","U")
  l4 = c("HOH","PMS","FOR","ALK","ANI")
  t = FALSE
  new_file = list()
  #name_new_file = "temp1.cln"
  con = file(name_new_file,"w")
  file = readLines(name_file)
  for(line in file){
    if(startsWith(line,"ATOM")){
      if((substr(line,13,13)!="H")&&(substr(line,14,14)!="H")&&(substr(line,14,16)!="OXT")&&(substr(line,14,15)!="2H")&&(substr(line,14,15)!="3H")&&(substr(line,14,14)!="D")&&(substr(line,13,13)!="D")){
        if(!(substr(line,17,17) %in% l1) && !(substr(line,27,27) %in% l2)){
          if(substr(line,17,17) %in% l3){
            line=paste(substr(line,1,16),substr(line,18,100))
          }
          for(item in l4){
            if (item %in% line){
              t=TRUE
            }
          }
          #browser()
          if(i==0){
            line=paste(substr(line,1,24)," 1",substr(line,27,100),sep="")
            i = 1
          }
          if("HSD" %in% line){
            line = gsub("HSD","HIS",line)
          }
          if("HSE" %in% line){
            line = gsub("HSE","HIS",line)
          }
          if("OT1" %in% line){
            line = gsub("OT1","O",line)
          }
          if("OT2" %in% line){
            line = gsub("OT2","OXT",line)
          }
          if("CD ILE" %in% line){
            line = gsub("CD ILE","CD1 ILE",line)
          }
          if (t==FALSE){
            writeLines(line,con)
          }
        }
      }
    }
  }
  writeLines("END",con)
  close(con)
}

# ALIGN TWO PDB SEQUENCES INTO A MULTIPLE SEQUENCE ALIGNMENT (DEFAULT CLUSTALW).
# OPTIONALLY SAVES ALIGNMENT PDF AND RETURNS THE MSA AAStringSet OBJECT.
align_pdbs <- function(pdb1, pdb2, method = "ClustalW",
                       folder_mca = "MCA", save_align = FALSE){
  
  #browser()
  PDB_ID1 <- attr(pdb1, "PDB")
  PDB_ID2 <- attr(pdb2, "PDB")
  seq1 <- pdb1 |> pdbseq() |> paste0(collapse = "")
  seq2 <- pdb2 |> pdbseq() |> paste0(collapse = "")
  aaligned_msa <- AAStringSet(c(seq1, seq2)) |>  msa(method = method) 
  #aaligned <- aaligned_msa |> as("AAStringSet")
  if (save_align){
    if (!dir.exists(folder_mca)) dir_create(folder_mca)
    file <- path(paste0(PDB_ID1,"_exp"), ext = "pdf")
    msaPrettyPrint(aaligned_msa, askForOverwrite=FALSE, output="pdf", file = file)
    file_move(file, folder_mca)
    file <- path(paste0(PDB_ID1,"_exp"), ext = "tex")
  file.remove(file)
  }
  return(aaligned_msa)
}

# GENERATES MSA PDF FOR GIVEN PDB_ID AND MOVES IT TO folder_mca.
# REMOVES TEMPORARY TEX FILE CREATED BY msaPrettyPrint.
print_alignment <- function(PDB_ID, aaligned_msa, folder_mca, suffix = ""){

  if (!dir.exists(folder_mca)) dir_create(folder_mca)
  file <- path(paste0(PDB_ID, suffix), ext = "pdf")
  msaPrettyPrint(aaligned_msa, askForOverwrite=FALSE, output="pdf", file = file)
  file_move(file, folder_mca)
  file <- path(paste0(PDB_ID, suffix), ext = "tex")
  file.remove(file)
}

# ALIGNS TWO PDB SEQUENCES, IDENTIFIES CONSENSUS REGION, AND TRIMS ATOMS ACCORDINGLY.
# VALIDATES MATCHED RESIDUES AND DIMENSIONS; RETURNS FALSE ON MISMATCH WITH OPTIONAL LOGGING.
# OPTIONALLY SAVES TRIMMED PDB FILES TO SPECIFIED DIRECTORIES.
align_and_save_pdb2 <- function(pdb1, pdb2, folder1, folder2, method = "ClustalW",
                             folder_mca = "MCA", log = FALSE, save = FALSE){
  
  #browser()
  PDB_ID1 <- attr(pdb1, "PDB")
  PDB_ID2 <- attr(pdb2, "PDB")
  #print(paste("ALIGNING AND SAVE", PDB_ID1))
  seq1 <- pdb1 |> pdbseq() |> paste0(collapse = "")
  seq2 <- pdb2 |> pdbseq() |> paste0(collapse = "")
  aaligned_msa <- AAStringSet(c(seq1, seq2)) |>  msa(method = method) 
  aaligned <- aaligned_msa |> as("AAStringSet")
  align1 <- as.character(aaligned[1])
  align2 <- as.character(aaligned[2])
  cons <- str_split(align1,"")[[1]] == str_split(align2,"")[[1]]
  alignc <- str_split(align1,"")[[1]][cons] |> paste(collapse = "")
  resn1 <- pdb1$atom$resno |> unique()
  resn2 <- pdb2$atom$resno |> unique()
  res1.loc <- str_locate(seq1, alignc) |> as.vector()
  res2.loc <- str_locate(seq2, alignc) |> as.vector()
  if (sum(is.na(res1.loc))==0){
    ids1 <- resn1[res1.loc[1]:res1.loc[2]]
    inds <- pdb1 |> atom.select(resno = ids1)
    pdb1.trim <- pdb1 |>  trim.pdb(inds = inds)
  } else {
    print(paste("WARNING: consensus delimitation failed in", PDB_ID1))
    if (log) print_alignment(PDB_ID1, aaligned_msa, folder_mca, suffix = "_csm")
    return(FALSE)
  }
  if (sum(is.na(res2.loc))==0){
    ids2 <- resn2[res2.loc[1]:res2.loc[2]]
    inds <- pdb2 |> atom.select(resno = ids2)
    pdb2.trim <- pdb2 |>  trim.pdb(inds = inds)
  } else {
    print(paste("WARNING: consensus delimitation failed in", PDB_ID1))
    if (log) print_alignment(PDB_ID1, aaligned_msa, folder_mca, suffix = "_exp")
    return(FALSE)
  }
  #browser()
  if (dim(pdb1.trim$atom)[1] != dim(pdb2.trim$atom)[1]){
    #browser()
    print(paste("WARNING: EXP and CSM sizes after aligning are different in", PDB_ID1))
    #if (log) print_alignment(PDB_ID1, aaligned_msa, folder_mca, suffix = "_size")
    return(FALSE)
  }
  if (sum(pdb1.trim$atom$resid != pdb2.trim$atom$resid)){
    print(paste("WARNING: residues matching are different in", PDB_ID1))
    #if (log) print_alignment(PDB_ID1, aaligned_msa, folder_mca, suffix = "_res")
    return(FALSE)
  }
  if (save){
    bio3d::write.pdb(pdb2.trim, file = path(folder2, PDB_ID2, ext = "pdb"))
    bio3d::write.pdb(pdb1.trim, file = path(folder1, PDB_ID1, ext = "pdb"))
  }
  if (log) print_alignment(PDB_ID1, aaligned_msa, folder_mca, suffix = "_ver")
  return(TRUE)
  #browser()
}

# SAVE ALIGNED PDBs
save_pdb_aligned <- function(alignment, pdb, folder){
  
  browser()
  align1 <- as.character(alignment[1])
  align2 <- as.character(alignment[2])
  cons <- str_split(align1,"")[[1]] == str_split(align2,"")[[1]]
  resn <- pdb$atom$resno |> unique()
  inds = pdb |> atom.select("protein")
  pdb <- pdb |>  trim.pdb(inds = inds)
}

# READ SRF FILES (OS BY ATOM)
read_prot2 = function(file, verbose=FALSE){
  if(verbose) print(file)
  #browser()
  dado = read_fwf(file,show_col_types = FALSE)
  dado = filter(dado, X1 == "INF")
  dado$X7 = NULL
  dado = rename(dado, INF = X1, ATOM = X2, NUMBER_POINTS = X3, AREA = X4, RAYLENGTH = X5, DISTANCE = X6)
  dado$NUMBER_POINTS = gsub("\\s\\pts","", dado$NUMBER_POINTS)
  dado$AREA = gsub("\\s\\A2","", dado$AREA)
  dado$RAYLENGTH = gsub("\\s\\Rlen","", dado$RAYLENGTH)
  dado$NUMBER_POINTS = as.integer(dado$NUMBER_POINTS)
  dado$AREA = as.double(dado$AREA)
  dado$RAYLENGTH = as.double(dado$RAYLENGTH)
  dado$DISTANCE = as.double(dado$DISTANCE)
  dado$INF = NULL
  return(dado)
}

# CALCULATES MEAN B-FACTOR PER RESIDUE FROM ATOM TABLE.
get_mean_b_factor = function(atom_tab, k = 3) {
  
  res = atom_tab |> as.tibble() |>  select(resid, chain, resno, b) |>
        group_by(resno, resid, chain) |> summarise(b = mean(b), .groups = "drop") |> 
        mutate(b = round(b, k)) |> 
        rename(Resnum = resno, Resname = resid, Chain = chain) |>
        ungroup()
  #browser()
  return(res)
}

general_plot <- function(m, var1, var2, var3, var4, labs, size = 0.1, method = NULL, 
                         ylim = NULL, brk = NULL, grd = 0, color.pt = "gray80", color.smo = "black",
                         type = "point", facet_grid = FALSE, pal = "Dark2"){
  
  #browser()
  if(type != "barplot"){
    if(is.null(var3)){
      p <- m |> ggplot(aes(x = .data[[var1]], y = .data[[var2]]))
    } else {
      p <- m |> ggplot(aes(x = .data[[var1]], y = .data[[var2]], color = .data[[var3]]))
    }
  }
  
  if (type == "point"){
    if(!is.null(var4)){
      p <- p + geom_point(aes(alpha = .data[[var4]]), size = size)
    } else {
      p <- p + geom_point(size = size, color = color.pt)
    }
  }
  
  if(type == "violin"){
    p <- p + geom_violin(aes(fill = .data[[var3]]), size = size, trim = TRUE) #, color = color.pt)
    # p <- p + stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", 
    #                       size=1, color=color.smo)
    #p <- p + geom_violin(aes(fill = size = size, color = color.pt)
  }
  #browser()
  if(type == "barplot"){
    
    df_counts <- osp.short %>%
      group_by(across(all_of(var1)), across(all_of(var2))) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(across(all_of(var1))) %>%
      mutate(prop = n / sum(n),
             label = paste0(!!sym(var2), "  ", percent(prop, accuracy = 0.1))) %>%
      ungroup()
    #browser()
    n_class <- (df_counts |> select(all_of(var2)) |> distinct() |> dim())[1]
    cores <- brewer.pal(n_class, pal)
    
    p <- ggplot(m, aes(x = .data[[var1]], fill = .data[[var2]])) +
      geom_bar(position = "fill") +
      geom_text(data = df_counts,
                aes(x = .data[[var1]], y = prop, label = label),
                position = position_fill(vjust = 0.5),
                color = "white",
                size = 8,
                fontface = "bold") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_fill_manual(values = cores)

  }
  
  if(!is.null(method)){
    p <- p + geom_smooth(method = method, se = TRUE, color = color.smo)
  }
  if(!is.null(ylim)){
    p <- p + ylim(ylim[1], ylim[2])
  }
  if(!is.null(brk)){
    p <- p + scale_x_continuous(breaks = brk)
  }
  p <- p + labs(
    x = labs[1],
    y = labs[2],
  )
  if(!is.null(var3)){
    if (facet_grid){
      p <- p + facet_grid(~ .data[[var3]])
    }
  }
  p <- p + theme_minimal()
  p <- p + theme(legend.position="none",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 #axis.title.y = element_text(size = 14, color = "black",
                 #                           vjust = 1.1, hjust = 0.5), #vjust = 0.9
                 axis.text.x = element_text(face="bold", size = 20, color = "black"),
                 axis.text.y = element_blank() #,
                 #axis.text.y = element_text(size = 14, color = "black") #,
                 #plot.margin = margin(t =20, r = 5, b = 20, l = 5)
                 ) #,
                 #strip.text.x = element_text(size = 12, color = "black")) 

  
  if(!is.null(grd)){
    if(grd == 0){
      p <- p + theme(panel.grid = element_blank())
    } else {
      p <- p + theme(panel.grid.minor = element_blank())
    }
  }
  
  return(p)
  
  
}


boot_fun <- function(split, fun = var, form = c("diff", "ratio")){
  
  #browser()
  d <- analysis(split)  # extrai a amostra bootstrap
  varA <- fun(d[[1]])
  varB <- fun(d[[2]])
  
  if (form[1] == "diff"){
    res <- varA - varB
  }
  if (form[1] == "ratio"){
    res <- varA/varB
  }

  return(res)
  
}

# CALCULATES SQUARED DEVIATIONS OF OSP.EXP AND OSP.CSM FROM VALUES IN V.
add_stat_square_diff <- function(m, v, k=3){
  
  #browser()
  m <- m |> mutate(OSP.exp.sqdev = round((OSP.exp-v[1])^2, k)) |> 
            relocate(OSP.exp.sqdev, .after = OSP.exp) |> 
            mutate(OSP.csm.sqdev = round((OSP.csm-v[2])^2, k)) |> 
            relocate(OSP.csm.sqdev, .after = OSP.csm) |> 
            mutate(OSP.sqdiff = round((OSP.csm.sqdev - OSP.exp.sqdev), k)) |> 
            relocate(OSP.sqdiff, .after = OSP.diff)
  return(m)
}

# IDENTIFIES INFLUENTIAL RESIDUES BY ITERATIVELY COMPARING DIFF VALUES USING FUN.
set_influential <- function(m, v, fun, cols, rev, k = 3, flag = "INFLUENTIAL", verbose = FALSE){
  
  #browser()
  osp.type <- cols[1]
  osp.diff <- cols[2]
  osp.infl <- cols[3]
  if (rev){
    m.order <- m |> arrange(.data[[osp.diff]])
    v = rev(v)
  } else {
    m.order <- m |> arrange(desc(.data[[osp.diff]]))
  }
  if (verbose) print(v)
  ids <- c()
  while(v[2] > v[1]){
    ids <- c(ids, m.order$Resid[1])
    #m$Influential_score[id] <- flag
    # m <- m |> mutate({{col}} : = replace({{col}}, id, flag)
    m.order <- m.order |> slice(-1)
    var <- m.order |> select(all_of(osp.type)) |> unlist(use.names = FALSE)
    r <- fun(var)
    v[2] <- round(r, k)
    if (verbose) print(paste(v[1], v[2], r))
    #browser()
    if (dim(m.order)[1] < 1) {
      print(paste("WARNING: all is influential in", m$PDB_ID[1], "this is weird"))
      break
    }
  }
  #browser()
  m <- m |> mutate({{osp.infl}} := replace(.data[[osp.infl]], ids, flag))
  return(m)
}

# FLAGS INFLUENTIAL RESIDUES BASED ON COMPARISON OF OSP VALUES AND SQDIFF.
# CHOOSES OSP.CSM OR OSP.EXP COLUMN DEPENDING ON V COMPARISON.
# CALLS set_influential WITH REVERSE FLAG AS NEEDED TO MARK RESIDUES.
set_influential_res <- function(m, v, fun = sd, k = 3, flag = "INFLUENTIAL", 
                                   cols = c("OSP.sqdiff", "Influential_score")){
  
  v <- round(v, k)
  if(v[2] == v[1]) {
    return(m)
  }
  if (v[2] > v[1]){
    cols = c("OSP.csm", cols)
    m <- set_influential(m, v, fun, cols = cols, rev = FALSE, k = k , flag = flag)
    return(m)
  } else {
    cols = c("OSP.exp", cols)
    m <- set_influential(m, v, fun, cols = cols, rev = TRUE, k = k , flag = flag)
    return(m)
  }
}

# DEFINE LABEL ACCORDING TO DIFFERENCES IN METRICS V IN CSM AND EXP
set_pack_type <- function(m, v, k = 3){
  
  v <- round(v, k)
  if(v[2] == v[1]) {
    m$Packtype = "CSM == EXP"
    return(m)
  }
  if (v[2] > v[1]){
    m$Packtype = "CSM > EXP"
    return(m)
  } else {
    m$Packtype = "CSM < EXP"
    return(m)
  }
}

# COMPUTES INFLUENCE SCORES FOR OSP EXP AND CSM BASED ON MEAN AND SD VALUES.
add_stat_infl_sd <- function(m, v, k=3){
  
  #browser()
  mean.exp = v[1]
  mean.csm = v[2]
  sd.exp = v[3]
  sd.csm = v[4]
  m <- m |> mutate(OSP.exp.infl = round( ((OSP.exp - mean.exp)^2 - sd.exp^2)/(2*sd.exp) , k)) |> 
            relocate(OSP.exp.infl, .after = OSP.exp.sqdev) |> 
            mutate(OSP.csm.infl = round( ((OSP.csm - mean.csm)^2 - sd.csm^2)/(2*sd.csm) , k)) |> 
            relocate(OSP.csm.infl, .after = OSP.csm.sqdev) |> 
            mutate(OSP.infl.diff = round((OSP.csm.infl - OSP.exp.infl), k)) |> 
            relocate(OSP.infl.diff, .after = OSP.sqdiff) |> 
            mutate(OSP.infl.sum = round((OSP.csm.infl + OSP.exp.infl), k)) |>  
            relocate(OSP.infl.sum, .after = OSP.infl.diff)
  return(m)
}




