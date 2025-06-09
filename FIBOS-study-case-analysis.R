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
conflict_prefer("path", "fs")


# SOURCE OF AUXILIARY FUNCTIONS
source("FIBOS-case-study-expanded-more-fun.R")

# MAIN DATASET ENTRIES
if(0){
  if(0){
    folder <- "data" 
    if (!dir.exists(folder)) dir_create(folder)
    filename <- "FIBOS_extended-dataset-all-v1.csv" 
    file <- path(folder,filename)
    db.extend <- read.csv(file)
  }
  
  if(0){
    folder <- "data" 
    if (!dir.exists(folder)) dir_create(folder)
    filename <- "FIBOS_strict-dataset-all-v1.csv" 
    file <- path(folder,filename)
    db.strict <- read.csv(file)
  }
  
  if(0){
    folder <- "data" 
    if (!dir.exists(folder)) dir_create(folder)
    filename <- "FIBOS_strict-dataset-clean-v1.csv" 
    file <- path(folder,filename)
    db.strict.clean <- read.csv(file)
  }
  
  if(0){
    folder <- "data" 
    if (!dir.exists(folder)) dir_create(folder)
    filename <- "FIBOS_by_residue-dataset-all-v1.csv" 
    file <- path(folder,filename)
    osp.res <- read.csv(file)
  }
}


# SOME USEFUL STATISTICS
if(0){
  
  x <- db.strict$osp.csm.mean
  y <- db.strict$osp.exp.mean
  wilcox.test(x, y, paired = T, conf.int = T)
  #ks.test(x, y)
  cohen.d(x, y, paired = T)

  x <- db.strict$osp.csm.sd
  y <- db.strict$osp.exp.sd
  wilcox.test(x, y, paired = T, conf.int = T)
  #ks.test(x, y)
  cohen.d(x, y, paired = T)
  
  x <- db.strict.clean$osp.csm.mean
  y <- db.strict.clean$osp.exp.mean
  wilcox.test(x, y, paired = T, conf.int = T)
  #ks.test(x, y)
  
  x <- db.strict.clean$osp.csm.sd
  y <- db.strict.clean$osp.exp.sd
  wilcox.test(x, y, paired = T, conf.int = T)
  #ks.test(x, y)
  cohen.d(x, y, paired = T)
}

# BOOTSTRAP FOR MEAN AND SD CONFIDENCE INTERVALS
if(0){

  ### BOOTSTRAP MEAN ###
  wilcox.test(db.strict$osp.csm.mean, db.strict$osp.exp.mean, paired = T, conf.int = T)
  obs_ref <- mean(db.strict$osp.csm.mean) - mean(db.strict$osp.exp.mean)
  boot_split <- db.strict |> select(osp.csm.mean, osp.exp.mean) |> bootstraps(times = 1000)
  boot_vals <- boot_split$splits |> map_dbl(\(x) boot_fun(x, fun = mean, form = "diff"))
  quantile(boot_vals, c(0.025, 0.975))

  ### BOOTSTRAP SD ###
  wilcox.test(db.strict$osp.csm.sd, db.strict$osp.exp.sd, paired = T, conf.int = T)
  #obs_ref <- median(db.strict$osp.csm.sd) - median(db.strict$osp.exp.sd)
  obs_ref <- mean(db.strict$osp.csm.sd) - mean(db.strict$osp.exp.sd)
  boot_split <- db.strict |> select(osp.csm.sd, osp.exp.sd) |> bootstraps(times = 1000)
  boot_vals <- boot_split$splits |> map_dbl(\(x) boot_fun(x, fun = mean, form = "diff"))
  quantile(boot_vals, c(0.025, 0.975))

  ### BOOTSTRAP MEAN CLEAN ###
  wilcox.test(db.strict.clean$osp.csm.mean, db.strict.clean$osp.exp.mean, paired = T, conf.int = T)
  obs_ref <- mean(db.strict.clean$osp.csm.mean) - mean(db.strict.clean$osp.exp.mean)
  boot_split <- db.strict.clean |> select(osp.csm.mean, osp.exp.mean) |> bootstraps(times = 1000)
  boot_vals <- boot_split$splits |> map_dbl(\(x) boot_fun(x, fun = mean, form = "diff"))
  quantile(boot_vals, c(0.025, 0.975))

  ### BOOTSTRAP SD CLEAN ###
  wilcox.test(db.strict.clean$osp.csm.sd, db.strict.clean$osp.exp.sd, paired = T, conf.int = T)
  obs_ref <- mean(db.strict.clean$osp.csm.sd) - mean(db.strict.clean$osp.exp.sd)
  boot_split <- db.strict.clean |> select(osp.csm.sd, osp.exp.sd) |> bootstraps(times = 1000)
  boot_vals <- boot_split$splits |> map_dbl(\(x) boot_fun(x, fun = mean, form = "diff"))
  quantile(boot_vals, c(0.025, 0.975))
}

# RUN TIME FIBOS AND VORONOI BENCHMARK
if(0){
  fibos.bench <- c(17024.72, 16961.94, 16926.32, 16838.16)
  voronoia.bench <- c(4701.701, 4700.120, 4699.052, 4705.802)
  x <- fibos.bench/60/60
  round(x,3)
  round(c(mean(x),sd(x)),3)
  x <- voronoia.bench/60/60
  round(x,3)
  round(c(mean(x),sd(x)),3)
}


# MAKES PLOT OF SI FIGURES S6A AND S6B
if(0){

  stat <- db.strict |> select(PDB_ID, starts_with("osp")) 
  par.lines <- tibble(
    mean = c(mean(stat$osp.exp.mean),mean(stat$osp.csm.mean)),
    sd =  c(mean(stat$osp.exp.sd),mean(stat$osp.csm.sd)),
    mad =  c(mean(stat$osp.exp.mad),mean(stat$osp.csm.mad)),
    type = c("solid", "81"),
    color = c("#F8766D","#00BFC4")
  )
  tag_pos = c(0.2, 1.02)
  
  p1a <- stat |> 
    #rename(`Experimental PDB` = osp.exp.mean, `Predicted AF` = osp.csm.mean) |>
    rename(`EXP` = osp.exp.mean, `CSM` = osp.csm.mean) |>
    select(!PDB_ID,) |> select(!contains("mad")) |> 
    pivot_longer(!ends_with("sd"), names_to = "TYPE:",
                 values_to = "mean(OSP)") |> 
    ggplot(aes(x = `mean(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity",  linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) + 
    scale_x_continuous(breaks = seq(0.24, 0.48, by = 0.04)) +
    labs(tag = expression(bold("A)") ~ " MEAN - ALL")) +
    #geom_vline(xintercept = c(0.28,0.36)) +
    theme_classic() + 
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          #axis.title.x = element_text(vjust = -1.5),
          axis.title.x = element_blank(),
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="bold", size=12),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
  #print(p1)
  
  p1b <- stat |> 
    rename(`EXP` = osp.exp.sd, `CSM` = osp.csm.sd) |>
  #p2 <- stat |> rename(`Experimental PDB` = osp.exp.mad, `Predicted AF` = osp.csm.mad) |>
    select(!PDB_ID) |> select(!contains("mad")) |> 
    #select(!PDB_ID) |> select(!contains("sd")) |> 
    pivot_longer(!ends_with("mean"), names_to = "TYPE:",
                 values_to = "sd(OSP)") |>
    ggplot(aes(x = `sd(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) +
    #xlim(0.12, 0.19) + 
    scale_x_continuous(breaks = seq(0.12, 0.19, by = 0.01)) +
    #xlim(0.12, 0.22) + 
    coord_cartesian(xlim = c(0.12, 0.19)) + 
    labs(tag = expression(bold("B)") ~ " SD - ALL")) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          #axis.title.x = element_text(vjust = -1.5),
          axis.title.x = element_blank(),
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="bold", size=12),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
}

# MAKES PLOT OF SI FIGURES S6C AND S6D
if(0){
  tag_pos = c(0.3, 1.02)
  stat <- db.strict.clean |> select(PDB_ID, starts_with("osp")) 
  par.lines <- tibble(
    mean = c(mean(stat$osp.exp.mean),mean(stat$osp.csm.mean)),
    sd =  c(mean(stat$osp.exp.sd),mean(stat$osp.csm.sd)),
    mad =  c(mean(stat$osp.exp.mad),mean(stat$osp.csm.mad)),
    type = c("solid", "81"),
    color = c("#F8766D","#00BFC4")
  )
  p1c <- stat |> 
    rename(`EXP` = osp.exp.mean, `CSM` = osp.csm.mean) |>
    select(!PDB_ID,) |> select(!contains("mad")) |> 
    pivot_longer(!ends_with("sd"), names_to = "TYPE:",
                 values_to = "mean(OSP)") |> 
    ggplot(aes(x = `mean(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity",  linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) + 
    scale_x_continuous(breaks = seq(0.24, 0.48, by = 0.04)) +
    labs(tag = expression(bold("D)") ~ " MEAN - NON-INFLUENTIAL")) +
    #geom_vline(xintercept = c(0.28,0.36)) +
    theme_classic() + 
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust = -1.2, size = 11),
          #axis.title.x = element_blank(),
          plot.margin = margin(t = 20, r = 5, b = 20, l = 5), 
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="bold", size=12),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
  #print(p1)
  
  p1d <- stat |> 
    rename(`EXP` = osp.exp.sd, `CSM` = osp.csm.sd) |>
    #p2 <- stat |> rename(`Experimental PDB` = osp.exp.mad, `Predicted AF` = osp.csm.mad) |>
    select(!PDB_ID) |> select(!contains("mad")) |> 
    #select(!PDB_ID) |> select(!contains("sd")) |> 
    pivot_longer(!ends_with("mean"), names_to = "TYPE:",
                 values_to = "sd(OSP)") |>
    ggplot(aes(x = `sd(OSP)`, color = `TYPE:`, linetype = `TYPE:`)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.5) +
    scale_linetype_manual(values = par.lines$type) +
    scale_x_continuous(breaks = seq(0.12, 0.19, by = 0.01)) +
    #xlim(0.12, 0.19) + 
    coord_cartesian(xlim = c(0.12, 0.19)) + 
    labs(tag = expression(bold("D)") ~ " SD - NON-INFLUENTIAL")) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust = -1.2, size = 11),
          #axis.title.x = element_blank(),
          plot.margin = margin(t = 20, r = 5, b = 20, l = 5),
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="plain", size=12),
          legend.position = "none") +
    guides(linetype = guide_legend(keywidth = 3))
  #print(p1d)

  title <- "Protein Packing Density Distributions for CSM and EXP Models\n"
  subtitle <- " "
  pa2ab2 <- (p1a + p1b) / (p1c + p1d) +
    plot_annotation(
      title    = title,
      subtitle = subtitle, 
    ) +
    plot_layout(guides = "collect") &
    theme(
      plot.title    = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0.5,),
      legend.position = "bottom"
    )
  print(pa2ab2)
}

# MAKES PLOT OF SI FIGURE S8A
if(0){
  db.ana <- db.extend |> select(PDB_ID, !contains("path")) |> select(!contains(".id")) |> 
    select(!c(TITLE, n_SEQRES, n_SEQAF)) |> 
    mutate(R_FACTOR_M = rowMeans(across(R_FACTOR:R_FACTOR_W), na.rm = TRUE)) |> 
    relocate(R_FACTOR_M, .after = RESOLUTION) |> 
    select(!c(R_FACTOR, R_FACTOR_W))
    #select(!c(R_FACTOR, R_FACTOR_W, n_SEQDIFF))
  db.ana <- db.ana |> rename(strand = STRAND, helix = HELIX, resolution = RESOLUTION,
                             R.Factor.mean = R_FACTOR_M, size = n_SEQSTR, 
                             size_diff = n_SEQDIFF)
  db.ana <- db.ana |> left_join(res.stat, by = "PDB_ID") |> select(!ends_with("mad"))
  
  title <- "STRUCTURAL METRIC CORRELATIONS AT CHAIN LEVEL"
  subtitle <- paste("Extended Dataset - ",dim(db.ana)[1], "chains")
  corr <- cor(db.ana[-1], method = "pearson", use = "pairwise.complete.obs")
  corr.p <- cor_pmat(db.ana[-1], method = "pearson")
  pa2a <- ggcorrplot(corr, type = "lower", lab = "true", p.mat = corr.p, lab_size = 3.2)
  #print(pa2a)
  corr <- cor(db.ana[-1], method = "spearman", use = "pairwise.complete.obs")
  corr.p <- cor_pmat(db.ana[-1], method = "spearman")
  pa2b <- ggcorrplot(corr, type = "lower", lab = "true", p.mat = corr.p, lab_size = 3.2)
  #print(pa2b)
  pa2ab <- ggarrange(pa2a, pa2b, nrow = 1, ncol=2, 
                     labels=c("A) Pearson\n",
                              "B) Spearman\n"), 
                     common.legend = T,
                     legend = "bottom")
  #print(pa2ab)
  pa2ab <- annotate_figure(
    annotate_figure(pa2ab, top = text_grob(subtitle, color = "black",  face = "plain", 
                                           size = 12, vjust = 1)),
    top =  text_grob(title, color = "black",  face = "plain", size = 18, vjust = 0.7)
  )
  print(pa2ab)
}
 
# MAKE PLOT OF SI FIGURE S8B
if(0){ 
  db.ana <- db.strict |> select(PDB_ID, !contains("path")) |> select(!contains(".id")) |> 
    select(!c(TITLE, n_SEQRES, n_SEQAF)) |> select(!contains("mad")) |> 
    mutate(R_FACTOR_M = rowMeans(across(R_FACTOR:R_FACTOR_W), na.rm = TRUE)) |> 
    relocate(R_FACTOR_M, .after = RESOLUTION) |> 
    select(!c(R_FACTOR, R_FACTOR_W))
  #select(!c(R_FACTOR, R_FACTOR_W, n_SEQDIFF))
  
  db.ana <- db.ana |> rename(strand = STRAND, helix = HELIX, resolution = RESOLUTION,
                             R.Factor.mean = R_FACTOR_M, size = n_SEQSTR, 
                             size_diff = n_SEQDIFF)
  
  #db.ana <- db.ana |> left_join(res.stat, by = "PDB_ID") |> select(!ends_with("mad"))
  
  title <- "STRUCTURAL METRIC CORRELATIONS AT CHAIN LEVEL"
  subtitle <- paste("Strict Dataset - ",dim(db.ana)[1], "chains")
  
  corr <- cor(db.ana[-1], method = "pearson", use = "pairwise.complete.obs")
  corr.p <- cor_pmat(db.ana[-1], method = "pearson")
  pb2a <- ggcorrplot(corr, type = "lower", lab = "true", p.mat = corr.p, lab_size = 3.2)
  #print(pb2a)
  
  corr <- cor(db.ana[-1], method = "spearman", use = "pairwise.complete.obs")
  corr.p <- cor_pmat(db.ana[-1], method = "spearman")
  pb2b <- ggcorrplot(corr, type = "lower", lab = "true", p.mat = corr.p, lab_size = 3.2)
  #print(pb2b)
  
  pb2ab <- ggarrange(pb2a, pb2b, nrow = 1, ncol=2, 
                      labels=c("A) Pearson\n",
                               "B) Spearman\n"), 
                      common.legend = T,
                      legend = "bottom")
  #print(pa2ab)
  
  pb2ab <- annotate_figure(
    annotate_figure(pb2ab, top = text_grob(subtitle, color = "black",  face = "plain", 
                                            size = 12, vjust = 1)),
    top =  text_grob(title, color = "black",  face = "plain", size = 18, vjust = 0.7)
  )
  print(pb2ab)
  
}

# SEVERAL PLOTS INVOLVING COMPLEMENTARY METRICS
if(0){
  
  labels <- c(
    "OSP.csm" = "OSP.csm",
    "ASA.csm" = "ASA.csm",
    "b" = "B-FACTOR",
    "pLDDT" = "pLDDT",
    "chi1" = expression("abs("*Delta*"CHI1)"),
    "chi2" = expression("abs("*Delta*"CHI2)"),
    "chi3" = expression("abs("*Delta*"CHI3)"),
    "chi4" = expression("abs("*Delta*"CHI4)"),
    "chi5" = expression("abs("*Delta*"CHI5)"),
    "Position" = "Position",
    "dOSP"  = expression(Delta*"OSP"),
    "dASA" = expression(Delta*"ASA"),
    "OSP.diff"  = expression(Delta*"OSP"),
    "OSP.exp" = "OSP.exp",
    "Surface_exp" = "SURF_EXP",
    "Surface_csm" = "SURF_CSM",
    "Influential" = "INFLUENTIAL",
    "adOSP"  = expression("abs("*Delta*"OSP<0)")
    )
  
  x <- "Influential_if"
  y <- "OSP.csm"
  #y <- "OSP.diff"
  z <- NULL
  z <- x
  w <- NULL
  i = 1
  labs = labels[c(16,i)]
  pp1 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp1)
  
  y <- "ASA.csm"
  i = i + 1
  labs = labels[c(16,i)]
  pp2 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp2)
  
  y <- "b"
  i = i + 1
  labs = labels[c(16,i)]
  pp3 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp3)
  
  y <- "pLDDT"
  i = i + 1
  labs = labels[c(16,i)]
  pp4 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp4)
  
  y <- "chi1"
  i = i + 1
  labs = labels[c(16,i)]
  pp5 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp5)
  
  y <- "chi2"
  i = i + 1
  labs = labels[c(16,i)]
  pp6 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp6)
  
  y <- "chi3"
  i = i + 1
  labs = labels[c(16,i)]
  pp7<- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp7)
  
  y <- "chi4"
  i = i + 1
  labs = labels[c(16,i)]
  pp8 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp8)
  
  y <- "chi5"
  i = i + 1
  labs = labels[c(16,i)]
  pp9 <- general_plot(osp.res, x, y, z, w, labs, type = "violin")
  #print(pp9)
  
  #source("FIBOS-case-study-expanded-more-fun.R")
  y <- "Position.exp"
  i = i + 1
  i = 10
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  pp10 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot", pal = "Set1")
  #print(pp10)
  
  #source("FIBOS-case-study-expanded-more-fun.R")
  y <- "Position.csm"
  i = i + 1
  i = 11
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  pp11 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot", pal = "Set1")
  #print(pp11)
  
  #source("FIBOS-case-study-expanded-more-fun.R")
  y <- "SSE.exp"
  i = i + 1
  i = 12
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  pp12 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot")
  #print(pp12)
  
  #source("FIBOS-case-study-expanded-more-fun.R")
  y <- "SSE.csm"
  i = i + 1
  i = 13
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  pp13 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot")
  #print(pp13)
}

# MAKES PLOT OF SI FIGURE S7
if(0){
  v <- c(1.2,0,1.2,0.4)
  pp1 <- pp1 + theme(plot.margin = unit(v, "cm"))
  pp2 <- pp2 + theme(plot.margin = unit(v, "cm"))
  pp3 <- pp3 + theme(plot.margin = unit(v, "cm"))
  pp4 <- pp4 + theme(plot.margin = unit(v, "cm")) 
  ppa <- ggarrange(pp1, pp2, pp3, pp4, nrow = 2, ncol=2, labels=c("A)","B)","C)","D)"), 
                   common.legend = T, legend = "none") #, vjust = -0.4)
  print(ppa)
}

# MAKES PLOT OF MAIN TEXT FIGURE 1B (PART OF IT)
if(0){
  v <- c(0.8,0.0,0.5,0.0)
  pp11 <- pp11 + theme(plot.margin = unit(v, "cm"))
  pp13 <- pp13 + theme(plot.margin = unit(v, "cm")) 
  ppc1 <- ggarrange(pp11, pp13, nrow = 1, ncol=2, widths = c(1,1),
                   labels=c("A) CSM - Position","B) CSM - SSE"), 
                   common.legend = FALSE, legend = "none")
  print(ppc1)
}

# MAKES PLOT OF SI FIGURE S10
if(0){
  df <- osp.res |> filter(Position.exp == "CORE") |> group_by(PDB_ID) |> summarise(OSP.exp = mean(OSP.exp))
  y.stat <- t.test(df$OSP.exp, conf.level = 0.95)
  y.med <- y.stat$estimate |> round(3)
  y.int2 <- y.stat$conf.int[2] |> round(3)
  y.err <- (y.int2 - y.med) |> round(3)
  
  y.med2 <- mean(df$OSP.exp) |> round(2)
  y.sd <- sd(df$OSP.exp) |> round(2)
  y.mod <- 0.557
  
  tag_pos = c(0.4, 1.05)
  # title <-  paste("CORE OSP.exp residue distribution for", dim(df)[1], "chains\n")
  # subtitle <- " "
  pb8f <- df |> ggplot(aes(x = OSP.exp)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.5) +
    scale_x_continuous(breaks = seq(0.45, 0.60, by = 0.05)) + 
    coord_cartesian(xlim = c(0.45, 0.62)) + 
    geom_vline(xintercept = y.med2, linetype = "dotted") +
    labs(tag = expression(bold("A)") ~ " CORE OSP.EXP DISTRIBUTION")) +
    #geom_vline(xintercept = y.mod, linetype = "dashed") +
    labs() +
    annotate("text", x = y.mod, y = 27,
             label = y.mod,
             color = "black", size = 5, hjust = 0.1) +
    annotate("text", x = y.med2, y = 27,
             label = bquote(.(y.med2) %+-% .(y.sd)), 
             color = "black", size = 5, hjust = 1.1) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(vjust = -1.2, size = 16),
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="plain", size=16),
          plot.margin = margin(t = 20, r = 5, b = 20, l = 5),
          axis.ticks.y = element_blank())
  
  
  df <- osp.res.infl.if |> group_by(PDB_ID) |> summarise(n_INFLU = n())
  df <- db.strict |> left_join(df, by = "PDB_ID") |> relocate(n_INFLU, .after = n_SEQDIFF) |> 
        select(PDB_ID, n_SEQSTR, n_INFLU) |> mutate(n_INFLUp = round(n_INFLU/n_SEQSTR,3)*100)
  df <- df |> mutate(across(everything(), ~ replace_na(., 0)))
  
  y.med <- median(df$n_INFLUp, na.rm = T)
  y.med2 <- mean(df$n_INFLUp, na.rm = T) |> round(2)
  y.sd <- sd(df$n_INFLUp, na.rm = T) |> round(2)
  y.mod <- 1.2
  
  tag_pos = tag_pos #c(0.3, 1.01)
  title <-  paste("AUXILIARY PLOTS")
  title <- " "
  subtitle <- " "
  pb8g <- df |> ggplot(aes(x = n_INFLUp)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.5) +
    scale_x_continuous(breaks = seq(0, 20, by = 2)) +
    coord_cartesian(xlim = c(-1, 19)) + 
    geom_vline(xintercept = y.med, linetype = "dotted") +
    #geom_vline(xintercept = y.mod, linetype = "dashed") +
    #xlim(-2, 20) +
    labs(tag = expression(bold("B)") ~ " INFLUENTIAL RESIDUES DISTRIBUTION"),
         x = "Influential residues (%)") +
    annotate("text", x = y.mod, y = 0.26,
             label = paste(y.mod,"%"),
             color = "black", size = 5, hjust = 1.05) +
    annotate("text", x = y.med, y = 0.26,
             label = bquote(bold(.(y.med)) ~ "%"), #bquote(.(y.med2) %+-% .(y.sd)),
             color = "black", size = 5, hjust = -0.2) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(vjust = -1.2, size = 16),
          plot.tag.position = tag_pos,
          plot.tag = element_text(face="plain", size=16),
          plot.margin = margin(t =20, r = 5, b = 20, l = 5),
          axis.ticks.y = element_blank())
  #print(pb8g)
  
  pb8fg <- (pb8f + plot_spacer() + pb8g) +
    plot_annotation(
      title    = title,
      subtitle = subtitle, 
    ) +
    plot_layout(guides = "collect",
                widths = c(1, 0.2, 1)) &
    theme(
      plot.title    = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0.5, vjust =0)#,
     # legend.position = "bottom"
    )
  
  print(pb8fg)
}

# FORMATS IMAGE OF SI FIGURE S9
if(0){
  imgA <- image_read("Figs/PUB/1PZ4-pub-v1h.png")
  imgB <- image_read("Figs/PUB/1PZ4-pub-v2h.png")
  labelA <- image_annotate(imgA, "A) CSM", size = 50, gravity = "northwest",
                           location = "+20+20", color = "black",
                           font = "Arial_Bold")
  labelB <- image_annotate(imgB, "B) EXP", size = 50, gravity = "northwest",
                           location = "+20+20", color = "black",
                           font = "Arial_Bold")
  final <- image_append(c(labelA, labelB))
  image_write(final, "Figs/PUB/1PZ4-pub-v12h.png")
}

# READJUSTS PLOT OF MAIN TEXT FIGURE 1B
if(0){
  source("FIBOS-case-study-expanded-more-fun.R")
  y <- "Position.csm"
  i = i + 1
  i = 11
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  pp11 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot", pal = "Set1")

  title <- "CSM - Position"
  pp11 <- pp11 + labs(title = title) + 
          theme(plot.title = element_text(face="bold", size = 24, hjust = 0.5, vjust = -1))
  #print(pp11)
  
  ggsave("Figs/PUB/bar-pos-v1.png", plot = pp11,
         width = 6.857, height = 4.571, dpi = 350, units = "in")
  
  y <- "SSE.csm"
  i = i + 1
  i = 13
  labs = labels[c(16,i)]
  labs = c("Position (%)", "")
  
  pp13 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot")
  pp13 <- general_plot(osp.res, x, y, z, w, labs, type = "barplot", pal = "Set1")

  title <- "CSM - SSE"
  pp13 <- pp13 + labs(title = title) + 
    theme(plot.title = element_text(face="bold", size = 24, hjust = 0.5, vjust = -1))
  #print(pp13)
  
  ggsave("Figs/PUB/bar-pos-v2.png", plot = pp13,
         width = 6.857, height = 4.571, dpi = 350, units = "in")
}

# FORMATS ALL PANELS OF FIGURE 1 FROM THE MAIN TEXT

if(0){
  img1 <- image_read("Figs/PUB/1ubq-v1a.png")
  img2 <- image_read("Figs/PUB/1ubq-v2a.png")
  img3 <- image_read("Figs/PUB/bar-pos-v1.png")
  img4 <- image_read("Figs/PUB/bar-pos-v2.png")
  img5 <- image_read("Figs/PUB/1NG6-pub-v2e.png")
  img6 <- image_read("Figs/PUB/1NG6-pub-v1e.png")
  
  img1 <- image_background(img1, color = "white", flatten = TRUE)
  img2 <- image_background(img2, color = "white", flatten = TRUE)
  img3 <- image_background(img3, color = "white", flatten = TRUE)
  img4 <- image_background(img4, color = "white", flatten = TRUE)
  img5 <- image_background(img5, color = "white", flatten = TRUE)
  img6 <- image_background(img6, color = "white", flatten = TRUE)
  

  img.dim <- "600x400"
  img1 <- image_resize(img1, img.dim)
  img2 <- image_resize(img2, img.dim)
  img3 <- image_resize(img3, img.dim)
  img4 <- image_resize(img4, img.dim)
  img5 <- image_resize(img5, img.dim)
  img6 <- image_resize(img6, img.dim)
  
  img.size = 44
  img.gravity = "north" #"northwest" 
  img.loc = "+20+20"
  img.font = "Arial-Bold"
  img.color = "black"
  
  img.geo <- "0x75"

  img1 <- image_border(img1, color = "white", geometry = img.geo)
  #img2 <- image_border(img2, color = "white", geometry = img.geo)
  img3 <- image_border(img3, color = "white", geometry = img.geo)
  #img4 <- image_border(img4, color = "white", geometry = img.geo)
  img5 <- image_border(img5, color = "white", geometry = img.geo)
  #img6 <- image_border(img6, color = "white", geometry = img.geo)

  
  img1 <- image_annotate(img1, "A", size = img.size, gravity = img.gravity, 
                         location = img.loc, font = img.font, color = img.color)
  # img2 <- image_annotate(img2, "B) NORMALS", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)
  img3 <- image_annotate(img3, "B", size = img.size, gravity = img.gravity, 
                         location = img.loc, font = img.font, color = img.color)
  # img4 <- image_annotate(img4, "D) SSE", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)
  img5 <- image_annotate(img5, "C", size = img.size, gravity = img.gravity, 
                         location = img.loc, font = img.font, color = img.color)
  # img6 <- image_annotate(img6, "F) CSM", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)
  
  # img4 <- image_annotate(img4, "D)", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)
  # img5 <- image_annotate(img5, "E)", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)
  # img6 <- image_annotate(img6, "F)", size = img.size, gravity = img.gravity, 
  #                        location = img.loc, font = img.font, color = img.color)

  linha1 <- image_append(c(img1, img3, img5))
  linha2 <- image_append(c(img2, img4, img6))
  
  figura_final <- image_append(c(linha1, linha2), stack = TRUE)
  image_write(figura_final, "Figs/PUB/main-fig-v1.png")
}


