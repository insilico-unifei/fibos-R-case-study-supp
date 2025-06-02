###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# MAIN R SCRIPT FOR: images, graphs, tables, comparative validations used 
# in the article's supplementary material
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
# See more here:  https://github.com/insilico-unifei/fibos-R 

library(FIBOS)
library(tidyverse)
library(fs)
library(httr)
library(bio3d)
library(furrr)
library(tictoc)
library(ggpubr)
library(tictoc)
library(rgl)
library(htmlwidgets)
library(sphereTessellation)
library(pagedown)
library(plotly)
library(pracma)
library(geosphere)
library(rPDBapi)

# SOURCE OF AUXILIARY FUNCTIONS
source("FIBOS-supp-fun.R")

# DEFINE PARAMETERS FOR OS ON SPHERE, FOR 212 DOTS
if(0){
  dotden = 5 # dot density
  r = 1.9 # radius
  N = round(4*pi*r^2*dotden) # number of dots
  i = 1 
  ss = list()
}

# GENERATE CLASSICAL OS DOT DISTRIBUTION (212 DOTS)
if(0){
  ss[[i]] = list()
  res_c = GENUN(N=N)
  m_os = t(res_c$U)
  colnames(m_os) = c("x","y","z")
  m_os = m_os %>% as_tibble()
  f.os = m_os %>% dplyr::transmute(sum = x+y+z != 0) %>% as.vector() %>% unlist()
  m_os = m_os[f.os,]
  m_os = round(m_os,3)
  ss[[i]]$m_os = m_os; print(ss[[i]]$m_os)
  ss[[i]]$dotden.eff = dim(ss[[i]]$m_os)[1]/(4*pi*r^2)
}

# GENERATE NEW FIBOS DOT DISTRIBUTION (212 DOTS)
if(0){
  m_fb2 = fibonaccisphere3(n = dim(ss[[i]]$m_os)[1], r = 1, eps = 0.5)
  m_fb2 = m_fb2 %>% as_tibble()
  f.fb2 = m_fb2 %>% dplyr::transmute(sum = x+y+z != 0) %>% as.vector() %>% unlist()
  m_fb2 = m_fb2[f.fb2, ]
  ss[[i]]$m_fb2 = m_fb2
}

# DEFINE PARAMETERS FOR OS ON SPHERE, FOR 417 DOTS
if(0){
  dotden = 10
  r = 1.9
  N = round(4*pi*r^2*dotden)
  i = 2
}

# GENERATE CLASSICAL OS DOT DISTRIBUTION (417 DOTS)
if(0){
  ss[[i]] = list()
  res_c = GENUN(N=N)
  m_os = t(res_c$U)
  colnames(m_os) = c("x","y","z")
  m_os = m_os %>% as_tibble()
  f.os = m_os %>% dplyr::transmute(sum = x+y+z != 0) %>% as.vector() %>% unlist()
  m_os = m_os[f.os,]
  m_os = round(m_os,3)
  ss[[i]]$m_os = m_os; print(ss[[i]]$m_os)
  ss[[i]]$dotden.eff = dim(ss[[i]]$m_os)[1]/(4*pi*r^2)
}

# GENERATE NEW FIBOS DOT DISTRIBUTION (417 DOTS)
if(0){
  m_fb2 = fibonaccisphere3(n = dim(ss[[i]]$m_os)[1], r = 1, eps = 0.5)
  m_fb2 = m_fb2 %>% as_tibble()
  f.fb2 = m_fb2 %>% dplyr::transmute(sum = x+y+z != 0) %>% as.vector() %>% unlist()
  m_fb2 = m_fb2[f.fb2, ]
  ss[[i]]$m_fb2 = m_fb2
}

# GENERATE VORONOI TILING ON SPHERE (212 DOTS)
# MAKE IMAGENS OF SUPPLEMENTARY FIGURE S1a OF THE ARTCILE
#
# OBS 01: In RStudio, you may need reset it by pressing Ctrl+Shift+F10 or Cmd+Shift+F10 (Mac)
# between each Voronoi image generation
# OBS 02: Depending on the computer, it can be very slow
if(0){

  folder = "images" 
  if (!dir.exists(folder)) dir_create(folder)
  
  i = 1
  # VORONOI FOR OS (212 DOTS)
  ss[[i]]$m_os.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_os))
  plotVoronoiOnSphere(ss[[i]]$m_os.vor, edges = T, specular = "black")
  filename = "os_vor"
  file_html = path(folder, filename , ext = "html")
  saveWidget(rglwidget(), file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
}

if(0){

  # VORONOI FOR OS ROTATED BY 90 AT X AXIS (212 DOTS)
  ss[[i]]$m_os.rot = rotate3d(as.matrix(ss[[i]]$m_os), pi/2, 1, 0, 0)
  ss[[i]]$m_os.vor.rot = VoronoiOnSphere(as.matrix(ss[[i]]$m_os.rot))
  plotVoronoiOnSphere(ss[[i]]$m_os.vor.rot, edges = T, specular = "black")
  filename = "os_vor_rot"
  file_html = path(folder, filename , ext = "html")
  saveWidget(rglwidget(), file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
}

if(0){
  # VORONOI FOR FIBOS (212 DOTS)
  ss[[i]]$m_fb2.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_fb2))
  plotVoronoiOnSphere(ss[[i]]$m_fb2.vor, edges = T, specular = "black")
  filename = "fb2_vor"
  file_html = path(folder, filename , ext = "html")
  saveWidget(rglwidget(), file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
}

if(0){
  # VORONOI FOR FIBOS ROTATED BY 90 AT X AXIS (212 DOTS)
  ss[[i]]$m_fb2.rot = rotate3d(as.matrix(ss[[i]]$m_fb2), pi/2, 1, 0, 0)
  ss[[i]]$m_fb2.vor.rot = VoronoiOnSphere(as.matrix(ss[[i]]$m_fb2.rot))
  plotVoronoiOnSphere(ss[[i]]$m_fb2.vor.rot, edges = T, specular = "black")
  filename = "fb2_vor_rot"
  file_html = path(folder, filename , ext = "html")
  saveWidget(rglwidget(), file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
}

# MAKE IMAGENS OF SUPPLEMENTARY FIGURE S1b OF THE ARTCILE
if(0){
  
  folder = "images" 
  if (!dir.exists(folder)) dir_create(folder)
  
  i = 2
  
  m = as_tibble(ss[[i]]$m_os)
  colnames(m) = c("x","y","z")
  col = abs(m$z)
  scene = list(camera = list(eye = list(x = -1.5, y = 1.5, z = 0)))
  p3 = plot_ly()
  p3 = p3 %>% add_trace(
    data = m,
    x = ~`x`,
    y = ~`y`,
    z = ~`z`,
    color = col,
    type = "scatter3d",
    mode = 'markers',
    marker = list(size = 6),
    showlegend = F
  ) %>% layout(title = list(text = "<b>OS</b>", x = 0.23, y = 0.85, font = list(size = 24)), 
               scene = scene) %>% hide_colorbar()
  #print(p3)
  
  filename = "OS_scatter_3d"
  file_html = path(folder, filename , ext = "html")
  saveWidget(p3, file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
  
  m = as_tibble(ss[[i]]$m_fb2)
  colnames(m) = c("x","y","z")
  col = abs(m$z)
  scene = list(camera = list(eye = list(x = -1.5, y = 1.5, z = 0)))
  p4 = plot_ly()
  p4 = p4 %>% add_trace(
    data = m,
    x = ~`x`,
    y = ~`y`,
    z = ~`z`,
    color = col,
    type = "scatter3d",
    mode = 'markers',
    marker = list(size = 6),
    showlegend = F
  ) %>% layout(title = list(text = "<b>FIBOS</b>", x = 0.23, y = 0.85, font = list(size = 24)), 
               scene = scene) %>% hide_colorbar()
  #print(p4)
  
  filename = "FIBOS_scatter_3d"
  file_html = path(folder, filename , ext = "html")
  saveWidget(p4, file_html)
  file_png = path(folder, filename, ext = "png")
  chrome_print(file_html, output = file_png, format = "png", timeout = 600)
  dir_delete(paste0(path(folder, filename),"_files"))
  file_delete(file_html)
  
}

# CALCULATE THE AREAS OF VORONOI CELLS ON OS AND FIBOS SPHERES
if(0){
  
  i = 1 # with 212 dots
  
  if(is.null(ss[[i]]$m_os.vor)) ss[[i]]$m_os.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_os))
  ss[[i]]$os.vor.area = ss[[i]]$m_os.vor %>% map(Vor_Cell_Area) %>% unlist()
  
  if(is.null(ss[[i]]$m_fb2.vor)) ss[[i]]$m_fb2.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_fb2))
  ss[[i]]$fb2.vor.area = ss[[i]]$m_fb2.vor %>% map(Vor_Cell_Area)%>% unlist()
  
  i = 2  # with 417 dots
  
  if(is.null(ss[[i]]$m_os.vor)) ss[[i]]$m_os.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_os))
  ss[[i]]$os.vor.area = ss[[i]]$m_os.vor %>% map(Vor_Cell_Area) %>% unlist()
  
  if(is.null(ss[[i]]$m_fb2.vor)) ss[[i]]$m_fb2.vor = VoronoiOnSphere(as.matrix(ss[[i]]$m_fb2))
  ss[[i]]$fb2.vor.area = ss[[i]]$m_fb2.vor %>% map(Vor_Cell_Area)%>% unlist()
  
}

# MAKE IMAGEN OF SUPPLEMENTARY FIGURE S2a OF THE ARTCILE
if(0){

  dotarea = c()
  dn = c()
  i = 1
  aux = ss[[i]]$os.vor.area*r^2
  d = paste0("Dots = ", length(aux))
  dotarea = c(dotarea, round(1/ss[[i]]$dotden.eff,3))
  df1 = tibble(Type = "OS", Area = aux, Dots = d)

  aux = ss[[i]]$fb2.vor.area*r^2
  d = paste0("Dots = ", length(aux))
  dn = c(dn, d)
  df2 = tibble(Type = "FIBOS", Area = aux, Dots = d)

  df = df1 %>% bind_rows(df2)
  
  i = 2
  aux = ss[[i]]$os.vor.area*r^2
  d = paste0("Dots = ", length(aux))
  dn = c(dn, d)
  dotarea = c(dotarea, round(1/ss[[i]]$dotden.eff,3))
  df1 = tibble(Type = "OS", Area = aux, Dots = d)
  
  aux = ss[[i]]$fb2.vor.area*r^2
  d = paste0("Dots = ", length(aux))
  df2 = tibble(Type = "FIBOS", Area = aux, Dots = d)
  
  df = df %>% bind_rows(df1) %>% bind_rows(df2)
  
  p5 = df %>% ggplot(aes(x = Area, fill = Type)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Dots, scales = "free") +
    ggtitle("VORONOI CELL AREA DISTRIBUTION") + #,
    theme_minimal() + 
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = rel(1.1)),
          axis.title.y = element_blank(),
          strip.text = element_text(size=12),
          legend.position = "bottom")
  
  
  df = tibble(label = dotarea,
              density = c(350, 2300),
              Dots = dn)
  p5 = p5 + geom_text(data = df,
                      mapping = aes(x = dotarea, y = density, label = dotarea),
                      inherit.aes=F,
                      size = 3,
                      nudge_x = c(0.003, 0.0015))
  
  print(p5)
  
}

# MAKE PLOT OF SUPPLEMENTARY FIGURE S3 OF THE ARTCILE
if(0){
  i = 1
  cap_cut = 0.75
  aux = ss[[i]]$m_fb2
  d = paste0("Dots = ", dim(aux)[1])
  
  dfx = tibble(Type = "FIBOS", Axes = "x", Dist = as.vector(dist(aux %>% filter(x >= cap_cut))), Dots = d)
  dfy = tibble(Type = "FIBOS", Axes = "y", Dist = as.vector(dist(aux %>% filter(y >= cap_cut))), Dots = d)
  dfz = tibble(Type = "FIBOS", Axes = "z", Dist = as.vector(dist(aux %>% filter(z >= cap_cut))), Dots = d)
  
  df = dfx %>% bind_rows(dfy) %>% bind_rows(dfz)
  
  aux = ss[[i]]$m_os
  dfx = tibble(Type = "OS", Axes = "x", Dist = as.vector(dist(aux %>% filter(x >= cap_cut))), Dots = d)
  dfy = tibble(Type = "OS", Axes = "y", Dist = as.vector(dist(aux %>% filter(y >= cap_cut))), Dots = d)
  dfz = tibble(Type = "OS", Axes = "z", Dist = as.vector(dist(aux %>% filter(z >= cap_cut))), Dots = d)
  
  df = df %>% bind_rows(dfx) %>% bind_rows(dfy) %>% bind_rows(dfz)
  
  i = 2
  aux = ss[[i]]$m_fb2
  d = paste0("Dots = ", dim(aux)[1])
  dfx = tibble(Type = "FIBOS", Axes = "x", Dist = as.vector(dist(aux %>% filter(x >= cap_cut))), Dots = d)
  dfy = tibble(Type = "FIBOS", Axes = "y", Dist = as.vector(dist(aux %>% filter(y >= cap_cut))), Dots = d)
  dfz = tibble(Type = "FIBOS", Axes = "z", Dist = as.vector(dist(aux %>% filter(z >= cap_cut))), Dots = d)
  
  df = df %>% bind_rows(dfx) %>% bind_rows(dfy) %>% bind_rows(dfz)
  
  aux = ss[[i]]$m_os
  dfx = tibble(Type = "OS", Axes = "x", Dist = as.vector(dist(aux %>% filter(x >= cap_cut))), Dots = d)
  dfy = tibble(Type = "OS", Axes = "y", Dist = as.vector(dist(aux %>% filter(y >= cap_cut))), Dots = d)
  dfz = tibble(Type = "OS", Axes = "z", Dist = as.vector(dist(aux %>% filter(z >= cap_cut))), Dots = d)
  
  df = df %>% bind_rows(dfx) %>% bind_rows(dfy) %>% bind_rows(dfz)
  
  p6 = df %>% ggplot(aes(x = Dist, color = Axes)) + 
    geom_density(key_glyph = draw_key_path) + 
    facet_grid(Type ~ Dots) + 
    ggtitle("SPHERICAL CAP DISTANCE DISTRIBUTION ACCORDING TO AXES",
            subtitle = paste0("Cap cut = ", cap_cut)) +
    xlab(expression("Distance (" * ring(A) * ")")) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = rel(1.1)),
          plot.subtitle = element_text(hjust=0.5, size = rel(1)),
          axis.title.y = element_blank(),
          strip.text = element_text(size=12),
          legend.position = "bottom")
  
  print(p6)
  
}

# DATA INPUT FOR VALIDATION OF R AND PYTHON OS VERSIONS
if(0){
  folder = "data" 
  filename = "osp_values_all_languages-v2"
  file_csv = path(folder, filename , ext = "csv")
  m = read_csv(file_csv)
  m = m %>% mutate(across(where(is.character) & !ID_PDB, 
                          ~ parse_double(.x, locale = locale(decimal_mark = ","))
                         )
                  )
}

# MAKE PLOT OF SUPPLEMENTARY FIGURE S4 OF THE ARTCILE
if(0){
  
  df = m %>% pivot_longer(-1, names_to = "Type", values_to = "OSP") %>% 
    separate(Type, c("Fun", "Code", "Type")) %>% 
    pivot_wider(names_from = "Code", values_from = "OSP") %>% 
    rename(OSP_FORTRAN = FORTRAN) %>% 
    rename(OSP_PYTHON = PYTHON) %>% 
    rename(OSP_R = R)

  tx = text_grob(paste("FIBOS OSP IN FORTRAN, PYTHON AND R"), size = 14, face = "bold")
  pt = as_ggplot(tx) + 
    theme(plot.margin = margin(0,0,0,0, "cm"))
  
  p7 = df %>% filter(Type == "FIBOS") %>% 
    ggscatter(x = "OSP_FORTRAN", y = "OSP_PYTHON", add = "reg.line", shape = 21, 
              color = "black", fill = "#FAA521", size = 4) +
    stat_cor(aes(label = after_stat(r.label)), label.x = 0.36, label.y = 0.45, 
             digits = 3, size = 5) #+

  p8 = df %>% filter(Type == "FIBOS") %>% 
    ggscatter(x = "OSP_FORTRAN", y = "OSP_R", add = "reg.line", shape = 21, 
              color = "black", fill = "#FAA521", size = 4) +
    stat_cor(aes(label = after_stat(r.label)), label.x = 0.36, label.y = 0.45, 
             digits = 3, size = 5) #+

  p9 = df %>% filter(Type == "FIBOS") %>% 
    ggscatter(x = "OSP_PYTHON", y = "OSP_R", add = "reg.line", shape = 21, 
              color = "black", fill = "#FAA521", size = 4) +
    stat_cor(aes(label = after_stat(r.label)), label.x = 0.36, label.y = 0.45, 
             digits = 3, size = 5) #+
  
  p789 = ggarrange(pt, ggarrange(p7, p8, p9, ncol=3), ncol = 1, heights = c(1,8))
  print(p789)

}

# GET AND EXPORT EXPERIMENTAL PDB INFORMATIONS TO CSV
if(0){
  tit = m$ID_PDB %>% map_dfr(get_pdb_info) ### API TO PDB 
  m.tab = bind_cols(tibble(ID_PDB = m$ID_PDB), tit) %>%
    mutate(TITLE = toupper(TITLE))
  
  folder = "data" 
  filename = "fleming-base-pdb-exp-title"
  file_csv = path(folder, filename , ext = "csv")
  
  write_csv(m.tab %>% select(ID_PDB, TITLE), file_csv)
  
  m.tab1 = m.tab %>% 
    mutate(R_FACTOR = if_else(is.na(R_FACTOR), R_FACTOR_W, R_FACTOR)) %>% 
    mutate(R_FACTOR = round(R_FACTOR,2)) %>% 
    select(!c(TITLE,R_FACTOR_W))
  
  m.pdb = m$ID_PDB %>% map(get_pdb_bio3d)
  m.sse = m.pdb %>% map_dfr(set_pdb_sse)
  m.sse = bind_cols(tibble(ID_PDB = m$ID_PDB), m.sse)
  
  m.all = m.tab1 %>% full_join(m.sse, by = "ID_PDB") 
  m.all = m.all %>% full_join(m, by = "ID_PDB")
  m.all = m.all %>% mutate(across(starts_with("OSP") , ~ round(.x, 3)))
  
  folder = "data" 
  filename = "fleming-base-pdb-exp"
  file_csv = path(folder, filename , ext = "csv")
  
  write_csv(m.all, file_csv)
}

  

  





