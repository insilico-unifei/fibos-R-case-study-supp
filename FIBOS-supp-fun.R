###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# R SCRIPT OF AUXILIARY FUNCTIONS USED BY FIBOS-supp.main.R
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
# library(devtools)
# install_github("https://github.com/insilico-unifei/FIBOS-R.git") 


# R VERSION FOR THE DISTRIBUTION OF DOTS ON THE SPHERE, ORIGINALLY CODED IN FORTRAN
GENUN = function(N){

  #browser()
  U = matrix(rep(0,3*N), nrow = 3, ncol = N)
  AR = rep(0,N)
  area = 0
  
  res = list()
  
  PI = pi
  
  onehalf = 1/2
  two = 2
  
  twoPI = two*PI
  
  #area = null 
  
  NEQUAT = floor(sqrt(N * PI))
  
  NVERT = floor(onehalf * NEQUAT)
  
  if (NVERT < 2) NVERT = 2
  
  NU = 0
  
  dNVERT = NVERT
  
  dtet = PI/dNVERT
  dtet2= dtet*onehalf
  sdtet =two* sin(dtet2)
  cdtet = cos(dtet2)
  
  for (I in 0:NVERT){
    #DO 100 I = 0,NVERT
    
    FI = dtet*I
    
    Z = cos(FI)
    
    XY = sin(FI)
    
    if(I == 0 | I == NVERT)
      aat = 1.0 - cdtet
    else{
      #browser()
      aat = XY*sdtet
      aat2 = cdtet
      #print(paste("aat= ", aat," aat2= ",aat2))
    }
    NHOR = floor(NEQUAT * XY)
    
    res$NHOR = c(res$NHOR, NHOR)
    
    if (NHOR < 1) NHOR = 1
    
    dNHOR = NHOR
    
    aaf = twoPI/dNHOR
    
    for (J in 0:(NHOR-1)){ 
    #DO 50 J = 0,NHOR-1
    
      FJ = aaf*J
      
      X = cos(FJ) * XY
      
      Y = sin(FJ) * XY
      
      if (NU >= N) break
      
      NU = NU + 1
      
      U[1,NU] = X
      
      U[2,NU] = Y
      
      U[3,NU] = Z
      
      AR[NU] = aaf * aat
      
      area = area + AR[NU]
      
      #50         CONTINUE
    }
    #100     CONTINUE
    if (NU >= N) break
  }
  #150     CONTINUE
  
  res = list(NHOR = res$NHOR, N = NU, U = U, AR = AR, area = area)
  
  return(res)

}

# R VERSION FOR THE DISTRIBUTION OF DOTS ON THE SPHERE, USING FIBONACCI SPIRALS
fibonaccisphere3 = function(n=1000, eps = 0.5, r=1, out.xyz=TRUE, out.sph=FALSE) {
  
  if (n<1 | round(n)!=n) stop("n must be a positive integer")
  if (!out.xyz & !out.sph) stop("either out.xyz and/or out.sph must be TRUE")
  
  goldenratio = (1+sqrt(5))/2

  i = seq(0, n-1)
  theta = 2*pi*i/goldenratio
  phi = acos(1 - 2*(i+eps)/(n-1+2*eps))

  if (out.xyz) {
    x = r * cos(theta) * sin(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(phi)
    if (out.sph) {
      out = cbind(x=x,y=y,z=z,theta=theta,phi=phi)
    } else {
      out = cbind(x=x,y=y,z=z)
    }
  } else {
    out = cbind(theta=theta,phi=phi)
  }
  return(out)
}

# CALCULATE VORONOI CELL AREA
Vor_Cell_Area = function(v, k = 180/pi, a = 1, f = 0){
  

  v.sph = cart2sph(t(v$cell))
  v.area = areaPolygon(v.sph[,1:2]*k, a = a, f = f)
  
  return(v.area)
  
}

# GET SOME EXPERIMENTAL PDB PARAMETERS
get_pdb_info = function(v, k = 2){
  
  aux = get_info(v)
  res = tibble(TITLE = aux$struct$title,
               RESIDUES = aux$rcsb_entry_info$deposited_polymer_monomer_count,
               RESOLUTION = aux$rcsb_entry_info$resolution_combined,
               R_FACTOR = aux$refine[1,"ls_rfactor_obs"],
               R_FACTOR_W = aux$refine[1,"ls_rfactor_rwork"]
               )
  return(res)
  
}

# GET PDB
get_pdb_bio3d = function(v){
  
  pdb = bio3d::read.pdb(v)
  return(pdb)
}

# GET SSE INFORMATION
set_pdb_sse = function(pdb, a = c("G","H","I"), b = c("B","E"), k = 0){
  
  sse = dssp(pdb)
  types = table(sse$sse)
  all = sum(pdb$calpha)
  #browser()
  f = names(types) %in% a
  alpha = round(sum(types[f])/all*100,k)
  f = names(types) %in% b
  beta = round(sum(types[f])/all*100,k)
  res = tibble(HELIX = alpha,
               STRAND = beta)
  return(res)
  
}
