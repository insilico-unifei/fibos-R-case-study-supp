###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# PYTHON SCRIPT ACCESSORY FOR: a case study using FIBOS, comparing packing densities 
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

# '''
# COMPUTES RELATIVE SASA FOR EACH RESIDUE USING SHRAKE RUPLEY.
# NORMALIZES BY MAX AMINO ACID SASA FROM TIEN ET AL. 2013.
# WRITES RESNUM, RESNAME, AND FRACTION TO REL_SASA.DAT.
# '''
# Maximum Amino Acid SASA
# Tien, M.Z. et al., Maximum Allowed Solvent Accessibilites of Residues in
#  Proteins, PLOSone 8:e80635, 2013

import sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
p = PDBParser(QUIET=1)
file = sys.argv[1]
struct = p.get_structure("prot", file)

AA_SASA = {
  'ALA': ( 129.0),
  'ARG': ( 274.0),
  'ASN': ( 195.0),
  'ASP': ( 193.0),
  'CYS': ( 167.0),
  'GLU': ( 223.0),
  'GLN': ( 225.0),
  'GLY': ( 104.0),
  'HIS': ( 224.0),
  'HSD': ( 224.0),
  'HSE': ( 224.0),
  'ILE': ( 197.0),
  'LEU': ( 201.0),
  'LYS': ( 236.0),
  'MET': ( 224.0),
  'PHE': ( 240.0),
  'PRO': ( 159.0),
  'SER': ( 155.0),
  'THR': ( 172.0),
  'TRP': ( 285.0),
  'TYR': ( 263.0),
  'VAL': ( 174.0)
}


sr = ShrakeRupley()
sr.compute(struct, level="R")

ffo = open('REL_SASA.dat', 'w')
for model in struct:
    for chain in model:
        for residue in chain:
            resname = residue.get_resname()
            if resname in AA_SASA:
              max_sasa = AA_SASA[resname]
              frac_sasa = residue.sasa/max_sasa
              #print(residue.get_resname(), residue.get_id()[1], frac_sasa) 
              ffo.write('%6d %3s %4.2f\n' % (residue.get_id()[1], resname, frac_sasa))

ffo.close()
