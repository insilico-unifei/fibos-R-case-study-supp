#!/usr/bin/env python3
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
# READS PDB FILE, EXTRACTS FIRST MODEL AS ATOM ARRAY, AND ANNOTATES SECONDARY STRUCTURE 
# VIA P-SEA. PRINTS LIST OF SSE ANNOTATIONS IN JSON FORMAT.
# '''

import sys
import json
from biotite.structure.io.pdb import PDBFile
import biotite.structure as struc
from biotite.structure import AtomArrayStack

def main():
    # Verifica se foi fornecido caminho do PDB
    if len(sys.argv) < 2:
        print(json.dumps([]))
        sys.exit(0)
    pdb_path = sys.argv[1]

    # Leitura do PDB
    pdb_file  = PDBFile.read(pdb_path)
    structure = pdb_file.get_structure()

    # Se for AtomArrayStack, extrai primeiro modelo
    if isinstance(structure, AtomArrayStack):
        atom_array = structure[0]
    else:
        atom_array = structure

    # Anota SSE usando P-SEA
    sse = struc.annotate_sse(atom_array)

    # Imprime lista P-SEA em JSON
    print(json.dumps(list(sse)))

if __name__ == "__main__":
    main()
