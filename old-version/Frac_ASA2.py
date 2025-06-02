import sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBIO

def remove_lower_occupancy_alternates(structure):
    """
    For each residue in the structure, keep only the atom with the highest occupancy
    (for each atom name) and remove the others.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                best_atoms = {}
                # Iterate over atoms in the residue (using list() to safely modify the residue)
                for atom in list(residue):
                    atom_name = atom.get_name()
                    occupancy = atom.get_occupancy()
                    # If occupancy is None, assume 0.0
                    if occupancy is None:
                        occupancy = 0.0
                    # Record this atom if it's the first seen or if it has a higher occupancy
                    if atom_name not in best_atoms or occupancy > best_atoms[atom_name].get_occupancy():
                        best_atoms[atom_name] = atom
                # Remove atoms that are not the best (i.e., have lower occupancy)
                for atom in list(residue):
                    if best_atoms[atom.get_name()] != atom:
                        residue.detach_child(atom.get_id())


p = PDBParser(QUIET=1)
file = sys.argv[1]
#print(f"Calculating ASA for {file}")
struct = p.get_structure("prot", file)

remove_lower_occupancy_alternates(struct)

io = PDBIO()
io.set_structure(struct)
io.save("test.pdb")

# Maximum Amino Acid SASA
# Tien, M.Z. et al., Maximum Allowed Solvent Accessibilites of Residues in
#  Proteins, PLOSone 8:e80635, 2013
#
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
