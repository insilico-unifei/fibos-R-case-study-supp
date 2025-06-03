# FIBOS Case Study

## Description

FIBOS is an R and Python package designed to estimate the packing density of biomolecules.

In this case study, we perform in R a comparative analysis of packing densities between 
experimentally determined structures (EXP) and those predicted by AlphaFold (AF).

Additionally, we include the codes used to generate the supplementary materials for the article.

## R Files

1.  **`FIBOS-study-case-main-v2.R`**: Main R script used in the case study (version 2) 

2.  **`FIBOS-study-case-fun-v2.R`**: R file containing auxiliary functions used by `FIBOS-study-case.main.R` (version 2)
3.  **`FIBOS-study-case-analysis.R`**: Main R script used for statistical and plotting analysis

4.  **`FIBOS-supp-main.R`**: Main R script used for producing supplementary material 

5.  **`FIBOS-supp-fun.R`**: R file containing auxiliary functions used by `FIBOS-supp.main.R` 

## Python scripts (accessory)

1.  **`get_raw_pdbs.py`**: download raw pdb ids

2.  **`get_torsions.py`**: obtain the side chain chi angles in separated files 

3.  **`compare_chis2.py`**: Compares two output files from `get_torsions.py`

4.  **`Frac_ASA.py`**: computes relative SASA for each residue using Shrake Rupley

5.  **`sse_psea.py`**: annotates secondary structure  via P-SEA

## Packages and Instalation Details

FIBOS R: <https://github.com/insilico-unifei/fibos-R.git> \
FIBOS Python: <https://github.com/insilico-unifei/fibos-py.git>

## Authors

-   Carlos Silveira ([carlos.silveira\@unifei.edu.br](mailto:carlos.silveira@unifei.edu.br))\
    Herson Soares ([d2020102075\@unifei.edu.br](mailto:d2020102075@unifei.edu.br))\
    Institute of Technological Sciences,\
    Federal University of Itajubá,\
    Campus Itabira, Brazil.

-   João Romanelli ([joaoromanelli\@unifei.edu.br](mailto:joaoromanelli@unifei.edu.br)) \
    Institute of Applied and Pure Sciences, \
    Federal University of Itajubá, \
    Campus Itabira, Brazil.

-   Patrick Fleming ([Pat.Fleming\@jhu.edu](mailto:Pat.Fleming@jhu.edu)) \
    Thomas C. Jenkins Department of Biophysics, \
    Johns Hopkins University, \
    Baltimore, MD, USA

## References

Fleming PJ, Richards FM. Protein packing: Dependence on protein size, secondary structure and amino acid composition. J Mol Biol 2000;299:487–98.

Pattabiraman N, Ward KB, Fleming PJ. Occluded molecular surface: Analysis of protein packing. J Mol Recognit 1995;**8**:334–44.
