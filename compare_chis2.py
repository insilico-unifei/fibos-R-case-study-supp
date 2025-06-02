#!/usr/bin/env python3
#Pat Fleming, 2025
#Carlos Silveira, 2025

# """
# COMPUTES CHI ANGLE DIFFERENCES BETWEEN TWO PDB TORSION FILES AND CALCULATES ANGLE DIFFS.
# WRITES RESULTS FOR CHI1–CHI5 AND ALT CHI ANGLES TO SEPARATE OUTPUT DAT FILES.
# 
# Compares two output files from get_torsion.py
# 
# python3 chi_diffs.py PDB_torsions.dat AF_torsions.dat
# where PDB_torsions.dat is,
# id resn resi chi1 chi2 altchi2 chi3 chi4 chi5
# 7rsa_renum LYS 27 162  161  -168 82 
# 7rsa_renum GLU 28 -170  -178  33  
# 7rsa_renum THR 29 68      
# 7rsa_renum ALA 30       
# 
# 
# Output in chi*_diffs.dat
# Example: 
# head chi1_diffs.dat
#    27 133.00
#    28   5.00
#    29   1.00
#    30   0.00
# 
# """


import sys
import math
import csv
import numpy as np


def write_rotamers(filename, resid, chix, chiy, pdbid='', resn=''):
    with open(filename, 'w', newline='', encoding='utf-8') as fo:
        writer = csv.writer(fo)
        
        for i in range(len(chiy)):
            #i = 82
            x = chix[i]
            y = chiy[i]
            if (x != '') and (y != ''):
                diff = math.cos(np.deg2rad(float(x)) - np.deg2rad(float(y)))
                diff_ang = np.rad2deg(math.acos(diff))
                diff_ang = round(diff_ang) 
            else:
                diff_ang = "NA"
            if (pdbid != '') and (resn != ''):
                writer.writerow([pdbid[i], resn[i], int(resid[i]), diff_ang])
            else:
                writer.writerow([int(resid[i]), diff_ang])


def main():
    exp_data1 = open(sys.argv[1]).readlines()
    exp_data2 = open(sys.argv[2]).readlines()
    #exp_data1 = open("data_exp_tor_cls/1AWD.dat").readlines()
    #exp_data2 = open("data_csm_tor/1AWD.dat").readlines()
    numlines_exp1 = len(exp_data1)
    numlines_exp2 = len(exp_data2)
    if numlines_exp1 != numlines_exp2:
        print('Number residues are not equal')

    pdbid = []
    resid = []
    resn = []
    chi11 = []
    chi1a1 = []
    chi21 = []
    chi2a1 = []
    chi31 = []
    chi41 = []
    chi51 = []

    chi12 = []
    chi1a2 = []
    chi22 = []
    chi2a2 = []
    chi32 = []
    chi42 = []
    chi52 = []

    for i in range(1, numlines_exp1):
        line1 = exp_data1[i].strip("\n").split(" ")
        pdbid.append(line1[0])
        resn.append(line1[1])
        resid.append(line1[2])
        line2 = exp_data2[i].strip("\n").split(" ")

        chi11.append(line1[3])
        chi1a1.append(line1[4])
        chi21.append(line1[5])
        chi2a1.append(line1[6])
        chi31.append(line1[7])
        chi41.append(line1[8])
        chi51.append(line1[9])

        chi12.append(line2[3])
        chi1a2.append(line2[4])
        chi22.append(line2[5])
        chi2a2.append(line2[6])
        chi32.append(line2[7])
        chi42.append(line2[8])
        chi52.append(line2[9])

    write_rotamers('chi1_diffs.dat', resid, chi12, chi11, pdbid, resn)
    write_rotamers('achi1_diffs.dat', resid, chi1a2, chi1a1)
    write_rotamers('chi2_diffs.dat', resid, chi22, chi21)
    write_rotamers('achi2_diffs.dat', resid, chi2a2, chi2a1)
    write_rotamers('chi3_diffs.dat', resid, chi32, chi31)
    write_rotamers('chi4_diffs.dat', resid, chi42, chi41)
    write_rotamers('chi5_diffs.dat', resid, chi52, chi51)


if __name__ == '__main__':
    main()
