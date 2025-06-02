#!/usr/bin/env python3
#Pat Fleming, 2025
"""
Compares two output files from get_torsion.py

python3 chi_diffs.py PDB_torsions.dat AF_torsions.dat
where PDB_torsions.dat is,
id resn resi chi1 chi2 altchi2 chi3 chi4 chi5
7rsa_renum LYS 27 162  161  -168 82 
7rsa_renum GLU 28 -170  -178  33  
7rsa_renum THR 29 68      
7rsa_renum ALA 30       


Output in chi*_diffs.dat
Example: 
head chi1_diffs.dat
   27 133.00
   28   5.00
   29   1.00
   30   0.00

"""
import sys
import math
import csv
import numpy as np

def write_rotamer(filename, resid, chix, chiy):
    with open(filename, 'w', newline='', encoding='utf-8') as fo:
        writer = csv.writer(fo)
        # Opcional: escrever cabeçalho no CSV
        # writer.writerow(['resid', 'diff_ang'])
        
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
                
            writer.writerow([int(resid[i]), diff_ang])

#exp_data1 = open(sys.argv[1]).readlines()
#exp_data2 = open(sys.argv[2]).readlines()
exp_data1 = open("data_exp_tor_cls/1AWD.dat").readlines()
exp_data2 = open("data_csm_tor/1AWD.dat").readlines()
numlines_exp1 = len(exp_data1)
numlines_exp2 = len(exp_data2)
if numlines_exp1 != numlines_exp2:
    print('Number residues are not equal')

## Put experimental data in lists
resid = []
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

for i in range(1,numlines_exp1):
    line1 = exp_data1[i].strip("\n").split(" ")
    resid.append(line1[2])
    line2 = exp_data2[i].strip("\n").split(" ")
    
    '''
    if len(line1) >= 4:
        chi11.append(line1[3])
    else:
        chi11.append(0)
    if len(line2) >= 4:
        chi12.append(line2[3])
    else:
        chi12.append(0)

    if len(line1) >= 5:
        chi21.append(line1[4])
    else:
        chi21.append(0)
    if len(line2) >= 5:
        chi22.append(line2[4])
    else:
        chi22.append(0)

    if len(line1) >= 6:
        chi31.append(line1[5])
    else:
        chi31.append(0)
    if len(line2) >= 6:
        chi32.append(line2[5])
    else:
        chi32.append(0)

    if len(line1) >= 7:
        chi41.append(line1[6])
    else:
        chi41.append(0)
    if len(line2) >= 7:
        chi42.append(line2[6])
    else:
        chi42.append(0)

    if len(line1) >= 8:
        chi51.append(line1[7])
    else:
        chi51.append(0)
    if len(line2) >= 8:
        chi52.append(line2[7])
    else:
        chi52.append(0)
    '''

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


write_rotamer('chi1a_diffs.dat', resid, chi1a2, chi1a1)
'''
fo = open('chi1a_diffs.dat', 'w')
for i in range(0,len(chi1a1)):
    x = chi1a2[i]
    y = chi1a1[i]
    if (x != '') and (y != ''):
        diff = math.cos(np.deg2rad(float(x)) - np.deg2rad(float(y)))
        diff_ang = np.rad2deg(math.acos(diff))
    else:
        diff_ang = ''
    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()

with open('chi1a_diffs.dat', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[x] for x in chi1a1])

fo = open('chi1_diffs.dat', 'w')
for i in range(0,len(chi1a1)):
    diff = math.cos(np.deg2rad(float(chi1a2[i])) - np.deg2rad(float(chi1a1[i])))
    diff_ang = np.rad2deg(math.acos(diff))

    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()

fo = open('chi2_diffs.dat', 'w')
for i in range(0,len(chi21)):
    diff = math.cos(np.deg2rad(float(chi22[i])) - np.deg2rad(float(chi21[i])))
    diff_ang = np.rad2deg(math.acos(diff))

    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()

fo = open('chi3_diffs.dat', 'w')
for i in range(0,len(chi31)):
    diff = math.cos(np.deg2rad(float(chi32[i])) - np.deg2rad(float(chi31[i])))
    diff_ang = np.rad2deg(math.acos(diff))

    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()

fo = open('chi4_diffs.dat', 'w')
for i in range(0,len(chi41)):
    diff = math.cos(np.deg2rad(float(chi42[i])) - np.deg2rad(float(chi41[i])))
    diff_ang = np.rad2deg(math.acos(diff))

    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()

fo = open('chi5_diffs.dat', 'w')
for i in range(0,len(chi51)):
    diff = math.cos(np.deg2rad(float(chi52[i])) - np.deg2rad(float(chi51[i])))
    diff_ang = np.rad2deg(math.acos(diff))

    fo.write('%5d %6.2f\n' % (int(resid[i]),diff_ang))
fo.close()
'''


