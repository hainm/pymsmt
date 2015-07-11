#!/usr/bin/env python
# Filename: CartHess2FC.py
"""
This program is to get the force constant from Hessian Matrix
"""
from __future__ import absolute_import
from mcpb.gene_final_frcmod_file import (get_bond_fc_with_sem, get_ang_fc_with_sem,
     get_dih_fc_with_sem, get_imp_fc_with_sem)
from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtmol.getlist import get_blist, get_all_list
from pymsmtmol.gauio import get_crds_from_fchk, get_matrix_from_fchk
from pymsmtmol.gmsio import get_crds_from_gms, get_matrix_from_gms
from optparse import OptionParser

parser = OptionParser("usage: -i pdbfile -f Hess_file [--scalef scale_factor] "
                      "[-v software]")
parser.set_defaults(softversion='g03', scalef=1.000, bondavg=0, angavg=0)
parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-f", dest="hessfile", type='string',
                  help="Hessian file name")
parser.add_option("--scalef", dest="scalef", type='float',
                  help="Scale factor")
parser.add_option("--bavg", dest="bondavg", type='int',
                  help="Bond force constant average? No (0, default) or Yes (1).")
parser.add_option("--aavg", dest="angavg", type='int',
                  help="Angle force constant average? No (0, default) or Yes (1).")
parser.add_option("-v", dest="softversion", type='string',
                  help="Software version [Default is g03 (means Gaussian03), \n"
                       "           other option are, g09 (means Gaussian09), \n"
                       "                       and gms (means GAMESS-US)]")
(options, args) = parser.parse_args()

mol, atids, resids = get_atominfo_fpdb(options.inputfile)
blist = get_blist(mol, atids)
all_list = get_all_list(mol, blist, atids, 8.0)

natids = {}
for i in range(0, len(atids)):
  natids[atids[i]] = i + 1

#crds after optimization
if options.softversion in ['g03', 'g09']:
    crds = get_crds_from_fchk(options.hessfile, len(atids))
elif options.softversion == 'gms':
    crds = get_crds_from_gms(options.hessfile)

#Whole Hessian Matrix
if options.softversion in ['g03', 'g09']:
    fcmatrix = get_matrix_from_fchk(options.hessfile, 3*len(atids))
elif options.softversion == 'gms':
    fcmatrix = get_matrix_from_gms(options.hessfile, 3*len(atids))

for i in all_list.bondlist:
    at1 = i[0]
    at2 = i[1]
    nat1 = natids[at1]
    nat2 = natids[at2]

    if options.bondavg == 1:
        dis, fcfinal, stdv = get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2, options.scalef, options.bondavg)
    elif options.bondavg == 0:
        dis, fcfinal = get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2, options.scalef, options.bondavg)

    print '### Bond force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + ' and ' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + ' :',

    if options.bondavg == 1:
        print str(round(fcfinal, 1)) + ' with StandardDev ' + str(round(stdv, 1)) + \
              ' with bond distance : ' + str(round(dis, 4))
    elif options.bondavg == 0:
        print str(round(fcfinal, 1)) + ' with bond distance : ' + str(round(dis, 4))

for i in all_list.anglist:
    at1 = i[0]
    at2 = i[1]
    at3 = i[2]
    nat1 = natids[at1]
    nat2 = natids[at2]
    nat3 = natids[at3]

    if options.angavg == 1:
        angval, fcfinal, stdv = get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, options.scalef, options.angavg)
    elif options.angavg == 0:
        angval, fcfinal = get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, options.scalef, options.angavg)

    print '### Angle force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + ', ' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + ' and ' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + ' :',

    if options.angavg == 1:
        print str(round(fcfinal, 2)) + ' with StandardDev ' + str(round(stdv, 2)) + \
              ' with angle value : ' + str(round(angval, 2))
    elif options.angavg == 0:
        print str(round(fcfinal, 2)) + ' with angle value : ' + str(round(angval, 2))

for i in all_list.dihlist:
    at1 = i[0]
    at2 = i[1]
    at3 = i[2]
    at4 = i[3]

    n1 = 0
    n2 = 0
    for j in all_list.bondlist:
        if (at2, at3, 1) != j and (at3, at2, 1) != j:
            if at2 in j:
                n1 = n1 + 1
            elif at3 in j:
                n2 = n2 + 1

    nat1 = natids[at1]
    nat2 = natids[at2]
    nat3 = natids[at3]
    nat4 = natids[at4]
    dihval, fcfinal1, fcfinal = get_dih_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, nat4, n1, n2, options.scalef)
    print '### Dihedral force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + ', ' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + ', ' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + ' and ' + \
          mol.atoms[at4].resname + str(mol.atoms[at4].resid) + '@' + mol.atoms[at4].atname + ' : ' + \
          str(round(fcfinal, 2)) + ' with angle value ' + str(round(dihval, 2)) + \
          ' *or* ' + str(round(fcfinal1, 2)) + ' (as kpi in kpi(pi-pi0)^2), pi is the dihedral.'

for i in all_list.implist:
    at1 = i[0]
    at2 = i[1]
    at3 = i[2] #Central atom
    at4 = i[3]

    nat1 = natids[at1]
    nat2 = natids[at2]
    nat3 = natids[at3]
    nat4 = natids[at4]

    if mol.atoms[at1].element == mol.atoms[at2].element:
        fcfinal1, fcfinal = get_imp_fc_with_sem(crds, fcmatrix, nat3, nat1, nat2, nat4, options.scalef) #Treat the central atom first
        print '### Improper torsion force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + '-' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + '-' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + '(central atom)-' + \
          mol.atoms[at4].resname + str(mol.atoms[at4].resid) + '@' + mol.atoms[at4].atname + ' : ' + \
          str(round(fcfinal, 2)) 
        print '    *or* ' + str(round(fcfinal1, 2)) + ' (as kb in kb(b)^2), b is the distance of central atom to the plane of the other three atoms.'
    elif mol.atoms[at1].element == mol.atoms[at4].element:
        fcfinal1, fcfinal = get_imp_fc_with_sem(crds, fcmatrix, nat3, nat1, nat4, nat2, options.scalef) #Treat the central atom first
        print '### Improper torsion force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + '-' + \
          mol.atoms[at4].resname + str(mol.atoms[at4].resid) + '@' + mol.atoms[at4].atname + '-' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + '(central atom)-' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + ' : ' + \
          str(round(fcfinal, 2))
        print '    *or* ' + str(round(fcfinal1, 2)) + ' (as kb in kb(b)^2), b is the distance of central atom to the plane of the other three atoms.'
    elif mol.atoms[at2].element == mol.atoms[at4].element:
        fcfinal1, fcfinal = get_imp_fc_with_sem(crds, fcmatrix, nat3, nat2, nat4, nat1, options.scalef) #Treat the central atom first
        print '### Improper torsion force constant between ' + \
          mol.atoms[at4].resname + str(mol.atoms[at4].resid) + '@' + mol.atoms[at4].atname + '-' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + '-' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + '(central atom)-' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + ' : ' + \
          str(round(fcfinal, 2))
        print '    *or* ' + str(round(fcfinal1, 2)) + ' (as kb in kb(b)^2), b is the distance of central atom to the plane of the other three atoms.'
    else:
        fcfinal1, fcfinal = get_imp_fc_with_sem(crds, fcmatrix, nat3, nat1, nat2, nat4, options.scalef) #Treat the central atom first
        print '### Improper torsion force constant between ' + \
          mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + '-' + \
          mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + '-' + \
          mol.atoms[at3].resname + str(mol.atoms[at3].resid) + '@' + mol.atoms[at3].atname + '(central atom)-' + \
          mol.atoms[at4].resname + str(mol.atoms[at4].resid) + '@' + mol.atoms[at4].atname + ' : ' + \
          str(round(fcfinal, 2))
        print '    *or* ' + str(round(fcfinal1, 2)) + ' (as kb in kb(b)^2), b is the distance of central atom to the plane of the other three atoms.'
 
quit()
