#!/usr/bin/env python
# Filename: CartHess2FC.py
"""
This program is to get the force constant from Hessian Matrix
"""
from mcpb.gene_final_frcmod_file import get_bond_fc_with_sem, get_ang_fc_with_sem
from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtmol.getlist import get_blist, get_alist
from pymsmtmol.gauio import get_crds_from_fchk, get_matrix_from_fchk
from pymsmtmol.gmsio import get_crds_from_gms, get_matrix_from_gms
from optparse import OptionParser

parser = OptionParser("usage: -i pdbfile -f Hess_file -v software")
parser.set_defaults(softversion='g03')
parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-f", dest="hessfile", type='string',
                  help="Hessian file name")
parser.add_option("-v", dest="softversion", type='string',
                  help="Software version [Default is g03 (means Gaussian03), \n"
                       "           other option are, g09 (means Gaussian09), \n"
                       "                       and gms (means GAMESS-US)]")
(options, args) = parser.parse_args()

mol, atids, resids = get_atominfo_fpdb(options.inputfile)
blist = get_blist(mol, atids)
alist = get_alist(mol, blist)

natids = {}
for i in range(0, len(atids)):
  natids[atids[i]] = i + 1

#crds after optimization
if options.softversion == 'g03':
    crds = get_crds_from_fchk(options.hessfile, 'Current cartesian coordinates',
                            'Int Atom Types')
elif options.softversion == 'g09':
    crds = get_crds_from_fchk(options.hessfile, 'Current cartesian coordinates',
                            'Force Field')
elif options.softversion == 'gms':
    crds = get_crds_from_gms(options.hessfile)

#Whole Hessian Matrix
if options.softversion in ['g03', 'g09']:
    fcmatrix = get_matrix_from_fchk(options.hessfile, 'Cartesian Force Constants',
                                'Dipole Moment', 3*len(atids))
elif options.softversion == 'gms':
    fcmatrix = get_matrix_from_gms(options.hessfile, 3*len(atids))

for i in blist:
    at1 = i[0]
    at2 = i[1]
    nat1 = natids[at1]
    nat2 = natids[at2]
    dis, fcfinal = get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2)
    dis = round(dis, 4)
    fcfinal = round(fcfinal, 1)
    resid1 = mol.atoms[at1].resid
    resname1 = mol.atoms[at1].resname
    atname1 = mol.atoms[at1].atname
    resid2 = mol.atoms[at2].resid
    resname2 = mol.atoms[at2].resname
    atname2 = mol.atoms[at2].atname
    rep1 = resname1 + str(resid1) + '-' + atname1
    rep2 = resname2 + str(resid2) + '-' + atname2
    print "%4s %15s %15s %6.1f %9.4f" %("BOND", rep1, rep2, fcfinal, dis)

for i in alist:
    at1 = i[0]
    at2 = i[1]
    at3 = i[2]
    nat1 = natids[at1]
    nat2 = natids[at2]
    nat3 = natids[at3]
    angval, fcfinal = get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3)
    angval = round(angval, 2)
    fcfinal = round(fcfinal, 2)
    resid1 = mol.atoms[at1].resid
    resname1 = mol.atoms[at1].resname
    atname1 = mol.atoms[at1].atname
    resid2 = mol.atoms[at2].resid
    resname2 = mol.atoms[at2].resname
    atname2 = mol.atoms[at2].atname
    resid3 = mol.atoms[at3].resid
    resname3 = mol.atoms[at3].resname
    atname3 = mol.atoms[at3].atname
    rep1 = resname1 + str(resid1) + '-' + atname1
    rep2 = resname2 + str(resid2) + '-' + atname2
    rep3 = resname3 + str(resid3) + '-' + atname3
    print "%4s %15s %15s %15s %7.2f %11.2f" %("ANGL", rep1, rep2, rep3, fcfinal, angval)

quit()
