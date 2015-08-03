#!/usr/bin/env python
# Filename: AFQMMM_NMR.py
"""
This is the AFQMMM_NMR.py program written by Pengfe Li from Merz Research
Group from Michigan State University. All Rights Reserved. The program is
not gurantee to work due to bug may exist. Suggestions and bug reports are
welcome to send to Pengfei Li (Email address: ldsoar1990@gmail.com).

Please cite the following paper if you use the software to perform the
modeling:

The Automated Fragmentation QM/MM approach:
** X. He, B. Wang and K. M. Merz. JPCB, 2009, 113, 10380-10388.
"""
from pymsmtapi.AmberParm import read_amber_prm
from pymsmtmol.cal import calc_bond
from pymsmtmol.element import Atnum, CoRadiiDict, ResChgDict, resnamel
from pymsmtexp import *

from chemistry.amber.readparm import AmberParm
from chemistry.amber.mask import AmberMask
from optparse import OptionParser
import os

#==============================================================================
# Define related functions
#==============================================================================
def get_spin(fname, totchg):
    fpr = open(fname, 'r')
    elmt = 0
    for line in fpr:
      if line[0:2] in [' C', ' O', ' S']:
        elmt = elmt + 2
      elif line[0:2] in [' H', ' N']:
        elmt = elmt + 1
    fpr.close()
    elmt = elmt - int(round(totchg, 0))
    if elmt%2 == 0:
      spin = 1
    elif elmt%2 == 1:
      spin = 2
    return spin

def afqmmm_nmr(cresids, mol, atids, resids):
    #==========================================================================
    # For each core residue, do the analysis and print out work
    #==========================================================================
    bdld = {'CH': 1.090, 'NH': 1.010}

    #determine the amino acid residues
    aaresids = []
    for j in resids:
      hasatoms = []
      for k in mol.residues[j].resconter:
        hasatoms.append(mol.atoms[k].atname)
      if set(['C', 'O', 'CA', 'N']) <= set(hasatoms):
        aaresids.append(j)

    for i in cresids:
      if i == 0:
        raise pymsmtError('Residue number could not be %s.' % i)
      if i not in aaresids:
        raise pymsmtError('The core residue %s is not an amino acid.' % i)

      fp = open('number' + str(i) + '.com', 'w') #Gaussian input file
      fpf = open('number' + str(i) + '.fpf', 'w') #finger print file
      infof = open('number' + str(i) + '.info', 'w') #information file
      xyzf = open('number' + str(i) + '.xyz', 'w') #xyz file

      print >> xyzf, 'NNNNN'
      print >> xyzf, 'Title'

      print >> infof, 'This molecule has ' + str(len(resids)) + ' residues.'
      print >> infof, 'Core residue = ', i

      print >> fp, '%%chk=number%d.chk' % i
      print >> fp, '%mem=22GB'
      print >> fp, '%nproc=8'
      print >> fp, '#b3lyp/6-31G* charge nmr'
      print >> fp, ' '
      print >> fp, 'Have a nice day'

      bresids = [] #buffer region
      ccaps = [] #C terminal capped residue
      ncaps = [] #N terminal capped residue

      catoms = []
      im1 = i - 1

      if (im1 not in aaresids) and (im1 != 0):
        raise pymsmtError('Could not find the C, O backbone atoms in '
                          'the residue %s, which is connect to the core '
                          'residue %s' % (im1, i))
      ##=======================================================================
      ##  get the core atoms
      ##=======================================================================
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname not in ['C', 'O']:
          catoms.append(j)
      if i != 1:
        for j in mol.residues[im1].resconter:
          atname = mol.atoms[j].atname
          if atname in ['C', 'O']:
            catoms.append(j)
      ##=======================================================================
      ##  get the buffer residues
      ##=======================================================================
      #1. distance creteria
      for j in catoms:
        hnum1 = 0
        elmt1 = mol.atoms[j].element
        crd1 = mol.atoms[j].crd
        if elmt1 == 'H':
          hnum1 = 1
        for k in atids:
          if k != j:
            hnum2 = 0
            elmt2 = mol.atoms[k].element
            crd2 = mol.atoms[k].crd
            dis = calc_bond(crd1, crd2)
            at2 = mol.atoms[k].atname
            resid2 = mol.atoms[k].resid
            if elmt2 == 'H':
              hnum2 = 1
            hnums = hnum1 + hnum2
            #first criterion
            if hnums == 2 and dis <= 3.0:
              if mol.atoms[k].atname not in ['C', 'O']:
                if mol.atoms[k].resid not in bresids:
                  bresids.append(mol.atoms[k].resid)
              else:
                if (mol.atoms[k].resid+1 not in bresids):
                  bresids.append(mol.atoms[k].resid+1)
            #second criterion
            elif hnums == 1 and dis <= 4.0:
              if mol.atoms[k].atname not in ['C', 'O']:
                if mol.atoms[k].resid not in bresids:
                  bresids.append(mol.atoms[k].resid)
              else:
                if mol.atoms[k].resid+1 not in bresids:
                  bresids.append(mol.atoms[k].resid+1)
            #third criterion
            elif (mol.atoms[k].resname in ['HIE', 'HIP', 'HID', 'HIS', 'PHE', \
                 'TYR', 'TRP']) and (at2 not in ['O', 'C', 'N', 'H', 'CA', \
                 'HA', 'CB', 'HB2', 'HB3', 'OH', 'HO']) and (dis <= 5.0):
              if mol.atoms[k].atname not in ['C', 'O']:
                if mol.atoms[k].resid not in bresids:
                  bresids.append(mol.atoms[k].resid)
              else:
                if mol.atoms[k].resid+1 not in bresids:
                  bresids.append(mol.atoms[k].resid+1)
      print >> infof, 'After the distance criterion: Buffer residues =       '\
               , sorted(bresids)
      #2. core creteria
      for j in range(i-2, i+3):
        if (j > 0) and (j <= max(resids)):
          if (j in aaresids) and (j not in bresids):
            bresids.append(j)
      print >> infof, 'After the core criterion: Buffer residues =           '\
               , sorted(bresids)
      #3. disulfur bond creteria
      for j in bresids:
        if j<= max(resids):
          if mol.residues[j].resname == 'CYX':
            for k in mol.residues[j].resconter:
              if mol.atoms[k].atname == 'SG':
                for l in atids:
                  if (l != k) and (mol.atoms[l].atname == 'SG') and \
                  (mol.atoms[l].resname == 'CYX'):
                    dis = calc_bond(mol.atoms[l].crd, mol.atoms[k].crd)
                    if (dis <= 2.5) and (mol.atoms[l].resid not in bresids):
                      bresids.append(mol.atoms[l].resid)
      print >> infof, 'After the disulfide bond criterion: Buffer residues = '\
               , sorted(bresids)
      ####=====================================================================
      ####  For special situation
      ####=====================================================================
      #for N terminal
      if (resids[0] in bresids) and (resids[1] not in bresids):
        if (resids[1] in aaresids):
          bresids.append(resids[1])

      #for C terminal
      if (max(aaresids)+1 in bresids):
        hasct = 1
        bresids.remove(max(aaresids)+1)
        if (max(aaresids) not in bresids):
          bresdis.append(max(aaresids))
      else:
        hasct = 0

      #get C-terminal and N-terminal caps
      for j in bresids:
        if (j-1 not in bresids) and (j-1 > 0) and (j-1 in aaresids):
          ncaps.append(j-1)#ncaps

      if hasct == 1:
        for j in bresids:
          if (j+1 not in bresids) and (j != max(aaresids)):
            ccaps.append(j)
      else:
        for j in bresids:
          if (j+1 not in bresids):
            ccaps.append(j)#ccaps

      bresids = list(set(bresids) - set(ccaps))
      bresids.sort()

      print >> infof, 'After the special case criterion: Buffer residues =   '\
               , sorted(bresids)
      print >> infof, 'After the special case criterion: N-terimal caps =    '\
               , sorted(ncaps)
      print >> infof, 'After the special case criterion: C-terimal caps =    '\
               , sorted(ccaps)
      ####=====================================================================
      ####  get the charge and spin of the two structures
      ####=====================================================================
      totchg = 0.0
      chgres = [i]
      chgres = list(set(chgres + bresids))

      for j in chgres:
        chglist = [mol.atoms[k].charge for k in mol.residues[j].resconter]
        chgj = sum(chglist)
        totchg = totchg + float(chgj)

      print >> fp, ' '
      print >> fp, '%-5d %5s' %(int(round(totchg, 0)), 'SSSSS')
      ####=====================================================================
      ####  Print the core region
      ####=====================================================================
      #print >> fp, "Core region residue: ", i

      for j in catoms:
        element = mol.atoms[j].element
        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        print >> fpf, str(mol.atoms[j].resid) + '-' + \
                 mol.atoms[j].resname + '-' + mol.atoms[j].atname
      ####=====================================================================
      ####  Print the N-terminal cap
      ####=====================================================================
      hcapnum = 0
      batoms = [] #buffer atoms

      for j in ncaps:
        #print >> infof, "N-terminal cap residues: ", j
        for k in mol.residues[j].resconter:
          if mol.atoms[k].atname == 'C':
            ccrd = mol.atoms[k].crd
        for k in mol.residues[j].resconter:
          if mol.atoms[k].atname in ['O', 'C']:
            batoms.append(k)
            element = mol.atoms[k].element
            crdx = mol.atoms[k].crd[0]
            crdy = mol.atoms[k].crd[1]
            crdz = mol.atoms[k].crd[2]
            print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
            print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
          elif mol.atoms[k].atname == 'CA': #change to H
            hcapnum += 1
            element = 'H'
            bvec = calc_bond(ccrd, mol.atoms[k].crd)
            crdx = ccrd[0] + bdld['CH'] * (mol.atoms[k].crd[0] - ccrd[0])/bvec
            crdy = ccrd[1] + bdld['CH'] * (mol.atoms[k].crd[1] - ccrd[1])/bvec
            crdz = ccrd[2] + bdld['CH'] * (mol.atoms[k].crd[2] - ccrd[2])/bvec
            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)
            print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
            print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
      ####=====================================================================
      ####  print the C-terminal cap
      ####=====================================================================
      for j in ccaps:
        #print >> fp, "C-terminal cap residues: ", j
        for k in mol.residues[j].resconter:
          if mol.atoms[k].atname == 'CA':
            cacrd = mol.atoms[k].crd

        for k in mol.residues[j].resconter:
          if mol.atoms[k].atname != 'O':
            if mol.atoms[k].atname == 'C': #change to H
              hcapnum += 1
              element = 'H'
              bvec = calc_bond(cacrd, mol.atoms[k].crd)
              crdx = cacrd[0] + bdld['CH'] * (mol.atoms[k].crd[0] - cacrd[0])\
                     /bvec
              crdy = cacrd[1] + bdld['CH'] * (mol.atoms[k].crd[1] - cacrd[1])\
                     /bvec
              crdz = cacrd[2] + bdld['CH'] * (mol.atoms[k].crd[2] - cacrd[2])\
                     /bvec
              crdx = round(crdx, 3)
              crdy = round(crdy, 3)
              crdz = round(crdz, 3)
              print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                       crdz)
              print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                       crdz)
            else:
              batoms.append(k)
              element = mol.atoms[k].element
              crdx = mol.atoms[k].crd[0]
              crdy = mol.atoms[k].crd[1]
              crdz = mol.atoms[k].crd[2]
              print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                       crdz)
              print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                       crdz)

      ##for the speicial case
      if hasct == 1:
        j = max(aaresids)
        for k in mol.residues[j].resconter:
          if mol.atoms[i].atname in ['C', 'O', 'OXT']:
            batoms.append(k)
            element = mol.atoms[k].element
            crdx = mol.atoms[k].crd[0]
            crdy = mol.atoms[k].crd[1]
            crdz = mol.atoms[k].crd[2]
            print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
            print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                     crdz)
      ####=====================================================================
      ####  print the other buffer resids
      ####=====================================================================
      for j in bresids:
        #print >> fp, "Buffer residues: ", j
        for k in mol.residues[j].resconter:
          if k not in catoms:
            batoms.append(k)
            element = mol.atoms[k].element
            crdx = mol.atoms[k].crd[0]
            crdy = mol.atoms[k].crd[1]
            crdz = mol.atoms[k].crd[2]
            print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, \
                     crdz)
            print >> xyzf, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy,\
                     crdz)

      print >> fp, ' '
      ####=====================================================================
      ####  print the other parts, background charge
      ####=====================================================================
      #print len(atids), len(catoms), len(batoms)
      nums = 0
      for j in atids:
        if (j not in catoms) and (j not in batoms):
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          resname = mol.atoms[j].resname
          atname = mol.atoms[j].atname
          chg = mol.atoms[j].charge
          print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
          nums += 1
          #elif mol.atoms[j].resid == resids[-2]:
          #  crdx = mol.atoms[j].crd[0]
          #  crdy = mol.atoms[j].crd[1]
          #  crdz = mol.atoms[j].crd[2]
          #  resname = mol.atoms[j].resname
          #  atname = mol.atoms[j].atname
          #  chg = mol.atoms[j].charge
          #  print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
          #elif mol.residues[mol.atoms[j].resid].resname != 'LIG':
          #  crdx = mol.atoms[j].crd[0]
          #  crdy = mol.atoms[j].crd[1]
          #  crdz = mol.atoms[j].crd[2]
          #  resname = mol.atoms[j].resname
          #  atname = mol.atoms[j].atname
          #  chg = mol.atoms[j].charge
          #  print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
      #print nums
      if nums + len(catoms) + len(batoms) == len(atids):
        pass
      else:
        raise pymsmtError('The output atom numbers are not consistent with '
                          'former toplogy file.')

      print >> fp, ' '
      print >> fp, ' '
      print >> fp, ' '

      fp.close()
      fpf.close()
      infof.close()
      xyzf.close()

      #get spin
      spin = get_spin('number' + str(i) + '.com', totchg)
      os.system("sed -i '' 's/SSSSS/%-5d/g' %s" %(spin, 'number' + str(i) + \
                '.com'))

      #get the xyz number
      xyznum = \
        os.popen("awk 'END {print NR}' %s" %('number'+str(i)+'.xyz')).read()
      xyznum = int(xyznum) - 2
      os.system("sed -i '' 's/NNNNN/%-5d/g' %s" %(xyznum, \
        'number' + str(i) + '.xyz'))

      #check whether the atom number is the same
      cfatnum = \
        os.popen("awk 'END {print NR}' %s" %('number'+str(i)+'.com')).read()
      cfatnum = int(cfatnum) - 12 #com file atom number
      topatnum = len(atids) + hcapnum #top file atom number

      if cfatnum == topatnum:
        pass
      else:
        raise pymsmtError('The atnum are not consistent between the com file '
                          'and toplogy file.')

#==============================================================================
# Setting the variables
#==============================================================================
parser = OptionParser("usage: -p toplogy_file -r coordinate_file "
                      "-m amber_mask [-w (0 to 3)]")

parser.set_defaults(water=1)

parser.add_option("-p", dest="pfile", help="topology file name")
parser.add_option("-c", dest="cfile", help="coordinate file name")
parser.add_option("-m", dest="mask", help="amber mask of the core residue")
parser.add_option("-w", type='int', dest="water", help="whether delete "
                  "water and ion. 0: no, 1: delete water and ions, 2: "
                  "delete water, 3: delete ions. [default: 1]")
(options, args) = parser.parse_args()

#==============================================================================
# Get the molecule information
#==============================================================================
prmtop, mol, atids, resids = read_amber_prm(options.pfile, options.cfile)

if options.water == 0:
  pass
if options.water == 1:
  mol.delwaterion()
elif options.water == 2:
  mol.delwater()
elif options.water == 3:
  mol.delion()
elif options.water >= 4:
  raise pymsmtError('bad choice of the -w options')

atids = mol.atoms.keys()
resids = mol.residues.keys()
mask = AmberMask(prmtop, options.mask)

cresids = []

for i in mask.Selected():
  if mol.atoms[i+1].resid not in cresids:
    cresids.append(mol.atoms[i+1].resid)

if set(cresids) <= set(resids):
  pass
else:
  raise pymsmtError('bad choice of the -m options.')

afqmmm_nmr(cresids, mol, atids, resids)

quit()

