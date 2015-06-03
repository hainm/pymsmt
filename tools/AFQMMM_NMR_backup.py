"""
This is the AFQMMM_NMR.py program written by Pengfe Li from Merz Research
Group from Michigan State University. The method is from X. He, B. Wang
and K. M. Merz. JPCB, 2009, 113, 10380-10388.
"""

from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtlib.lib import get_lib_dict
from pymsmtmol.cal import calc_bond
from pymsmtmol.element import Atnum
from optparse import OptionParser

#==============================================================================
# Setting the usuages and input variables
#==============================================================================

parser = OptionParser("usage: -p pdb_file -r residue_ids \
                      -l contain_ligand(0/1)")
parser.add_option("-p", dest="fname", help="pdb file name")
parser.add_option("-r", dest="resid", help="resid which was treated as core \
                  residue")
parser.add_option("-l", type='int', \
                  dest="lyn", help="whether to contain the ligand")           
(options, args) = parser.parse_args()

#==============================================================================
# Deal with the input variables
#==============================================================================

if options.fname0 is None:
  pass
else:
  fname0 = options.fname0
  mol0, atids0 = get_atominfo_fpdb(fname0)
  resids0 = []
  for i in atids0:
    resid = mol0.atoms[i].resid
    if resid not in resids0:
      resids0.append(resid) 

fname = options.fname

cresids0 = []
cresids = []

bdld = {'CH': 1.090, 'NH': 1.010}

fp0 = open(options.residf, 'r')
for line in fp0:
  line = line.strip('\n')
  line = line.strip(' ')
  try:
    line = int(line)
  except:
    pass
  cresids0.append(line)
fp0.close()

parmdict, chargedict = get_lib_dict('ff99SB')

if options.mol2f is None:
  pass
else:
  parmdict1, chargedict1 = get_lib_dict(options.mol2f)
  parmdict.update(parmdict1)
  chargedict.update(chargedict1)

mol, atids = get_atominfo_fpdb(fname)

resids = []
for i in atids:
  resid = mol.atoms[i].resid
  if resid not in resids:
    resids.append(resid) 

if options.fname0 is None:
  cresids = cresids0
else:
  for i in range(0, len(resids)):
    if resids0[i] in cresids0:
      cresids.append(resids[i])
      
#==============================================================================
# Calculate the charge and spin number of the system
#==============================================================================

totchg = 0.0
elmt = 0

for i in resids:
  if i == 1:
    chgi = chargedict['N' + mol.residues[i].resname]
  elif i == resids[-2]:
    chgi = chargedict['C' + mol.residues[i].resname]
  else:
    chgi = chargedict[mol.residues[i].resname]
  totchg = totchg + chgi
  
  for j in mol.residues[i].resconter:
    element = mol.atoms[j].element
    atnum = Atnum[element]
    elmt = elmt + atnum

elmt = elmt - int(round(totchg, 0))

if elmt%2 == 0:
  spin = 1
else:
  spin = 2

#--------------------
totchg2 = 0.0
elmt2 = 0

for i in resids:
  if i == 1:
    chgi = chargedict['N' + mol.residues[i].resname]
  elif i == resids[-2]:
    chgi = chargedict['C' + mol.residues[i].resname]
  elif mol.residues[i].resname != 'LIG':
    chgi = chargedict[mol.residues[i].resname]
  elif mol.residues[i].resname == 'LIG':
    chgi = 0.0

  totchg2 = totchg2 + chgi
  
  for j in mol.residues[i].resconter:
    if mol.residues[mol.atoms[j].resid].resname != 'LIG':
      element = mol.atoms[j].element
      atnum = Atnum[element]
      elmt2 = elmt2 + atnum

elmt2 = elmt2 - int(round(totchg2, 0))

if elmt2%2 == 0:
  spin2 = 1
else:
  spin2 = 2

#==============================================================================
# For each core residues
#==============================================================================

for i in cresids:

  #print 'Core residue = ', i
  fp = open('number' + str(i) + '.com', 'w')
  fp2 = open('number' + str(i) + '_ligand.com', 'w') 
  fpf = open('number' + str(i) + '.fpf', 'w')

  print >> fp, '%%chk=%s.chk' %('number' + str(i))
  print >> fp, '%mem=22GB'
  print >> fp, '%nproc=8'
  print >> fp, '#b3lyp/6-31G* charge nmr'
  print >> fp, ' '
  print >> fp, 'Have a nice day'
  print >> fp, ' '
  print >> fp, '%-5d %-5d' %(int(round(totchg, 0)), spin)

  print >> fp2, '%%chk=%s.chk' %('number' + str(i) + '_ligand')
  print >> fp2, '%mem=22GB'
  print >> fp2, '%nproc=8'
  print >> fp2, '#b3lyp/6-31G* charge nmr'
  print >> fp2, ' '
  print >> fp2, 'Have a nice day'
  print >> fp2, ' '
  print >> fp2, '%-5d %-5d' %(int(round(totchg2, 0)), spin2)

  bresids = [] #buffer region
  ccaps = [] #C terminal capped residue
  ncaps = [] #N terminal capped residue

  catoms = []
  im1 = i - 1

  if im1 != 0:
    for j in mol.residues[im1].resconter:
      atname = mol.atoms[j].atname
      if atname in ['C', 'O']:
        catoms.append(j)
  for j in mol.residues[i].resconter:
    atname = mol.atoms[j].atname
    if atname not in ['C', 'O']:
      catoms.append(j)

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
        res2 = mol.residues[mol.atoms[k].resid].resname
        if elmt2 == 'H':
          hnum2 = 1
        hnums = hnum1 + hnum2

        if hnums == 2 and dis <= 3.0:
          if mol.atoms[k].resid not in bresids:
            print mol.atoms[k].resid, "3.0"
            bresids.append(mol.atoms[k].resid)
        elif hnums == 1 and dis <= 4.0:
          if mol.atoms[k].resid not in bresids:
            print mol.atoms[k].resid, "4.0"
            bresids.append(mol.atoms[k].resid)
        elif (res2 in ['HIE', 'HIP', 'HID', 'HIS', 'PHE', 'TYR', 'TRP']) and \
          (at2 not in ['O', 'C', 'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OH', 'HO']) and \
          (dis <= 5.0):
          if mol.atoms[k].resid not in bresids:
            print mol.atoms[k].resid, "5.0"
            bresids.append(mol.atoms[k].resid)
  bresids.sort()
  print "bresids: ", bresids
  #Get the buffer region
  for j in range(i-2, i+3):
    if (j > 0) and (j not in bresids):
      bresids.append(j)

  #Get the disulfur bond
  for j in bresids:
    if mol.residues[j].resname == 'CYX':
      for k in mol.residues[j].resconter:
        if mol.atoms[k].atname == 'SG':
          for l in atids:
            if (l != k) and (mol.atoms[l].atname == 'SG') and \
              (mol.residues[mol.atoms[l].resid].resname == 'CYX'):
              dis = calc_bnd(mol.atoms[l].crd, mol.atoms[k].crd)
              if (dis <= 2.5) and (mol.atoms[l].resid not in bresids):
                bresids.append(mol.atoms[l].resid)

  #1. print the core region
  #print >> fp, "Core region residue: ", i
  for j in catoms:
    element = mol.atoms[j].element
    crdx = mol.atoms[j].crd[0]
    crdy = mol.atoms[j].crd[1]
    crdz = mol.atoms[j].crd[2]
    print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
    print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
    print >> fpf, str(mol.atoms[j].resid) + '-' + \
             mol.residues[mol.atoms[j].resid].resname + '-' + \
             mol.atoms[j].atname

  #0. print ligand
  for j in atids:
    if mol.residues[mol.atoms[j].resid].resname == 'LIG':
      element = mol.atoms[j].element
      crdx = mol.atoms[j].crd[0]
      crdy = mol.atoms[j].crd[1]
      crdz = mol.atoms[j].crd[2]
      print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)

  bresids.sort()
  #print 'Buffer residues = ', bresids

  for j in bresids:
    if (j-1 not in bresids) and (j-1 > 0):
      ncaps.append(j-1)

  for j in bresids:
    if (j+1 not in bresids):
      ccaps.append(j)

  bresids = list(set(bresids) - set(ccaps))
  bresids.sort()
 
  #for j in bresids:
  #  if (j+1 not in bresids):
  #    bresids.remove(j)

  batoms = []

  #print 'N-terminal Cap = ', ncaps, 'C-terminal Cap = ', ccaps, 
  #'Buffer residues = ', bresids

  #2. print the buffer region
  for j in ncaps:
    #print >> fp, "N-terminal cap residues: ", j
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
        print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
      elif mol.atoms[k].atname == 'CA': #change to H
        batoms.append(k)
        element = 'H'
        bvec = calc_bond(ccrd, mol.atoms[k].crd)
        crdx = ccrd[0] + bdld['CH'] * (mol.atoms[k].crd[0] - ccrd[0])/bvec
        crdy = ccrd[1] + bdld['CH'] * (mol.atoms[k].crd[1] - ccrd[1])/bvec
        crdz = ccrd[2] + bdld['CH'] * (mol.atoms[k].crd[2] - ccrd[2])/bvec
        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)
        print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)

  for j in ccaps:
    #print >> fp, "C-terminal cap residues: ", j
    for k in mol.residues[j].resconter:
      if mol.atoms[k].atname == 'CA':
        cacrd = mol.atoms[k].crd

    for k in mol.residues[j].resconter:
      if mol.atoms[k].atname != 'O':
        batoms.append(k)
        if mol.atoms[k].atname == 'C': #change to H
          element = 'H'
          bvec = calc_bond(cacrd, mol.atoms[k].crd)
          crdx = cacrd[0] + bdld['CH'] * (mol.atoms[k].crd[0] - cacrd[0])/bvec
          crdy = cacrd[1] + bdld['CH'] * (mol.atoms[k].crd[1] - cacrd[1])/bvec
          crdz = cacrd[2] + bdld['CH'] * (mol.atoms[k].crd[2] - cacrd[2])/bvec
          crdx = round(crdx, 3)
          crdy = round(crdy, 3)
          crdz = round(crdz, 3)
          print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
          print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        else:
          element = mol.atoms[k].element
          crdx = mol.atoms[k].crd[0]
          crdy = mol.atoms[k].crd[1]
          crdz = mol.atoms[k].crd[2]
          print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
          print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)

  for j in bresids:
    #print >> fp, "Buffer residues: ", j
    for k in mol.residues[j].resconter:
      if k not in catoms:
        batoms.append(k)
        element = mol.atoms[k].element
        crdx = mol.atoms[k].crd[0]
        crdy = mol.atoms[k].crd[1]
        crdz = mol.atoms[k].crd[2]
        print >> fp, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)
        print >> fp2, '%2s %10.3f %10.3f %10.3f' %(element, crdx, crdy, crdz)

  print >> fp, ' '
  print >> fp2, ' '

  #3. print the other parts:
  for j in atids:
    if (j not in catoms) and (j not in batoms):
      if mol.atoms[j].resid == 1:
        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        resname = mol.residues[mol.atoms[j].resid].resname
        atname = mol.atoms[j].atname
        chg = parmdict['N' + resname + '-' + atname][1]
        print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
        print >> fp2,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
      elif mol.atoms[j].resid == resids[-2]:
        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        resname = mol.residues[mol.atoms[j].resid].resname
        atname = mol.atoms[j].atname
        chg = parmdict['C' + resname + '-' + atname][1]
        print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
        print >> fp2,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
      elif mol.atoms[j].resid == resids[-1]:
        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        resname = mol.residues[mol.atoms[j].resid].resname
        atname = mol.atoms[j].atname
        try:
          chg = parmdict['C' + resname + '-' + atname][1]
        except:
          chg = parmdict[resname + '-' + atname][1]
        print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
        print >> fp2,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
      elif mol.residues[mol.atoms[j].resid].resname != 'LIG':
        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        resname = mol.residues[mol.atoms[j].resid].resname
        atname = mol.atoms[j].atname
        chg = parmdict[resname + '-' + atname][1]
        print >> fp,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)
        print >> fp2,'%10.3f %10.3f %10.3f %10.5f' %(crdx, crdy, crdz, chg)

  #print empty lines in Gaussian input files
  print >> fp, ' '
  print >> fp, ' '
  print >> fp, ' '
  print >> fp2, ' '
  print >> fp2, ' '
  print >> fp2, ' '

