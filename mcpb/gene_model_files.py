"""
This module was used to generate the pdb, fingerprint files of the sidechain,
standard and large models and Gaussian input file of sidechain (for bond and
angle parameter fitting) and large models(for RESP charge fitting).
"""
from __future__ import absolute_import
from pymsmtmol.readpdb import get_atominfo_fpdb, writepdbatm
from pymsmtmol.cal import calc_bond
from pymsmtmol.mol import pdbatm, gauatm, get_reslist
from pymsmtmol.element import (Atnum, CoRadiiDict,
                              get_ionljparadict, AtnumRev, bdld)
from pymsmtmol.gauio import (write_gauatm, write_gauatm_opth, write_gau_optf,
                             write_gau_fcf, write_gau_mkf)
from pymsmtmol.gmsio import (write_gmsatm, write_gms_optf, write_gms_fcf,
                             write_gms_mkf)
from pymsmtmol.sqmio import get_crdinfo_from_sqm, write_sqm_optf
from pymsmtlib.lib import get_lib_dict
import os

H_NAMES = ['HH31', 'HH32', 'HH33']  #hydrogen names for ACE and NME methyl group
SH_NAMES = ['H1', 'H2', 'H3'] #The names of the three Hs in the methyl group
GH_NAMES = ['HA2', 'HA3'] #GLY hydrogen names
SH_NAMES2 = ['H1', 'H2']

def del_files(fnamel):
    for fname in fnamel:
      if os.path.exists(fname):
        os.system("rm %s" %fname)

def count_lines(fname):
    ln = 0
    fp = open(fname, 'r')
    for line in fp:
      ln = ln + 1
    return ln
    fp.close()

#-------------------Get metal center residue names-----------------------------
def get_ms_resnames(pdbfile, ionids, cutoff):
    mol, atids, resids = get_atominfo_fpdb(pdbfile)
    ionids = ionids #metal ion atom id
    metresids = [] #metal ion residue id

    #Get the metal ion id
    for i in ionids:
      resid = mol.atoms[i].resid
      metresids.append(resid)

    msresids = [] #metal site residues
    msresids = msresids + metresids

    #Get the atoms which is in the cutoff of metal ion
    for met in ionids:
      for i in atids:
        if (i != met):
          dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
          if (dis <= cutoff) and mol.atoms[i].element != 'H':
            if (mol.atoms[i].resid not in msresids):
              msresids.append(mol.atoms[i].resid)

    msresids.sort()
    mcresnames = []
    tmpl = []
    for i in msresids:
      if len(mol.residues[i].resname) == 3:
        nresname = mol.residues[i].resname[0:3:2]
      elif len(mol.residues[i].resname) == 2:
        nresname = mol.residues[i].resname

      counter = 1
      for j in tmpl:
        if j == nresname:
          counter = counter + 1

      tmpl.append(nresname)
      nresname = nresname + str(counter)
      mcresnames.append(nresname)

    mcresnames0 = [mol.residues[i].resname for i in \
                   list(set(msresids)-set(metresids))]
    return mcresnames0, mcresnames

#--------------------Get metal site bonded atom ids--------------------------
def get_ms_ids(mol, atids, ionids, cutoff):

    bdatmids = []
    bdatnams = []

    #Get the atoms which is in the cutoff of metal ion
    for met in ionids:
      for i in atids:
        if (i != met):
          dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
          if (dis <= cutoff) and mol.atoms[i].element in ['O', 'N', 'S']:
            if i not in bdatmids:
              bdatmids.append(i)
              bdatnams.append(mol.atoms[i].atname)

    return bdatmids, bdatnams

#---------------------Write ACE residue into the PDB file---------------------
def write_ace(mol, i, gatms, pdbf, fpf=None):

    """
    ACE group
    a. C, O atoms are kept.
    b. CA --> CH3
    c. HA/HA2, CB/HA3, N --> Three Hs bonds to CH3 with adapting the \
       bond length.
    d. other atoms are deleted.

    Speical case:
    If there is a PRO was treated as ACE, there will be no influence.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    print "Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname +  " into ACE..."

    #get the coordinates of the CA atom
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if atname == 'CA':
        cacrd = mol.atoms[j].crd
    
    #rename the atom names to get the large model
    atnames = []
    hdict = {}
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname in ['HA', 'CB', 'N', 'HA2', 'HA3']):
        atnames.append(atname)

    for j in range(0, len(atnames)):
      hdict[atnames[j]] = H_NAMES[j]

    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      #C, O atoms will be still while CA is CH3, HA, CB and N
      #are three HH3s

      #Only for the backbone things
      if (atname in ['C', 'O', 'CA', 'HA', 'CB', 'N', 'HA2', 'HA3']):
        if (atname == 'C') or (atname == 'O'):
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = atname
        elif (atname == 'CA'):
          atname = 'CH3'
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'C'
        else: #left were HA, CB, C, HA2, HA3
          atname = hdict[atname]
          bvec = calc_bond(cacrd, mol.atoms[j].crd)
          crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
          crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
          crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
          element = 'H'

        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)

        #gausssian file
        gatms.append(gauatm(element, crdx, crdy, crdz))
    
        #assign other parameters to it
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = 'ACE'
        occp = 1.00
        tempfac = 0.00

        #pdb file
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, pdbf)

        #fingerprint file
        if fpf is not None:
          fpff = open(fpf, 'a')
          print >> fpff, str(resid) + '-' + 'ACE-' + atname
          fpff.close()

#---------------------Write CH3CO2- residue into the PDB file---------------------
def write_act(mol, i, gatms, pdbf, fpf=None):

    """
    ACT group
    a. C, O, OXT atoms are kept.
    b. CA --> CH3
    c. N, HA/HA2, CB/HA3 --> Three Hs bonds to CH3 with adapting the bond length.
    d. other atoms are deleted

    If there is a PRO was treated as ACT, there will be no influence.
    """

    global H_NAMES, SH_NAMES, GH_NAMES
   
    print "Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into ACT..."
   
    #get the coordinates of the CA atom
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if atname == 'CA':
        cacrd = mol.atoms[j].crd
   
    atnames = []
    hdict = {}
  
    #rename the atom names to get the large model
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname in ['N', 'HA', 'HA2', 'CB', 'HA3']):
        atnames.append(atname)
   
    for j in range(0, len(atnames)):
      hdict[atnames[j]] = H_NAMES[j]
   
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname in ['N', 'HA', 'HA2', 'CB', 'HA3', 'C', 'O', 'OXT', 'CA']):
        if (atname == 'C') or (atname == 'O'):
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = atname
        elif (atname == 'OXT'):
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'O'
        elif (atname == 'CA'):
          atname = 'CH3'
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'C'
        else: #left were HA, CB, C, HA2, HA3
          atname = hdict[atname]
          bvec = calc_bond(cacrd, mol.atoms[j].crd)
          crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
          crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
          crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
          element = 'H'
   
        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)
    
        #gaussian file
        gatms.append(gauatm(element, crdx, crdy, crdz))
    
        #assign other parameters to it
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = 'ACT'
        occp = 1.00
        tempfac = 0.00
   
        #pdb file
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, pdbf)

        #fingerprint file
        if fpf is not None:
          fpff = open(fpf, 'a')
          print >> fpff, str(resid) + '-' + 'ACT-' + atname
          fpff.close()

#---------------------Write NME residue into the PDB file---------------------
def write_nme(mol, i, gatms, pdbf, fpf=None):

    """
    NME group
    a. N, H atoms are kept.
    b. CA --> CH3
    c. HA/HA2, CB/HA3, C --> Three Hs bonds to CH3 with adapting the \
       bond length.
    d. other atoms are deleted.

    If there is a PRO was treated as NME, the atom CD need to change to \
       atom H.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    print "Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into NME..."

    #If the resname is not PRO
    if mol.residues[i].resname != 'PRO':

      #get the coordinates of the CA atom
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
          cacrd = mol.atoms[j].crd

      atnames = []
      hdict = {}
    
      #rename the atom names to get the large model
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['HA', 'CB', 'C', 'HA2', 'HA3']):
          atnames.append(atname)
 
      for j in range(0, len(atnames)):
        hdict[atnames[j]] = H_NAMES[j]
 
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'H', 'HN', 'CA', 'HA', 'CB', 'C', 'HA2', 'HA3']):
          #N, H atoms will be still 
          #while CA is CH3, HA, CB and C are three HH3s
          if (atname == 'N') or (atname == 'H'):
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = atname
          elif (atname == 'HN'):
            atname = 'H'
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = 'H'
          elif (atname == 'CA'):
            atname = 'CH3'
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = 'C'
          else: #left were HA, CB, C, HA2, HA3
            atname = hdict[atname]
            bvec = calc_bond(cacrd, mol.atoms[j].crd)
            crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
            crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
            crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
            element = 'H'
 
          crdx = round(crdx, 3)
          crdy = round(crdy, 3)
          crdz = round(crdz, 3)
    
          #gaussian file
          gatms.append(gauatm(element, crdx, crdy, crdz))
    
          #assign other parameters to it
          tiker = mol.atoms[j].gtype
          atid = mol.atoms[j].atid
          chainid = 'A'
          resid = mol.atoms[j].resid
          resname = 'NME'
          occp = 1.00
          tempfac = 0.00

          #pdb file
          atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                        crdx, crdy, crdz, occp, tempfac)
          writepdbatm(atmi, pdbf)

          #fingerprint file
          if fpf is not None:
            fpff = open(fpf, 'a')
            print >> fpff, str(resid) + '-' + 'NME-' + atname
            fpff.close()

    #If the resname is PRO
    else:

      #get the coordinates of the CA atom
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
          cacrd = mol.atoms[j].crd
        if atname == 'N':
          ncrd = mol.atoms[j].crd

      atnames = []
      hdict = {}
  
      #rename the atom names to get the large model
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['HA', 'CB', 'C']):
          atnames.append(atname)
 
      for j in range(0, len(atnames)):
        hdict[atnames[j]] = H_NAMES[j]
 
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'CA', 'HA', 'CB', 'C', 'CD']):
          #N atom will be still, CD will change to H
          #while CA is CH3, HA, CB and C are three HH3s
          if (atname == 'N'):
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = atname
          elif (atname == 'CA'):
            atname = 'CH3'
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = 'C'
          elif (atname == 'CD'):
            atname = 'H'
            bvec = calc_bond(ncrd, mol.atoms[j].crd)
            crdx = ncrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - ncrd[0])/bvec
            crdy = ncrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - ncrd[1])/bvec
            crdz = ncrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - ncrd[2])/bvec
            element = 'H'
          else: #left were HA, CB, C
            atname = hdict[atname]
            bvec = calc_bond(cacrd, mol.atoms[j].crd)
            crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
            crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
            crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
            element = 'H'
 
          crdx = round(crdx, 3)
          crdy = round(crdy, 3)
          crdz = round(crdz, 3)

          #gaussian file
          gatms.append(gauatm(element, crdx, crdy, crdz))
    
          #assign other parameters to it
          tiker = mol.atoms[j].gtype
          atid = mol.atoms[j].atid
          chainid = 'A'
          resid = mol.atoms[j].resid
          resname = 'NME'
          occp = 1.00
          tempfac = 0.00

          #pdb file
          atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                        crdx, crdy, crdz, occp, tempfac)
          writepdbatm(atmi, pdbf)

          #fingerprint file
          if fpf is not None:
            fpff = open(fpf, 'a')
            print >> fpff, str(resid) + '-' + 'NME-' + atname
            fpff.close()

#---------------------Write GLY residue into the PDB file---------------------
def write_gly(mol, i, gatms, pdbf, fpf=None):

    """
    GLY group
    a. N, H, C, O, CA are kept.
    b. HA, CB --> Two Hs bond to CA in GLY.
    c. other atoms are deleted.
    If there is a PRO was treated as GLY, the atom CD need to change to \
    atom H as well.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    #get the coordinates of the CA atom
    print "Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into GLY..."

    #If the resname is not PRO
    if mol.residues[i].resname != 'PRO':

      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
          cacrd = mol.atoms[j].crd

      #rename the atom names to get the large model
      atnames = []
      hdict = {}

      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['HA', 'CB']):
          atnames.append(atname)
 
      for j in range(0, len(atnames)):
        hdict[atnames[j]] = GH_NAMES[j]
 
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'H', 'HN', 'CA', 'HA', 'CB', 'C', 'O','HA2', 'HA3']):
          if atname in ['N', 'H', 'CA', 'C', 'O', 'HA2', 'HA3']:
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = mol.atoms[j].element
          elif (atname == 'HN'):
            atname = 'H'
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = 'H'
          else: #the left were HA, CB
            atname = hdict[atname]
            bvec = calc_bond(cacrd, mol.atoms[j].crd)
            crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
            crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
            crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
            element = 'H'
 
          crdx = round(crdx, 3)
          crdy = round(crdy, 3)
          crdz = round(crdz, 3)

          #gaussian file
          gatms.append(gauatm(element, crdx, crdy, crdz))
    
          #assign other parameters to it
          tiker = mol.atoms[j].gtype
          atid = mol.atoms[j].atid
          chainid = 'A'
          resid = mol.atoms[j].resid
          resname = 'GLY'
          occp = 1.00
          tempfac = 0.00

          #pdb file
          atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                        crdx, crdy, crdz, occp, tempfac)
          writepdbatm(atmi, pdbf)

          #fingerprint file
          if fpf is not None:
            fpff = open(fpf, 'a')
            print >> fpff, str(resid) + '-' + 'GLY-' + atname
            fpff.close()

    #If the resname is PRO
    else:

      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
          cacrd = mol.atoms[j].crd
        if atname == 'N':
          ncrd = mol.atoms[j].crd

      #rename the atom names to get the large model
      atnames = []
      hdict = {}
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['HA', 'CB']):
          atnames.append(atname)

      for j in range(0, len(atnames)):
        hdict[atnames[j]] = GH_NAMES[j]
 
      for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'CD', 'CA', 'HA', 'CB', 'C', 'O']):
          if atname in ['N', 'CA', 'C', 'O']:
            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            element = mol.atoms[j].element
          elif atname in ['HA', 'CB']:
            atname = hdict[atname]
            bvec = calc_bond(cacrd, mol.atoms[j].crd)
            crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
            crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
            crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
            element = 'H'
          else: #the only left is atom CD
            atname = 'H'
            bvec = calc_bond(ncrd, mol.atoms[j].crd)
            crdx = ncrd[0] + bdld['NH'] * (mol.atoms[j].crd[0] - ncrd[0])/bvec
            crdy = ncrd[1] + bdld['NH'] * (mol.atoms[j].crd[1] - ncrd[1])/bvec
            crdz = ncrd[2] + bdld['NH'] * (mol.atoms[j].crd[2] - ncrd[2])/bvec
            element = 'H'

          crdx = round(crdx, 3)
          crdy = round(crdy, 3)
          crdz = round(crdz, 3)

          #gaussian file
          gatms.append(gauatm(element, crdx, crdy, crdz))
    
          #assign other parameters to it
          tiker = mol.atoms[j].gtype
          atid = mol.atoms[j].atid
          chainid = 'A'
          resid = mol.atoms[j].resid
          resname = 'GLY'
          occp = 1.00
          tempfac = 0.00

          #pdb file
          atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                        crdx, crdy, crdz, occp, tempfac)
          writepdbatm(atmi, pdbf)
 
          #fingerprint file
          if fpf is not None:
            fpff = open(fpf, 'a')
            print >> fpff, str(resid) + '-' + 'GLY-' + atname
            fpff.close()

#---------------------Write normal residue into the PDB file---------------------
def write_normal(mol, reslist, i, gatms, pdbf, fpf=None):

    print "It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname

    for j in mol.residues[i].resconter:
      tiker = mol.atoms[j].gtype
      atid = mol.atoms[j].atid

      atname = mol.atoms[j].atname
      if (i in reslist.std) and (atname == 'HN'):
        atname = 'H'

      element = mol.atoms[j].element
      chainid = 'A'
      resid = mol.atoms[j].resid
      resname = mol.residues[resid].resname

      crdx = mol.atoms[j].crd[0]
      crdy = mol.atoms[j].crd[1]
      crdz = mol.atoms[j].crd[2]
      occp = 1.00
      tempfac = 0.00

      crdx = round(crdx, 3)
      crdy = round(crdy, 3)
      crdz = round(crdz, 3)

      #Gaussian file
      gatms.append(gauatm(element, crdx, crdy, crdz))

      #PDB file
      atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                    crdx, crdy, crdz, occp, tempfac)
      writepdbatm(atmi, pdbf)

      #Fingerprint file
      if fpf is not None:
        fpff = open(fpf, 'a')
        print >> fpff, str(resid) + '-' + resname + '-' + atname
        fpff.close()

#-----------------------Write Sidechain residues-------------------------------
def write_sc(mol, i, gatms, sidechf):

    global H_NAMES, SH_NAMES, GH_NAMES

    print "It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname

    #get the coordinates of the Ca atom
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if atname == 'CA':
        cacrd = mol.atoms[j].crd

    #N, C, HA are three Hs in the sidechain model
    atnames = []
    hdict = {}

    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname in ['N', 'C', 'HA']):
        atnames.append(atname)

    for j in range(0, len(atnames)):
      hdict[atnames[j]] = SH_NAMES[j]

    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname not in ['H', 'HN', 'O']): #HN is a alias of H
        #These two backbone atoms will be deleted in the sidechain modeling
        resname = mol.residues[i].resname
        if (atname == 'CA'):
          atname = 'CH3'
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'C'
        elif (atname in ['N', 'C', 'HA']):
          atname = hdict[atname]
          bvec = calc_bond(cacrd, mol.atoms[j].crd)
          crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
          crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
          crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
          element = 'H'
        else:
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = mol.atoms[j].element

        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)

        #Gaussian file
        gatms.append(gauatm(element, crdx, crdy, crdz))
    
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = mol.residues[resid].resname
        occp = 1.00
        tempfac = 0.00

        #PDB file
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, sidechf)

#-----------------------Write Sidechain residues2-------------------------------

#By keeping the N, H group inside the residue

def write_sc2(mol, i, gatms, sidechf):

    global H_NAMES, SH_NAMES, GH_NAMES, SH_NAMES2

    print "It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname

    #get the coordinates of the Ca atom
    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if atname == 'CA':
        cacrd = mol.atoms[j].crd

    #N, C, HA are three Hs in the sidechain model
    atnames = []
    hdict = {}

    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname in ['C', 'HA']):
        atnames.append(atname)

    for j in range(0, len(atnames)):
      hdict[atnames[j]] = SH_NAMES2[j]

    for j in mol.residues[i].resconter:
      atname = mol.atoms[j].atname
      if (atname not in ['O']):
        resname = mol.residues[i].resname
        if (atname == 'CA'):
          atname = 'CH3'
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'C'
        elif (atname in ['C', 'HA']):
          atname = hdict[atname]
          bvec = calc_bond(cacrd, mol.atoms[j].crd)
          crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
          crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
          crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
          element = 'H'
        elif atname == 'HN':
          atname = 'H'
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = 'H'
        else:
          crdx = mol.atoms[j].crd[0]
          crdy = mol.atoms[j].crd[1]
          crdz = mol.atoms[j].crd[2]
          element = mol.atoms[j].element

        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)

        #Gaussian file
        gatms.append(gauatm(element, crdx, crdy, crdz))
    
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = mol.residues[resid].resname
        occp = 1.00
        tempfac = 0.00

        #PDB file
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, sidechf)

#------------------------------Sidechain------------------------------------
def build_sidechain_model(mol, reslist, scresids, scresace, scresnme,
                          scresact, scresknh, totchg, outf, sqmopt):

    """
    For building the sidechain model
    1) Backbone O, H/HN was deleted.
    2) Backbone N, C, HA --> three Hs bonds to CA atom, with distance equals
       CH bond length.
    3) CA atom --> CH3 atom
    4) Other atoms is kept.

    Speical case: PRO does not have H atom, with CD atom instead, so there
    will be no influence since H was not considered in this modeling.
    """

    #Sidechain model file
    sidechf = outf + '_sidechain.pdb'

    #Gaussian
    goptf = outf + '_sidechain_opt.com'
    gfcf = outf + '_sidechain_fc.com'

    #GAMESS
    goptf2 = outf + '_sidechain_opt.inp'
    gfcf2 = outf + '_sidechain_fc.inp'

    #SQM
    siopf = outf + '_sidechain_sqm.in'
    soopf = outf + '_sidechain_sqm.out'

    #Delete the possible existing file
    del_files([sidechf, gfcf, goptf])

    #-------------------------------------------------------------------------
    ###############################Sidechain model############################
    #-------------------------------------------------------------------------

    print "***Creating the sidechain model..."

    gatms = [] #gaussian atom list

    for i in scresids:
      #1) For residue switching to ACE
      if i in scresace:
        write_ace(mol, i, gatms, sidechf)
      #2) For residue switching to NME
      elif i in scresnme:
        write_nme(mol, i, gatms, sidechf)
      #3) For residue switching to CH3CO2-
      elif i in scresact:
        write_act(mol, i, gatms, sidechf)
      #4) For residue which keep N and H in the model
      elif i in scresknh:
        write_sc2(mol, i, gatms, sidechf)
      #5) For normal amino acid residues
      elif i in reslist.std:
        write_sc(mol, i, gatms, sidechf)
      #6) For speical residue
      else:
        write_normal(mol, reslist, i, gatms, sidechf)

    ln = count_lines(sidechf)
    print "Totally there are " + str(ln) + " atoms in the sidechain model."

    #Calculate the spin number and print it into gaussian file
    gaelemts = 0
    for gatm in gatms:
      AtNum = Atnum[gatm.element]
      gaelemts = gaelemts + AtNum

    SpinNum = gaelemts - totchg
    SpinNum = int(round(SpinNum, 0))
    print "Totally there are " + str(SpinNum) + " electrons in the sidechain model."

    if SpinNum%2 == 0:
      SpinNum = 1
    else:
     SpinNum = 2

    #Gaussian
    write_gau_optf(outf, goptf, totchg, SpinNum, gatms)
    write_gau_fcf(outf, gfcf)

    #GAMESS
    write_gms_optf(goptf2, totchg, SpinNum, gatms)
    write_gms_fcf(gfcf2, totchg, SpinNum)

    #Perform the SQM calcualtion under PM6 first
    if (sqmopt == 1) or (sqmopt == 3):
      #Delete the possible existing file
      del_files([siopf, soopf])
      write_sqm_optf(siopf, totchg, gatms)
      if SpinNum == 1:
        print "Performing SQM optimization of sidechain model, please wait..."
        #Run SQM to optimize the coordinates
        os.system("sqm -i %s -o %s" %(siopf, soopf))
        gatms2 = get_crdinfo_from_sqm(soopf)
        write_gau_optf(outf, goptf, totchg, SpinNum, gatms2, 4)
        write_gms_optf(goptf2, totchg, SpinNum, gatms2, 4)
      else:
        print "Could not perform SQM optimization for the sidechain model " + \
              "with spin number not equal to 1."

#------------------------------------Standard model---------------------------
def build_standard_model(mol, reslist, cutoff, msresids, outf, ionids,
                         bdedatms, libdict, autoattyp):

    #Standard model file
    stf = outf + '_standard.pdb'
    stpf = outf + '_standard.fingerprint'
    del_files([stf, stpf])

    #-------------------------------------------------------------------------
    ###############################Standard model#############################
    #-------------------------------------------------------------------------

    print "***Creating the standard model..."

    for i in msresids:
      print "It contains the residue " + str(i) + '-' + \
            mol.residues[i].resname

      for j in mol.residues[i].resconter:

        atname = mol.atoms[j].atname
        if i in reslist.std and atname == 'HN':
          atname = 'H'

        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = mol.residues[resid].resname
        occp = 1.00
        tempfac = 0.00

        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)

        attype = libdict[resname + '-' + atname][0]
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, stf)
        stpff = open(stpf, 'a')

        #assign new atom types to atoms inside the metal site
        attype2 = attype

        if autoattyp == 1:
          if j in ionids:
            ionresn = mol.atoms[j].resname
            ionatn = mol.atoms[j].atname
            ionelmt = mol.atoms[j].element
            ionchg = libdict[ionresn + '-' + ionatn][1]
            ionchg = int(ionchg)
            if ionchg > 0:
              attype2 = ionelmt + str(ionchg) + '+'
            elif ionchg < 0:
              attype2 = ionelmt + str(ionchg) + '-'
            else:
              raise pymsmtError('Ion Should not have 0 charge on it.')
        elif autoattyp == 2:
          #for metal ion
          for k in range(0, len(ionids)):
            if j == ionids[k]:
              k2 = k + 1
              attype2 = 'Z' + str(k2)
          #for bonded atoms
          for k in range(0, len(bdedatms)):
            if j == bdedatms[k] and k < 9:
              k3 = k + 1
              attype2 = 'X' + str(k3)
            elif j == bdedatms[k] and k >= 9:
              k4 = k - 9 + 1
              attype2 = 'Y' + str(k4)

        print >> stpff, str(resid) + '-' + resname + '-' + atname, \
                 str(atid), attype, '->', attype2
        stpff.close()

    #Print the link information into sidechain fingerprint file
    stpff = open(stpf, 'a')
    for met in ionids:
      for i in bdedatms:
        dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
        if (dis <= cutoff):
          print >> stpff, "LINK", str(met)+'-'+mol.atoms[met].atname, \
                   str(i)+'-'+mol.atoms[i].atname
    stpff.close()

    ln = count_lines(stf)
    print "Totally there are " + str(ln) + " atoms in the standard model."

#---------------------------------Large model---------------------------------
def build_large_model(mol, reslist, lmsresids, lmsresace, lmsresnme,
                      lmsresgly, ionids, chargedict, totchg, outf,
                      watermodel, sqmopt):

    #Large model file
    largef = outf + '_large.pdb'
    lfpf = outf + '_large.fingerprint'
    gmkf = outf + '_large_mk.com'
    gmsf = outf + '_large_mk.inp'
    simkf = outf + '_large_sqm.in'
    somkf = outf + '_large_sqm.out'
    del_files([largef, lfpf, gmkf])

    #-------------------------------------------------------------------------
    ###############################Large model################################
    #-------------------------------------------------------------------------

    print "***Creating the large model..."
    gatms = []
    for i in lmsresids:
      #1) for atoms in ACE ---------------------------------------------------
      if i in lmsresace:
        write_ace(mol, i, gatms, largef, lfpf)
      #2) for atoms in NME ---------------------------------------------------
      elif i in lmsresnme:
        write_nme(mol, i, gatms, largef, lfpf)
      #3) for atoms in GLY ---------------------------------------------------
      elif i in lmsresgly:
        write_gly(mol, i, gatms, largef, lfpf)
      #4) for atoms in other residues ----------------------------------------
      else:
        write_normal(mol, reslist, i, gatms, largef, lfpf)

    ln = count_lines(largef)
    print "Totally there are " + str(ln) + " atoms in the large model."

    #Calculate the spin number and print it into gaussian file
    gaelemts = 0
    for gatm in gatms:
      AtNum = Atnum[gatm.element]
      gaelemts = gaelemts + AtNum

    SpinNum = gaelemts - totchg
    SpinNum = int(round(SpinNum, 0))
    print "Totally there are " + str(SpinNum) + " electrons in the large model."

    if SpinNum%2 == 0:
      SpinNum = 1
    else:
      SpinNum = 2

    #For Gaussian file
    IonLJParaDict = get_ionljparadict(watermodel)
    ionnames = [mol.atoms[i].atname for i in ionids]
    ionnames = list(set(ionnames))
    write_gau_mkf(outf, gmkf, totchg, SpinNum, gatms, ionnames, chargedict, IonLJParaDict)

    #For GAMESS file
    write_gms_mkf(gmsf, totchg, SpinNum, gatms)

    #-------------------------------------------------------------------------
    # Doing SQM Optimization
    #-------------------------------------------------------------------------
    if (sqmopt == 2) or (sqmopt == 3):
      del_files([simkf, somkf])
      write_sqm_optf(simkf, totchg, gatms)
      if SpinNum == 1:
        print "Performing SQM optimization of large model, please wait..."
        os.system("sqm -i %s -o %s" %(simkf, somkf))
        gatms2 = get_crdinfo_from_sqm(somkf)
        write_gau_mkf(outf, gmkf, totchg, SpinNum, gatms, ionnames, chargedict, IonLJParaDict, 4)
        write_gms_mkf(gmsf, totchg, SpinNum, gatms2, 4)
      else:
        print "Could not perform SQM optimization for the large model " + \
              "with spin number not equal to 1."

def gene_model_files(pdbfile, ionids, outf, ffchoice, naamol2f, cutoff, \
                     watermodel, autoattyp, sqmopt):

    mol, atids, resids = get_atominfo_fpdb(pdbfile)

    reslist = get_reslist(mol, resids)

    libdict, chargedict = get_lib_dict(ffchoice)

    for mol2f in naamol2f:
      libdict1, chargedict1 = get_lib_dict(mol2f)
      libdict.update(libdict1)
      chargedict.update(chargedict1)

    #-------------------------------------------------------------------------
    # Get the residues in the metal site
    #-------------------------------------------------------------------------

    print "******************************************************************"
    print "*                                                                *"
    print "*=======================Metal Site Information===================*"
    print "*                                                                *"
    print "******************************************************************"

    #1. Metal ions information
    metresids = [] #metal ion residue id
    for i in ionids:
      resid = mol.atoms[i].resid
      metresids.append(resid)
      print "***Selected Metal ion " + mol.atoms[i].atname + " is atom " + \
            str(i) + " in residue " + str(mol.atoms[i].resid) + '-' + \
            mol.residues[resid].resname
    ionids = list(set(ionids))
    ionids.sort()

    #2. Metal site residues information
    msresids = [] #metal site residues
    bdedatms, bdedatnams = get_ms_ids(mol, atids, ionids, cutoff)

    #3. Get the metal site containing residues
    for i in bdedatms:
      print str(mol.atoms[i].resid) + '-' + \
            mol.atoms[i].resname + \
            '@' + mol.atoms[i].atname + ' is in ' + str(cutoff) + \
            ' Angstrom of these metal ions' #+
            #str(mol.atoms[i].resid) + '-' + \
            #mol.atoms[i].resname + '@' + \
            #mol.atoms[i].atname
      if mol.atoms[i].resid not in msresids:
        msresids.append(mol.atoms[i].resid)

    msresids = msresids + metresids
    msresids.sort()

    print "***The following residues are in the Metal Site:"
    totchg = 0.0
    for i in msresids:
      print "Residue " + str(i) + '-' + mol.residues[i].resname
      totchg = totchg + chargedict[mol.residues[i].resname]
 
    #-------------------------------------------------------------------------
    # Get the residues for building the sidechain model and print
    #-------------------------------------------------------------------------

    scresids = msresids
    scresace = []
    scresnme = []
    scresact = []
    scresknh = [] #Residues to keep n and h, which is connect to the residue
                  #which has backbone oxygen bond to the ion and also bond to
                  #to the ion but with sidechain

    bdedresids = []
    bdedresdict = {}

    for i in range(0, len(bdedatms)):
      atm = bdedatms[i]
      atname = bdedatnams[i]
      resid = mol.atoms[atm].resid
      if resid not in bdedresids:
        bdedresids.append(resid)

    for resid in bdedresids:
      resatns = []
      for i in bdedatms:
        if mol.atoms[i].resid == resid:
          resatns.append(mol.atoms[i].atname)
      bdedresdict[resid] = resatns

    for resid in bdedresids:
      resatns = bdedresdict[resid]
      if resid in reslist.cterm:
        scresact.append(resid)
      elif (resid in reslist.std) and ('O' in resatns):
        scresace.append(resid)
        if (resid+1 in reslist.std) and (resid+1 not in scresids):
          scresids.append(resid+1)
          scresnme.append(resid+1)
        elif (resid+1 in reslist.std) and (resid+1 in scresids):
          scresknh.append(resid+1)

    print "***The sidechain model contains the following residues: "
    print scresids

    #-------------------------------------------------------------------------
    # Get the residues for building the large model and print
    #-------------------------------------------------------------------------

    lmsresids = [] #large model metal site residues
    lmsresace = [] #large model metal site ACE
    lmsresnme = [] #large model metal site NME
    lmsresgly = [] #large model metal site GLY

    #First several residues
    for i in range(0, len(msresids)-1):
      resi = msresids[i]
      resj = msresids[i+1]
      resnamei = mol.residues[resi].resname
      resnamej = mol.residues[resj].resname
      if resi in reslist.std:
        if (resj in reslist.std) and (resj - resi <= 5):
        #If two residues within 5 residues apart, treat the residue as GLY
          for j in resids:
            if (j > resi) and (j < resj) and (j not in msresids):
              lmsresgly.append(j) #GLY

    for i in range(0, len(msresids)):
      resid = msresids[i]
      resname = mol.residues[resid].resname
      if resid in reslist.std:
        if (resid-1 not in lmsresgly) and (resid-1 not in msresids):
          lmsresace.append(resid-1) #ACE
        if (resid+1 not in lmsresgly) and (resid+1 not in msresids):
          lmsresnme.append(resid+1) #NME

    lmsresids = msresids + lmsresace + lmsresnme + lmsresgly #Combine the residues
    lmsresids = list(set(lmsresids)) #Delete repeat elements
    lmsresids = sorted(lmsresids) #Sort the list

    print "***The large model contains the following residues: "
    print lmsresids

    #-------------------------------------------------------------------------
    # Generate model files
    #-------------------------------------------------------------------------

    print "******************************************************************"
    print "*                                                                *"
    print "*=======================Building models==========================*"
    print "*                                                                *"
    print "******************************************************************"

    build_sidechain_model(mol, reslist, scresids, scresace, scresnme,
                          scresact, scresknh, totchg, outf, sqmopt)

    build_standard_model(mol, reslist, cutoff, msresids, outf, ionids,
                         bdedatms, libdict, autoattyp)

    build_large_model(mol, reslist, lmsresids, lmsresace, lmsresnme, lmsresgly,
                      ionids, chargedict, totchg, outf, watermodel, sqmopt)

    #Using the automatically detect bond method for the backup
    #else:
    #  for met in ionids:
    #    Radiusmet = CoRadiiDict[mol.atoms[met].element]
    #    for i in atids:
    #      Radiusi = CoRadiiDict[mol.atoms[i].element]
    #      if (i != met):
    #       dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
    #        cutoff = Radiusmet + Radiusi + 0.40
    #        if (dis <= cutoff) and (dis > 0.1):
    #          print str(mol.atoms[i].resid) + '-' + \
    #                mol.atoms[i].resname + \
    #                '@' + mol.atoms[i].atname + ' is in ' + str(cutoff) + \
    #                ' Angstrom of ' + str(mol.atoms[met].resid) + '-' + \
    #                mol.atoms[met].resname + '@' + \
    #                mol.atoms[met].atname
    #          if (mol.atoms[i].resid not in msresids):
    #            msresids.append(mol.atoms[i].resid)
