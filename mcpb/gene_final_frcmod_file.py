"""
There are three methods in this module to get the force constants of the bond,
angle for the metal complexes.

1) The empirical method, which is developed by Pengfei Li and co-workers in
Merz Research Group, don't need to do the QM calculation to get the
parameters, is a convenient way to generate the meaningful parameters.

2) Seminario method, which is a method developed by Jorge M. Seminario, using
the sub-matrix of the Hessian Matrix to generate the force constants.

3) Z-matrix method, which uses the whole Hessian Matrix to generate the force
constants.
"""
from __future__ import absolute_import
from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtmol.getlist import get_alist, get_mc_blist
from pymsmtmol.gauio import (get_crds_from_fchk, get_matrix_from_fchk,
                             get_fc_from_log)
from pymsmtmol.gmsio import get_crds_from_gms, get_matrix_from_gms
from pymsmtmol.cal import calc_bond, calc_angle, calc_dih
from pymsmtmol.element import ionnamel
from pymsmtmol.constants import *
from pymsmtlib.lib import getfc
from pymsmtexp import *
from numpy import average, array, dot, cross, std
from numpy.linalg import eigvals, eig, norm
import math

#-----------------------------------------------------------------------------
# Related fuctions
#-----------------------------------------------------------------------------

def get_attypdict(stfpf, atids):

    #From standard fingerprint file to get the atomtype information
    attypdict = {}
    fpinfo = open(stfpf, 'r')
    for line in fpinfo:
      if line[0:4] != "LINK":
        line = line.strip('\n')
        line = line.split(' ')
        line = [i for i in line if i != '']
        if len(line[-1]) == 1:
          line[-1] = line[-1] + ' '
        longatname = line[0].split('-')
        atname = longatname[-1]
        if (int(line[1]) in atids):
          attypdict[int(line[1])] = line[-1]
    fpinfo.close()

    for i in atids:
      if i not in attypdict.keys():
        #print 'atomid %d Not find in standard fingerprint file.' %i
        attypdict[i] = '**'

    return attypdict

def get_misstyps(pref):

    #pre_frcmod file
    missbondtyps = []
    missangtyps = []

    prefrcmod = open(pref, 'r')
    for line in prefrcmod:
      line = line.strip('\n')
      if line[0:4] == 'NON ':
        if (line[6:7] == '-') and (line[9:10] == '-'):
          at1 = line[4:6]
          at2 = line[7:9]
          at3 = line[10:12]
          missangtyps.append((at1, at2, at3))
        else:
          at1 = line[4:6]
          at2 = line[7:9]
          missbondtyps.append((at1, at2))
    prefrcmod.close()
    return missbondtyps, missangtyps

def avg_bond_para(misbondat12, bondlen, bfconst):
    # bond length has 4 decimal places
    # bond force constant has 1 decimal place
    if (len(bondlen) == 0) or (len(bfconst) == 0):
      print 'For bondtype: ' + misbondat12 + '. There are ' + \
            ' 0 bond in this type: '
      print '  Could not find bond parameters for %s, treat it as 0.' \
            %misbondat12
      avgbfc = 0.0
      avgblen = 0.0
    else:
      bondlen = [round(i, 4) for i in bondlen]
      avgblen = round(average(bondlen), 4)
      bfconst = [round(i, 1) for i in bfconst]
      avgbfc = round(average(bfconst), 1)

      print 'For bondtype: ' + misbondat12 + '. There are ' + \
             str(len(bondlen)) + ' bond(s) in this type: '
      print '  The bond length(s) is(are): ' + str(bondlen)
      print '  The average bond length is: ' + str(avgblen) + \
            ' Angstrom.'
      print '  The force constant(s) is(are): ' + str(bfconst)
      print '  The average force constant is: ' + str(avgbfc) + \
            ' Kcal*mol^-1*A^-2.'

      if avgbfc < 0:
        print '  *The average force constant is negative, treat is as 0.0.'
        avgbfc = 0.0

    bond_para = (avgbfc, avgblen)
    return bond_para

def avg_angle_para(misangat123, angvals, afconst):
    # angle mangnitude has 2 decimal places
    # angle force constant has 2 decimal places
    if len(angvals) == 0 or len(afconst) == 0:
      print 'For angletype: ' + misangat123 + '. There are ' + \
            '0 angle in this type: '
      print '  Could not find angle parameters for %s, treat it as 0.' \
            %misangat123
      avgafc = 0.0
      avgangv = 0.0
    else:
      angvals = [round(i, 2) for i in angvals]
      avgangv = round(average(angvals), 2)
      afconst = [round(i, 2) for i in afconst]
      avgafc = round(average(afconst), 2)
 
      print 'For angletype: ' + misangat123 + '. There are ' + \
            str(len(angvals)) + ' angle(s) in this type: '
      print '  The angle value(s) is(are): ' + str(angvals)
      print '  The average angle value is: ' + str(avgangv) + \
            ' Degree.'
      print '  The force constant(s) is(are): '+ str(afconst)
      print '  The average force constant is: ' + str(avgafc) + \
            ' Kcal*mol^-1*Rad^-2.'

      if avgafc < 0:
        print '  *The average force constant is negative, treat is as 0.00.'
        avgafc = 0.00

    ang_para = (avgafc, avgangv)
    return ang_para

def print_frcmod_file(pref, finf, finalparmdict, method):

    ###print the final frcmod file
    finfrcmod = open(finf, 'w')
    prefrcmod = open(pref, 'r')
    note = 'Created by ' + method + ' method using MCPB.py'
    for line in prefrcmod:
      line = line.strip('\n')
      if line[0:4] == 'NON ':
        #Angle
        if (line[6:7] == '-') and (line[9:10] == '-'):
          at1 = line[4:6]
          at2 = line[7:9]
          at3 = line[10:12]
          angtyp = line[4:12]
          param = finalparmdict[angtyp]
          print >> finfrcmod, '%8s  %7.2f    %7.2f    %-s' %(angtyp,
                                             param[0], param[1], note)
        #Bond
        else:
          at1 = line[4:6]
          at2 = line[7:9]
          bondtyp = line[4:9]
          param = finalparmdict[bondtyp]
          print >> finfrcmod, '%5s  %5.1f   %7.4f      %-s' %(bondtyp,
                                              param[0], param[1], note)
      elif line[0:4] == 'YES ':
        print >> finfrcmod, line[4:]
      else:
        print >> finfrcmod, line
    prefrcmod.close()
    finfrcmod.close()

#-----------------------------------------------------------------------------
##Empirical method Ref: Li et al. In preparation.
#-----------------------------------------------------------------------------

def fcfit_ep_bond(dis, elmts):
    """Use the empirical method developed by Li et al. in Merz Research Group
    """

    if 'Zn' in elmts:
      if set(elmts) == set(['Zn', 'N']):
        fc = getfc('Zn-N.txt', dis)
      elif set(elmts) == set(['Zn', 'O']):
        fc = getfc('Zn-O.txt', dis)
      elif set(elmts) == set(['Zn', 'S']):
        fc = getfc('Zn-S.txt', dis)
      else:
        raise pymsmtError('Error: Could not deal with the atoms out of '
                          'N, O and S.')
    else:
      raise pymsmtError('Error: Not support metal ion other than Znic in '
                        'current version, more metal ions will be supported '
                        'in future version.')
    fc = round(fc, 1)
    return fc

def fcfit_ep_angle(elmts):
    #if metal is in the center
    if elmts[1] in ionnamel:
      return 35.00
    elif set(elmts[0:2]) == set(['Zn', 'S']) or set(elmts[1:3]) == set(['Zn', 'S']):
      return 70.00
    elif set(elmts[0:2]) == set(['Zn', 'N']) or set(elmts[1:3]) == set(['Zn', 'N']):
      return 50.00
    elif set(elmts[0:2]) == set(['Zn', 'O']) or set(elmts[1:3]) == set(['Zn', 'O']):
      return 50.00

def gene_by_empirical_way(scpdbf, ionids, stfpf, pref, finf):

    print "******************************************************************"
    print "*                                                                *"
    print "*===========Using the empirical way to solve the problem=========*"
    print "*                                                                *"
    print "******************************************************************"

    #Read from sidechain pdb
    mol, atids, resids = get_atominfo_fpdb(scpdbf)
    blist = get_mc_blist(mol, atids, ionids, stfpf)
    alist = get_alist(mol, blist)

    attypdict = get_attypdict(stfpf, atids)
    missbondtyps, missangtyps = get_misstyps(pref)

    finalparmdict = {}

    #for bond
    print "=======================Generate the bond parameters==============="

    for misbond in missbondtyps:
      bondlen = []
      bfconst = []
      for bond in blist:
        at1 = bond[0]
        at2 = bond[1]
        bondtyp = (attypdict[at1], attypdict[at2])
        if bondtyp == misbond or bondtyp[::-1] == misbond:
          crd1 = mol.atoms[at1].crd
          crd2 = mol.atoms[at2].crd
          dis = calc_bond(crd1, crd2)
          dis = round(dis, 4)
          bondlen.append(dis)

          #get the atom number
          elmt1 = mol.atoms[at1].element
          elmt2 = mol.atoms[at2].element
          elmts = [elmt1, elmt2]
          elmts = sorted(elmts)

          fcfinal = fcfit_ep_bond(dis, elmts)
          bfconst.append(fcfinal)

      #Get average bond parameters
      misbondat12 = misbond[0] + '-' + misbond[1]
      bond_para = avg_bond_para(misbondat12, bondlen, bfconst)
      #get the force constant
      finalparmdict[misbondat12] = bond_para

    #for angle
    print "=======================Generate the angle parameters=============="

    for misang in missangtyps:
      angvals = []
      afconst = []
      for ang in alist:
        at1 = ang[0]
        at2 = ang[1]
        at3 = ang[2]
        angtyp = (attypdict[at1], attypdict[at2], attypdict[at3])
        if angtyp == misang or angtyp[::-1] == misang:
          crd1 = mol.atoms[at1].crd
          crd2 = mol.atoms[at2].crd
          crd3 = mol.atoms[at3].crd
          angval = calc_angle(crd1, crd2, crd3)
          angvals.append(angval)

          elmt1 = mol.atoms[at1].element
          elmt2 = mol.atoms[at2].element
          elmt3 = mol.atoms[at3].element
          elmts = [elmt1, elmt2, elmt3]
          fcfinal = fcfit_ep_angle(elmts)
          afconst.append(fcfinal)

      #Get average angle parameters
      misangat123 = misang[0] + '-' + misang[1] + '-' + misang[2]
      ang_para = avg_angle_para(misangat123, angvals, afconst)
      finalparmdict[misangat123] = ang_para

    #Print out the final frcmod file
    print_frcmod_file(pref, finf, finalparmdict, 'empirical')

#-----------------------------------------------------------------------------

##Seminario method Ref: Calculation of intramolecular force fields from
##second-derivative tensors
##Seminario, International Journal of Quantum Chemistry 1996, 60(7), 1271-1277

#-----------------------------------------------------------------------------

def get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2, scalef, bondavg):

    crd1 = crds[3*nat1-3:3*nat1]
    crd2 = crds[3*nat2-3:3*nat2]
    disbohr = calc_bond(crd1, crd2) #unit is bohr
    dis = disbohr * B_TO_A #Transfer bohr to angstrom

    vec12 = array(crd2) - array(crd1) #vec12 is vec2 - vec1
    vec12 = [i/(disbohr) for i in vec12]
    vec12 = array(vec12)
 
    #bond force constant matrix, size 3 * 3
    bfcmatrix = array([[float(0) for x in range(3)] for x in range(3)])
   
    #1. First way to chose the matrix-----------------
    for i in range(0, 3):
      for j in range(0, 3):
        bfcmatrix[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat2-1)+j]
    eigval, eigvector = eig(bfcmatrix)
    fc = 0.0
    for i in range(0, 3):
      ev = eigvector[:,i]
      fc = fc + eigval[i] * abs(dot(ev, vec12))
    fcfinal1 = fc * HB2_TO_KCAL_MOL_A2 * 0.5

    if bondavg == 1:
      #2. Second way to chose the matrix-----------------
      for i in range(0, 3):
        for j in range(0, 3):
          bfcmatrix[i][j] = -fcmatrix[3*(nat2-1)+i][3*(nat1-1)+j]
      eigval, eigvector = eig(bfcmatrix)
      fc = 0.0
      for i in range(0, 3):
        ev = eigvector[:,i]
        fc = fc + eigval[i] * abs(dot(ev, vec12))
      fcfinal2 = fc * HB2_TO_KCAL_MOL_A2 * 0.5

      #Hatree/(Bohr^2) to kcal/(mol*angstrom^2)
      #Times 0.5 factor since AMBER use k(r-r0)^2 but not 1/2*k*(r-r0)^2

      fcfinal = average([fcfinal1, fcfinal2])
      stdv = std([fcfinal1, fcfinal2])
      fcfinal = fcfinal * scalef * scalef
      stdv = stdv * scalef * scalef
      return dis, fcfinal, stdv

    elif bondavg == 0:

      fcfinal = fcfinal1 * scalef * scalef
      return dis, fcfinal

def get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, scalef, angavg):

    #get the angle value
    crd1 = crds[3*nat1-3:3*nat1]
    crd2 = crds[3*nat2-3:3*nat2]
    crd3 = crds[3*nat3-3:3*nat3]
    dis12 = calc_bond(crd1, crd2) #unit is bohr
    dis32 = calc_bond(crd3, crd2) #unit is bohr

    #get the unit vector
    vec12 = array(crd2) - array(crd1) #vec12 is vec2 - vec1
    vec32 = array(crd2) - array(crd3)
    vec12 = array([i/dis12 for i in vec12])
    vec32 = array([i/dis32 for i in vec32])

    angval = calc_angle(crd1, crd2, crd3)

    #get the normalized vector
    vecUNp = cross(vec32, vec12)
    vecUN = array([i/norm(vecUNp) for i in vecUNp]) #vecUN is the vector
                                     #perpendicular to the plance of ABC

    vecPA = cross(vecUN, vec12)
    vecPC = cross(vec32, vecUN)

    afcmatrix12 = array([[float(0) for x in range(3)] for x in range(3)])
    afcmatrix32 = array([[float(0) for x in range(3)] for x in range(3)])

    #1. First way to chose the matrix----------------------------------
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix12[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat2-1)+j]
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix32[i][j] = -fcmatrix[3*(nat3-1)+i][3*(nat2-1)+j]
    eigval12, eigvector12 = eig(afcmatrix12)
    eigval32, eigvector32 = eig(afcmatrix32)
    contri12 = 0.0
    contri32 = 0.0
    for i in range(0, 3):
      ev12 = eigvector12[:,i]
      ev32 = eigvector32[:,i]
      contri12 = contri12 + eigval12[i] * abs(dot(vecPA, ev12))
      contri32 = contri32 + eigval32[i] * abs(dot(vecPC, ev32))
    contri12 = 1.0 / (contri12 * dis12 * dis12)
    contri32 = 1.0 / (contri32 * dis32 * dis32)
    fcfinal1 = (1.0 / (contri12 + contri32)) * H_TO_KCAL_MOL * 0.5

    if angavg == 1:
      #2. Second way to chose the matrix----------------------------------
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix12[i][j] = -fcmatrix[3*(nat2-1)+i][3*(nat1-1)+j]
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix32[i][j] = -fcmatrix[3*(nat3-1)+i][3*(nat2-1)+j]
      eigval12, eigvector12 = eig(afcmatrix12)
      eigval32, eigvector32 = eig(afcmatrix32)
      contri12 = 0.0
      contri32 = 0.0
      for i in range(0, 3):
        ev12 = eigvector12[:,i]
        ev32 = eigvector32[:,i]
        contri12 = contri12 + eigval12[i] * abs(dot(vecPA, ev12))
        contri32 = contri32 + eigval32[i] * abs(dot(vecPC, ev32))
      contri12 = 1.0 / (contri12 * dis12 * dis12)
      contri32 = 1.0 / (contri32 * dis32 * dis32)
      fcfinal2 = (1.0 / (contri12 + contri32)) * H_TO_KCAL_MOL * 0.5
 
      #Hatree to kcal/mol
      #Times 0.5 factor since AMBER use k(r-r0)^2 but not 1/2*k*(r-r0)^2
 
      #3. Third way to chose the matrix----------------------------------
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix12[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat2-1)+j]
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix32[i][j] = -fcmatrix[3*(nat2-1)+i][3*(nat3-1)+j]
      eigval12, eigvector12 = eig(afcmatrix12)
      eigval32, eigvector32 = eig(afcmatrix32)
      contri12 = 0.0
      contri32 = 0.0
      for i in range(0, 3):
        ev12 = eigvector12[:,i]
        ev32 = eigvector32[:,i]
        contri12 = contri12 + eigval12[i] * abs(dot(vecPA, ev12))
        contri32 = contri32 + eigval32[i] * abs(dot(vecPC, ev32))
      contri12 = 1.0 / (contri12 * dis12 * dis12)
      contri32 = 1.0 / (contri32 * dis32 * dis32)
      fcfinal3 = (1.0 / (contri12 + contri32)) * H_TO_KCAL_MOL * 0.5

      #4. Fourth way to chose the matrix----------------------------------
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix12[i][j] = -fcmatrix[3*(nat2-1)+i][3*(nat1-1)+j]
      for i in range(0, 3):
        for j in range(0, 3):
          afcmatrix32[i][j] = -fcmatrix[3*(nat2-1)+i][3*(nat3-1)+j]
      eigval12, eigvector12 = eig(afcmatrix12)
      eigval32, eigvector32 = eig(afcmatrix32)
      contri12 = 0.0
      contri32 = 0.0
      for i in range(0, 3):
        ev12 = eigvector12[:,i]
        ev32 = eigvector32[:,i]
        contri12 = contri12 + eigval12[i] * abs(dot(vecPA, ev12))
        contri32 = contri32 + eigval32[i] * abs(dot(vecPC, ev32))
      contri12 = 1.0 / (contri12 * dis12 * dis12)
      contri32 = 1.0 / (contri32 * dis32 * dis32)
      fcfinal4 = (1.0 / (contri12 + contri32)) * H_TO_KCAL_MOL * 0.5
 
      fcfinal = average([fcfinal1, fcfinal2, fcfinal3, fcfinal4])
      stdv = std([fcfinal1, fcfinal2, fcfinal3, fcfinal4])
      fcfinal = fcfinal * scalef * scalef
      stdv = stdv * scalef * scalef
      return angval, fcfinal, stdv

    elif angavg == 0:
      fcfinal = fcfinal1 * scalef * scalef
      return angval, fcfinal

def get_dih_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, nat4, n1, n2, scalef):

    crd1 = crds[3*nat1-3:3*nat1]
    crd2 = crds[3*nat2-3:3*nat2]
    crd3 = crds[3*nat3-3:3*nat3]
    crd4 = crds[3*nat4-3:3*nat4]

    dihval = calc_dih(crd1, crd2, crd3, crd4)

    dis12 = calc_bond(crd1, crd2) #unit is bohr
    dis23 = calc_bond(crd2, crd3) #unit is bohr
    dis34 = calc_bond(crd3, crd4) #unit is bohr

    #get the unit vector
    vec12 = array(crd2) - array(crd1) #vec12 is vec2 - vec1
    vec21 = - vec12
    vec23 = array(crd3) - array(crd2)
    vec32 = - vec23
    vec34 = array(crd4) - array(crd3)
    vec43 = - vec34

    vec12 = array([i/dis12 for i in vec12])
    vec21 = array([i/dis12 for i in vec21])
    vec23 = array([i/dis23 for i in vec23])
    vec32 = array([i/dis23 for i in vec32])
    vec34 = array([i/dis34 for i in vec34])
    vec43 = array([i/dis34 for i in vec43])

    #get the normalized vector
    vecUNABCp = cross(vec32, vec12)
    vecUNABC = array([i/norm(vecUNABCp) for i in vecUNABCp])

    vecUNBCDp = cross(vec43, vec23)
    vecUNBCD = array([i/norm(vecUNBCDp) for i in vecUNBCDp])

    afcmatrix12 = array([[float(0) for x in range(3)] for x in range(3)])
    afcmatrix43 = array([[float(0) for x in range(3)] for x in range(3)])
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix12[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat2-1)+j]
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix43[i][j] = -fcmatrix[3*(nat4-1)+i][3*(nat3-1)+j]
    eigval12, eigvector12 = eig(afcmatrix12)
    eigval43, eigvector43 = eig(afcmatrix43)

    contri12 = 0.0
    contri34 = 0.0
    for i in range(0, 3):
      ev12 = eigvector12[:,i]
      ev43 = eigvector43[:,i]
      contri12 = contri12 + eigval12[i] * abs(dot(vecUNABC, ev12))
      contri34 = contri34 + eigval43[i] * abs(dot(vecUNBCD, ev43))

    contri12 = contri12 * (norm(cross(vec12, vec23)) ** 2)
    contri34 = contri34 * (norm(cross(vec23, vec34)) ** 2)
    contri12 = 1.0 / (contri12 * dis12 * dis12)
    contri34 = 1.0 / (contri34 * dis34 * dis34)

    fcfinal1 = (1.0 / (contri12 + contri34)) * H_TO_KCAL_MOL * 0.5 #/ (math.degrees(1.0)**2)

    fcfinal = fcfinal1 * (math.radians(180.0/float(n2))) **2 / (float(n1) * float(n2))

    fcfinal = fcfinal * scalef * scalef

    return dihval, fcfinal1, fcfinal

def get_imp_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, nat4, scalef):

    crd1 = crds[3*nat1-3:3*nat1]
    crd2 = crds[3*nat2-3:3*nat2]
    crd3 = crds[3*nat3-3:3*nat3]
    crd4 = crds[3*nat4-3:3*nat4]

    dis12 = calc_bond(crd1, crd2) #unit is bohr
    dis13 = calc_bond(crd1, crd3) #unit is bohr
    dis23 = calc_bond(crd2, crd3) #unit is bohr
    dis34 = calc_bond(crd3, crd4) #unit is bohr

    #get the unit vector
    vec12 = array(crd2) - array(crd1) #vec12 is vec2 - vec1
    vec21 = - vec12
    vec23 = array(crd3) - array(crd2)
    vec32 = - vec23
    vec34 = array(crd4) - array(crd3)
    vec43 = - vec34

    vec12 = array([i/dis12 for i in vec12])
    vec21 = array([i/dis12 for i in vec21])
    vec23 = array([i/dis23 for i in vec23])
    vec32 = array([i/dis23 for i in vec32])
    vec34 = array([i/dis34 for i in vec34])
    vec43 = array([i/dis34 for i in vec43])

    #get the normalized vector
    vecUNp = cross(vec43, vec23)
    vecUN = array([i/norm(vecUNp) for i in vecUNp]) #vecUN is the vector
                                     #perpendicular to the plance of ABC

    afcmatrix12 = array([[float(0) for x in range(3)] for x in range(3)])
    afcmatrix13 = array([[float(0) for x in range(3)] for x in range(3)])
    afcmatrix14 = array([[float(0) for x in range(3)] for x in range(3)])
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix12[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat2-1)+j]
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix13[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat3-1)+j]
    for i in range(0, 3):
      for j in range(0, 3):
        afcmatrix14[i][j] = -fcmatrix[3*(nat1-1)+i][3*(nat4-1)+j]

    eigval12, eigvector12 = eig(afcmatrix12)
    eigval13, eigvector13 = eig(afcmatrix13)
    eigval14, eigvector14 = eig(afcmatrix14)
    contri12 = 0.0
    contri13 = 0.0
    contri14 = 0.0
    for i in range(0, 3):
      ev12 = eigvector12[:,i]
      ev13 = eigvector13[:,i]
      ev14 = eigvector14[:,i]
      contri12 = contri12 + eigval12[i] * abs(dot(vecUN, ev12))
      contri13 = contri13 + eigval13[i] * abs(dot(vecUN, ev13))
      contri14 = contri14 + eigval14[i] * abs(dot(vecUN, ev14))

    kAN = (contri12 + contri13 + contri14) * H_TO_KCAL_MOL * 0.5

    fcfinal1 = kAN/(B_TO_A**2)

    #Get hABCD
    pval = (dis12 + dis23 + dis13)/2.0
    sqABC = math.sqrt(pval*(pval-dis12)*(pval-dis23)*(pval-dis13))
    disAtoBC = sqABC * 2.0 / dis23 #distance between A and BC side, unit is bohr

    #get the normalized vector
    vecUNABCp = cross(vec32, vec12)
    vecUNABC = array([i/norm(vecUNABCp) for i in vecUNABCp])
    vecUNBCDp = cross(vec43, vec23)
    vecUNBCD = array([i/norm(vecUNBCDp) for i in vecUNBCDp])

    #dihang = math.acos(dot(vecUNABC, vecUNBCD))

    #hABCD = disAtoBC * dot(vecUNABC, vecUNBCD)

    fcfinal = (disAtoBC ** 2) * kAN / 2.0 #H_TO_KCAL_MOL #HB2_TO_KCAL_MOL_A2

    fcfinal1 = fcfinal1 * scalef * scalef
    fcfinal = fcfinal * scalef * scalef

    return fcfinal1, fcfinal

def gene_by_QM_fitting_sem(scpdbf, ionids, stfpf, pref, finf, chkfname,
                           logfile, g0x, scalef, bondavg, angavg):

    print "==================Using the Seminario method to solve the problem."

    mol, atids, resids = get_atominfo_fpdb(scpdbf)
    blist = get_mc_blist(mol, atids, ionids, stfpf)
    alist = get_alist(mol, blist)   

    #crds after optimization
    if g0x in ['g03', 'g09']:
      crds = get_crds_from_fchk(chkfname, len(atids))
    elif g0x == 'gms':
      crds = get_crds_from_gms(logfile)

    #Whole Hessian Matrix
    if g0x in ['g03', 'g09']:
      fcmatrix = get_matrix_from_fchk(chkfname, 3*len(atids))
    elif g0x == 'gms':
      fcmatrix = get_matrix_from_gms(logfile, 3*len(atids))

    natids = {}
    for i in range(0, len(atids)):
      natids[atids[i]] = i + 1

    attypdict = get_attypdict(stfpf, atids)
    missbondtyps, missangtyps = get_misstyps(pref)

    finalparmdict = {}

    #for bond
    print "=======================Generate the bond parameters==============="
    for misbond in missbondtyps:
      bondlen = []      
      bfconst = []
      for bond in blist:
        at1 = bond[0]
        at2 = bond[1]
        bondtyp = (attypdict[at1], attypdict[at2])
        "The unit in fchk file is a.u. so the distance is in Bohr."
        if bondtyp == misbond or bondtyp[::-1] == misbond:
          nat1 = natids[at1]
          nat2 = natids[at2]

          if bondavg == 1:
            dis, fcfinal, stdv = get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2, scalef, bondavg)
            print '### Bond force constant between ' + \
              mol.atoms[at1].resname + str(mol.atoms[at1].resid) + '@' + mol.atoms[at1].atname + ' and ' + \
              mol.atoms[at2].resname + str(mol.atoms[at2].resid) + '@' + mol.atoms[at2].atname + ' : ' + \
              str(round(fcfinal, 1)) + '±' + str(round(stdv, 1))
          elif bondavg == 0:
            dis, fcfinal = get_bond_fc_with_sem(crds, fcmatrix, nat1, nat2, scalef, bondavg)

          bondlen.append(dis)
          bfconst.append(fcfinal)

      #Get average bond parameters
      misbondat12 = misbond[0] + '-' + misbond[1]
      bond_para = avg_bond_para(misbondat12, bondlen, bfconst)
      #get the force constant
      finalparmdict[misbondat12] = bond_para

    #for angle
    print "=======================Generate the angle parameters=============="
    for misang in missangtyps:
      angvals = []
      afconst = []
      for ang in alist:
        at1 = ang[0]
        at2 = ang[1]
        at3 = ang[2]
        angtyp = (attypdict[at1], attypdict[at2], attypdict[at3])

        if angtyp == misang or angtyp[::-1] == misang:
          nat1 = natids[at1]
          nat2 = natids[at2]
          nat3 = natids[at3]

          if angavg == 1:
            angval, fcfinal, stdv = get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, scalef, angavg)
            print '### Angle force constant between ' + \
              mol.atoms[at1].resname + str(mol.atoms[at1].resid) +  '@' + mol.atoms[at1].atname + ', ' + \
              mol.atoms[at2].resname + str(mol.atoms[at2].resid) +  '@' + mol.atoms[at2].atname + ' and ' + \
              mol.atoms[at3].resname + str(mol.atoms[at3].resid) +  '@' + mol.atoms[at3].atname + ' : ' + \
              str(round(fcfinal, 2)) + '±' + str(round(stdv, 2))
          elif angavg == 0:
             angval, fcfinal = get_ang_fc_with_sem(crds, fcmatrix, nat1, nat2, nat3, scalef, angavg)

          angvals.append(angval)
          afconst.append(fcfinal)

      #Get average angle parameters
      misangat123 = misang[0] + '-' + misang[1] + '-' + misang[2]
      ang_para = avg_angle_para(misangat123, angvals, afconst)
      finalparmdict[misangat123] = ang_para

    #Print out the final frcmod file
    print_frcmod_file(pref, finf, finalparmdict, 'Seminario')

#-----------------------------------------------------------------------------
##Z-matrix method: obtain the force constant from the entire Hessian matrix
#-----------------------------------------------------------------------------
def gene_by_QM_fitting_zmatrix(scpdbf, ionids, stfpf, pref, finf, logfname,
                               scalef):

    print "=============Using the Z-matrix method to generate the parameters."

    sturefs, vals, fcs = get_fc_from_log(logfname)

    #pdb file
    mol, atids, resids = get_atominfo_fpdb(scpdbf)
    blist = get_mc_blist(mol, atids, ionids, stfpf)
    alist = get_alist(mol, blist)   

    #reverse new id dict
    rvnatids = {}
    for i in range(0, len(atids)):
      rvnatids[i+1] = atids[i]

    attypdict = get_attypdict(stfpf, atids)
    missbondtyps, missangtyps = get_misstyps(pref)

    #final parameter dicts
    finalparmdict = {}

    #for bond
    print "=======================Generate the bond parameters==============="
    for misbond in missbondtyps:
      bondlen = []      
      bfconst = []
      for i in range(0, len(sturefs)):
        if len(sturefs[i]) == 2:
          at1 = sturefs[i][0]
          at2 = sturefs[i][1]
          bondtyp = (attypdict[rvnatids[at1]], attypdict[rvnatids[at2]])
          "The unit in log file is Angs."
          if bondtyp == misbond or bondtyp[::-1] == misbond:
            dis = vals[i]
            fcfinal = fcs[i] * HB2_TO_KCAL_MOL_A2 * 0.5
            fcfinal = fcfinal * scalef * scalef
            #Hatree/(Bohr^2) to kcal/(mol*angstrom^2)
            #Times 0.5 factor since AMBER use k(r-r0)^2 but not 1/2*k*(r-r0)^2
            bondlen.append(dis)
            bfconst.append(fcfinal)

      #Get average bond parameters
      misbondat12 = misbond[0] + '-' + misbond[1]
      bond_para = avg_bond_para(misbondat12, bondlen, bfconst)
      #get the force constant
      finalparmdict[misbondat12] = bond_para

    #for angle
    print "=======================Generate the angle parameters=============="
    for misang in missangtyps:
      angvals = []
      afconst = []
      for i in range(0, len(sturefs)):
        if len(sturefs[i]) == 3:
          at1 = sturefs[i][0]
          at2 = sturefs[i][1]
          at3 = sturefs[i][2]
          angtyp = (attypdict[rvnatids[at1]], attypdict[rvnatids[at2]], attypdict[rvnatids[at3]])
          if angtyp == misang or angtyp[::-1] == misang:
            angval = vals[i]
            fcfinal = fcs[i] * H_TO_KCAL_MOL * 0.5
            fcfinal = fcfinal * scalef * scalef
            #Hatree to kcal/mol
            #Times 0.5 factor since AMBER use k(r-r0)^2 but not 1/2*k*(r-r0)^2
            angvals.append(angval)
            afconst.append(fcfinal)

      #Get average angle parameters
      misangat123 = misang[0] + '-' + misang[1] + '-' + misang[2]
      ang_para = avg_angle_para(misangat123, angvals, afconst)
      finalparmdict[misangat123] = ang_para

    #Print out the final frcmod file
    print_frcmod_file(pref, finf, finalparmdict, 'Z-matrix')

