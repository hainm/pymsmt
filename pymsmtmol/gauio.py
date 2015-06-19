"""
Module for writting a Gaussian file and read the coordinates and force
constants from Gaussian output file.
"""
import numpy
from constants import B_TO_A
from chemistry.periodic_table import AtomicNum

#------------------------------------------------------------------------------
#------------------------Write Gaussian input file-----------------------------
#------------------------------------------------------------------------------

def write_sdd_basis(gatms, gauf):

    atnames = [i.element for i in gatms]
    atnames = list(set(atnames))
    atnums = [AtomicNum[i] for i in atnames]

    w_gauf = open(gauf, 'a')
    print >> w_gauf, " "
    for i in atnames:
      if AtomicNum[i] <= 18:
        print >> w_gauf, i,
    print >> w_gauf, "0"
    print >> w_gauf, "6-31G*"
    print >> w_gauf, "****"

    for i in atnames:
      if Atomic_Num[i] >= 19:
        print >> w_gauf, i
    print >> w_gauf, "0"
    print >> w_gauf, "SDD"
    print >> w_gauf, "****"
    print >> w_gauf, " "

    for i in atnames:
      if Atomic_Num[i] >= 19:
        print >> w_gauf, i
    print >> w_gauf, "0"
    print >> w_gauf, "SDD"
    print >> w_gauf, "****"
    w_gauf.close()

def write_gau_optf(outf, goptf, totchg, SpinNum, gatms, signum=3):

    ##Geometry Optimization file
    optf = open(goptf, 'w')
    print >> optf, "$RunGauss"
    print >> optf, "%%Chk=%s_sidechain_opt.chk" %outf
    print >> optf, "%Mem=3000MB"
    print >> optf, "%NProcShared=2"
    print >> optf, "#N B3LYP/6-31G* Geom=PrintInputOrient " + \
                   "Integral=(Grid=UltraFine) Opt"
    print >> optf, "SCF=XQC"
    print >> optf, " "
    print >> optf, "CLR"
    print >> optf, " "
    print >> optf, "%d  %d" %(int(round(totchg, 0)), SpinNum)
    optf.close()

    if signum == 3:
      for gatmi in gatms:
        write_gauatm(gatmi, goptf)
    elif signum == 4:
      for gatmi in gatms:
        write_gauatm(gatmi, goptf, 4)

    ##Geometry Optimization file
    optf = open(goptf, 'a')
    print >> optf, " "
    print >> optf, " "
    optf.close()

def write_gau_fcf(outf, gfcf):

    ##Force constant calculation file
    fcf = open(gfcf, 'w')
    print >> fcf, "$RunGauss"
    print >> fcf, "%%Chk=%s_sidechain_opt.chk" %outf
    print >> fcf, "%Mem=3000MB"
    print >> fcf, "%NProcShared=2"
    print >> fcf, "#N B3LYP/6-31G* Freq=NoRaman Geom=AllCheckpoint Guess=Read"
    print >> fcf, "Integral=(Grid=UltraFine) SCF=XQC IOp(7/33=1)"
    print >> fcf, " "
    print >> fcf, " "
    fcf.close()

def write_gau_mkf(outf, gmkf, totchg, SpinNum, gatms, ionnames, chargedict, IonLJParaDict, signum=3):

    ##MK RESP input file
    mkf = open(gmkf, 'w')
    print >> mkf, "$RunGauss"
    print >> mkf, "%%Chk=%s_large_mk.chk" %outf
    print >> mkf, "%Mem=3000MB"
    print >> mkf, "%NProcShared=2"
    print >> mkf, "#N B3LYP/6-31G* Integral=(Grid=UltraFine) Opt " + \
                  "Pop(MK,ReadRadii) SCF=XQC"
    print >> mkf, "IOp(6/33=2)"
    print >> mkf, " "
    print >> mkf, "CLR"
    print >> mkf, " "
    print >> mkf, "%d  %d" %(int(round(totchg, 0)),SpinNum)
    mkf.close()

    #For Gaussian file
    if signum == 3:
      for gatmi in gatms:
        write_gauatm_opth(gatmi, gmkf)
    elif signum == 4:
      for gatmi in gatms:
        write_gauatm(gatmi, gmkf, signum)

    ##print the ion radius for resp charge fitting in MK RESP input file
    mkf = open(gmkf, 'a')
    print >> mkf, " "
    for i in ionnames:
      chg = str(int(chargedict[i]))
      if len(i) > 1:
        i = i[0] + i[1:].lower()
      vdwradius = IonLJParaDict[i + chg][0]
      print >> mkf, i, vdwradius
    print >> mkf, " "
    for i in ionnames:
      chg = str(int(chargedict[i]))
      if len(i) > 1:
        i = i[0] + i[1:].lower()
      vdwradius = IonLJParaDict[i + chg][0]
      print >> mkf, i, vdwradius
    print >> mkf, " "
    print >> mkf, " "
    mkf.close()

def write_gauatm(gauatm, fname, signum=3):
    wf = open(fname, 'a')
    if signum == 3:
      print >> wf, "%-6s   %8.3f %8.3f %8.3f" %(gauatm.element, \
                   gauatm.crdx, gauatm.crdy, gauatm.crdz)
    elif signum == 4:
      print >> wf, "%-6s   %9.4f %9.4f %9.4f" %(gauatm.element, \
                   gauatm.crdx, gauatm.crdy, gauatm.crdz)
    wf.close()

def write_gauatm_opth(gauatm, fname, signum=3):
    wf = open(fname, 'a')
    if gauatm.element == "H":
      print >> wf, "%-6s  0 %8.3f %8.3f %8.3f" %(gauatm.element, \
          gauatm.crdx, gauatm.crdy, gauatm.crdz)
    else:
      print >> wf, "%-6s -1 %8.3f %8.3f %8.3f" %(gauatm.element, \
          gauatm.crdx, gauatm.crdy, gauatm.crdz)
    wf.close()

#------------------------------------------------------------------------------
#-----------------------Read info from Gaussian output file--------------------
#------------------------------------------------------------------------------

def get_crds_from_fchk(fname, bstring, estring):

    #fchk file uses Bohr unit

    crds = []

    bl = len(bstring)
    blc = 0

    el = len(estring)
    elc = 0

    fp = open(fname, 'r')
    i = 1
    for line in fp:
      if (line[0:bl] == bstring) and (blc == 0):
        beginl = i
        blc = blc + 1
      elif (line[0:el] == estring) and (elc == 0):
        endl = i
        elc = elc + 1
      i = i + 1
    fp.close()

    fp1 = open(fname, 'r')
    i = 1
    for line in fp1:
      if (i > beginl) and (i < endl):
        crd = line.split()
        for j in crd:
          if j != ' ':
            crds.append(float(j))
      i = i + 1
    fp1.close()

    return crds

def get_matrix_from_fchk(fname, bstring, estring, msize):

    crds = []

    bl = len(bstring)
    blc = 0

    el = len(estring)
    elc = 0

    fp = open(fname, 'r')
    i = 1
    for line in fp:
      if (line[0:bl] == bstring) and (blc == 0):
        beginl = i
        blc = blc + 1
      elif (line[0:el] == estring) and (elc == 0):
        endl = i
        elc = elc + 1
      i = i + 1
    fp.close()

    fcmatrix = numpy.array([[float(0) for x in range(msize)] for x in range(msize)])

    fp1 = open(fname, 'r')
    i = 1
    for line in fp1:
      if (i > beginl) and (i < endl):
        crd = line.split()
        for j in crd:
          if j != ' ':
            crds.append(float(j))
      i = i + 1
    fp1.close()

    i = 0
    for j in range(0, msize):
      for k in range(0, j+1):
        fcmatrix[j][k] = crds[i]
        fcmatrix[k][j] = crds[i]
        i = i + 1

    return fcmatrix

def get_crds_from_log(logfname, g0x):

    #Log file uses angs. as unit

    if g0x == 'g03':
      fp = open(logfname)
      ln = 1
      for line in fp:
        if 'Redundant internal coordinates' in line:
          bln = ln + 3
        elif 'Recover connectivity data from disk' in line:
          eln = ln - 1
        ln = ln + 1
      fp.close()
    elif g0x == 'g09':
      fp = open(logfname)
      ln = 1
      for line in fp:
        if 'Redundant internal coordinates' in line:
          bln = ln + 1
        elif 'Recover connectivity data from disk' in line:
          eln = ln - 1
        ln = ln + 1
      fp.close()

    crds = []
    ln = 1
    fp1 = open(logfname)
    for line in fp1:
      if (ln >= bln) and (ln <= eln):
        line = line.strip('\n')
        line = line.split(',')
        line = line[-3:]
        line = [float(i) for i in line]
        crds += line
      ln = ln + 1
    fp1.close()

    return crds

def get_fc_from_log(logfname):

    stringle = 'calculate D2E/DX2 analytically'
    stringfc = ' Internal force constants:'

    sturefs = []
    vals = []

    ##Get the values for each bond, angle and dihedral
    fp = open(logfname, 'r')
    for line in fp:
      if stringle in line:
        line = line.strip('\n')
        line = line.strip('!')
        line = line.lstrip(' !') 
        line = line.split()
        val = float(line[2])
        vals.append(val)
        typ = line[0][0]
        line[1] = line[1].strip(typ)
        line[1] = line[1].lstrip('(')
        line[1] = line[1].strip(')')
        ats = line[1].split(',')
        ats = [int(i) for i in ats]
        sturefs.append(tuple(ats))
    fp.close()

    maxnum = len(sturefs)

    ##Get the force constant begin line
    fp = open(logfname, 'r')
    lnum = 1
    for line in fp:
      if stringfc in line:
        blnum = lnum + 2
      lnum = lnum + 1
    fp.close()

    blnum1 = blnum

    ##Get the line number list to read the force constant
    numl = []

    if (maxnum%5 != 0):
      for i in range(0, int(maxnum/5)+1):
        gap = maxnum - len(numl) + 1
        for j in range(0, 5):
          if len(numl) < maxnum:
            numl.append(blnum1 + j)
        blnum1 = blnum1 + gap
    else:
      for i in range(0, int(maxnum/5)):
        gap = maxnum - len(numl) + 1
        for j in range(0, 5):
          if len(numl) < maxnum:
            numl.append(blnum1 + j)
        blnum1 = blnum1 + gap

    fcs = []
    fp = open(logfname, 'r')
    lnum = 1
    for line in fp:
      if lnum in numl:
        line = line.strip('\n')
        line = line.split()
        numv = len(line)
        fc = float(line[-1].replace('D', 'e'))
        fcs.append(fc)
      lnum = lnum + 1
    fp.close()

    ##Return three lists: identifications, values, force constants
    return sturefs, vals, fcs

def get_esp_from_gau(logfile, espfile):
    #Gaussian log file uses Angstrom as unit, esp file uses Bohr

    #------------Coordinate List for the Atom and ESP Center--------------
    crdl1 = []
    crdl2 = []

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'Electrostatic Properties Using The SCF Density' in line:
            bln = ln
        ln = ln + 1
    fp.close()

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if ln >= bln:
            if '      Atomic Center' in line:
                line = line.strip('\n')
                line = line.split()
                crd = (float(line[-3])/B_TO_A, float(line[-2])/B_TO_A, float(line[-1])/B_TO_A)
                crdl1.append(crd)
            elif ('     ESP Fit Center' in line):
                line = line.strip('\n')
                line = line.split()
                crd = (float(line[-3])/B_TO_A, float(line[-2])/B_TO_A, float(line[-1])/B_TO_A)
                crdl2.append(crd)
        ln = ln + 1
    fp.close()

    #------------ESP values for Atom and ESP Center--------------------
    #Both log and esp files use Atomic Unit Charge

    espl1 = []
    espl2 = []

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'Electrostatic Properties (Atomic Units)' in line:
            bln = ln + 6
        ln = ln + 1
    fp.close()

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if ln >= bln:
            if ' Atom' in line:
                line = line.strip('\n')
                line = line.split()
                esp = float(line[-1])
                espl1.append(esp)
            elif (' Fit ' in line):
                line = line.strip('\n')
                line = line.split()
                esp = float(line[-1])
                espl2.append(esp)
        ln = ln + 1
    fp.close()

    #----------------Check and print-----------------------
    if (len(crdl1) == len(espl1)) and (len(crdl2) == len(espl2)):
        w_espf = open(espfile, 'w')
        print >> w_espf, "%5d%5d%5d" %(len(crdl1), len(crdl2), 0)
        for i in range(0, len(crdl1)):
            crd = crdl1[i]
            print >> w_espf, "%16s %15.7E %15.7E %15.7E" %(' ', crd[0], crd[1], crd[2])
        for i in range(0, len(crdl2)):
            crd = crdl2[i]
            esp = espl2[i]
            print >> w_espf, "%16.7E %15.7E %15.7E %15.7E" %(esp, crd[0], crd[1], crd[2])
        w_espf.close()
    else:
        raise ValueError("The length of coordinates and ESP charges are different!")

