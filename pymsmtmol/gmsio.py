"This module is for GAMESS"
import linecache
import numpy
from chemistry.periodic_table import AtomicNum
from constants import B_TO_A

#------------------------------------------------------------------------------
#--------------------------Write GAMESS input file-----------------------------
#------------------------------------------------------------------------------

def write_gms_optf(goptf2, totchg, SpinNum, gatms, signum=3):

    ##GAMESS OPT file
    optf2 = open(goptf2, 'w')
    print >> optf2, " $SYSTEM MEMDDI=400 MWORDS=200 $END"
    print >> optf2, " $CONTRL DFTTYP=B3LYP RUNTYP=OPTIMIZE ICHARG=%d MULT=%d $END" %(int(round(totchg, 0)), SpinNum)
    print >> optf2, " $STATPT NSTEP=1000 $END"
    print >> optf2, " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END"
    print >> optf2, " $DATA"
    print >> optf2, "Cluster/6-31G"
    print >> optf2, "C1"
    optf2.close()

    for gatmi in gatms:
      write_gmsatm(gatmi, goptf2, signum)

    ##Print the last line in GAMESS input file
    ##Geometry Optimization file
    optf2 = open(goptf2, 'a')
    print >> optf2, " $END"
    optf2.close()

def write_gms_fcf(gfcf2, totchg, SpinNum):

    ##GAMESS FC file
    fcf2 = open(gfcf2, 'w')
    print >> fcf2, " $SYSTEM MEMDDI=400 MWORDS=200 $END"
    print >> fcf2, " $CONTRL DFTTYP=B3LYP RUNTYP=HESSIAN ICHARG=%d MULT=%d $END" %(int(round(totchg, 0)), SpinNum)
    print >> fcf2, " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END"
    print >> fcf2, " $DATA"
    print >> fcf2, "Cluster/6-31G"
    print >> fcf2, "C1"
    print >> fcf2, " "
    print >> fcf2, " $END"
    fcf2.close()

def write_gms_mkf(gmsf, totchg, SpinNum, gatms, signum=3):

    #For GAMESS MK Charge file
    w_gmsf = open(gmsf, 'w')
    print >> w_gmsf, " $SYSTEM MEMDDI=400 MWORDS=200 $END"
    print >> w_gmsf, " $CONTRL DFTTYP=B3LYP ICHARG=%d MULT=%d $END" %(int(round(totchg, 0)), SpinNum)
    print >> w_gmsf, " $ELPOT IEPOT=1 WHERE=PDC $END"
    print >> w_gmsf, " $PDC PTSEL=CONNOLLY CONSTR=NONE $END"
    print >> w_gmsf, " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END"
    print >> w_gmsf, " $DATA"
    print >> w_gmsf, "Cluster/6-31G(d)"
    print >> w_gmsf, "C1"
    w_gmsf.close()

    #For GAMESS file
    for gatmi in gatms:
      write_gmsatm(gatmi, gmsf, signum)

    #Print the end character for GAMESS input file
    w_gmsf = open(gmsf, 'a')
    print >> w_gmsf, ' $END'
    w_gmsf.close()

def write_gmsatm(gmsatm, fname, signum=3):

    wf = open(fname, 'a')
    element = gmsatm.element
    nuchg = AtomicNum[element]
    nuchg = round(nuchg, 1)
    if signum == 3:
        print >> wf, "%-6s  %6.1f  %8.3f %8.3f %8.3f" %(gmsatm.element, \
                 nuchg, gmsatm.crdx, gmsatm.crdy, gmsatm.crdz)
    elif signum == 4:
        print >> wf, "%-6s  %6.1f  %9.4f %9.4f %9.4f" %(gmsatm.element, \
                 nuchg, gmsatm.crdx, gmsatm.crdy, gmsatm.crdz)
    wf.close()

#------------------------------------------------------------------------------
#----------------------Read info from GAMESS output file-----------------------
#------------------------------------------------------------------------------

def get_crds_from_gms(logfile):

    unit = 'angs' #Coordinates will use angs unit in default

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if ' ATOM      ATOMIC                      COORDINATES (BOHR)' in line:
            bln = ln + 2
            unit = 'bohr'    #Means using Bohr unit
        elif ' ATOM      ATOMIC                      COORDINATES (ANGS' in line:
            bln = ln + 2
            unit = 'angs'
        elif '          INTERNUCLEAR DISTANCES' in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    crdl = []
    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        if unit == 'bohr':
            crdl.append(float(line[2]))
            crdl.append(float(line[3]))
            crdl.append(float(line[4]))
        elif unit == 'angs':
            crdl.append(float(line[2])/B_TO_A)
            crdl.append(float(line[3])/B_TO_A)
            crdl.append(float(line[4])/B_TO_A)
    linecache.clearcache()

    return crdl

def get_matrix_from_gms(logfile, msize):

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'CARTESIAN FORCE CONSTANT MATRIX' in line:
            bln = ln + 6
        ln = ln + 1
    fp.close()

    fcmatrix = numpy.array([[float(0) for x in range(msize)] for x in range(msize)])

    for i in range(0, msize/6): #To see how many cycles need
        for j in range(0, msize-i*6): #How many lines in the cycle
            line = linecache.getline(logfile, bln)
            line = line.strip('\n')
            for k in range(0, 6): #There are 6 values in one line
                fcmatrix[j+i*6][k+i*6] = float(line[20+9*k:29+9*k])
            bln = bln + 1
        bln = bln + 4

    #If there is one more section remaining
    if (msize%6 == 3):
        for j in range(0, 3):
            line = linecache.getline(logfile, bln)
            line = line.strip('\n')
            for k in range(0, 3): #There are 6 values in one line
                fcmatrix[j + msize/6 * 6][k + msize/6 * 6] = float(line[20+9*k:29+9*k])
            bln = bln + 1

    linecache.clearcache()

    #To complete the matrix
    for j in range(0, msize):
        for k in range(0, j+1):
            fcmatrix[k][j] = fcmatrix[j][k]

    return fcmatrix

def get_esp_from_gms(logfile, espfile):

    #ESP file uses Bohr unit

    unit = 'bohr' #In default ESP coordinates will use bohr unit

    #---------For Coordinates--------
    crdl = []

    bln1 = 0
    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'COORDINATES OF ALL ATOMS ARE (ANGS)' in line:
            bln1 = ln + 3
            unit = 'angs'
        elif 'COORDINATES OF ALL ATOMS ARE (BOHR)' in line:
            bln1 = ln + 3
            unit = 'bohr'
        elif ' ATOM      ATOMIC                      COORDINATES (BOHR)' in line:
            bln2 = ln + 2
            unit = 'bohr'
        elif ' ATOM      ATOMIC                      COORDINATES (ANGS' in line:
            bln2 = ln + 2
            unit = 'angs'
        elif 'INTERNUCLEAR DISTANCES' in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    if bln1 == 0:
        bln = bln2
    else:
        bln = bln1

    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        if unit == 'angs':
            crd = (float(line[2])/B_TO_A, float(line[3])/B_TO_A, float(line[4])/B_TO_A)
        elif unit == 'bohr':
            crd = (float(line[2]), float(line[3]), float(line[4]))
        crdl.append(crd)

    linecache.clearcache()

    #---------For ESP points---------
    #ESP files use Bohr and Atomic Unit Charge

    espdict = {}

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'ELECTROSTATIC POTENTIAL' in line:
            bln = ln + 6
        ln = ln + 1
    fp.close()

    try:
        line = linecache.getline(logfile, bln)
        line = line.strip('\n')
        line = line.split()
        numps = int(line[-1])
    except:
        line = linecache.getline(logfile, bln)
        print "CAUTION: " + line.strip('\n') + " IN THE GAMESS CALCULATION."
        bln = bln + 1
        line = linecache.getline(logfile, bln)
        line = line.strip('\n')
        line = line.split()
        numps = int(line[-1])

    for i in range(bln+1, bln+1+numps):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        val = (float(line[1]), float(line[2]), float(line[3]), float(line[6]))
        espdict[int(line[0])] = val

    linecache.clearcache()

    w_espf = open(espfile, 'w')
    print >> w_espf, "%5d%5d%5d" %(len(crdl), len(espdict), 0)
    for i in crdl:
        print >> w_espf, "%16s %15.7E %15.7E %15.7E" %(' ', i[0], i[1], i[2])
    for i in range(1, len(espdict)+1):
        val = espdict[i]
        print >> w_espf, "%16.7E %15.7E %15.7E %15.7E" %(val[3], val[0], val[1], val[2])
    w_espf.close()


