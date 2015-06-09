"This module is for GAMESS"
import linecache
import numpy
from chemistry.periodic_table import AtomicNum

def write_gmsatm(gmsatm, fname):

    wf = open(fname, 'a')
    element = gmsatm.element
    nuchg = AtomicNum[element]
    nuchg = round(nuchg, 1)

    print >> wf, "%-6s  %6.1f  %8.3f%8.3f%8.3f" %(gmsatm.element, \
        nuchg, gmsatm.crdx, gmsatm.crdy, gmsatm.crdz)
    wf.close()

def get_crds_from_gms(logfile):

    B_TO_A = 0.529177249 #Bohr to Angstrom

    unit = 'angs' #Means using Angstrom unit

    fp = open(logfile, 'r')
    for line in fp:
        if ' ATOM      ATOMIC                      COORDINATES (BOHR)' in line:
            bln = ln + 2
            unit = 'au'    #means using Bohr unit
        elif ' ATOM      ATOMIC                      COORDINATES (ANGS.)' in line:
            bln = ln + 2
            unit = 'angs'
        elif '          INTERNUCLEAR DISTANCES' in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        if unit == 'bohr':
            crd = (float(line[2])*B_TO_A, float(line[3])*B_TO_A, float(line[4])*B_TO_A)
        elif unit == 'angs':
            crd = (float(line[2]), float(line[3]), float(line[4]))
        crdl.append(crd)

    linecache.clearcache()

    return crdl

def get_esp_from_gms(logfile, espfile):

    B_TO_A = 0.529177249 #Bohr to Angstrom

    #---------For Coordinates--------
    #Log file use Angstrom, esp file uses Bohr

    crdl = []

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'COORDINATES OF ALL ATOMS ARE' in line:
            bln = ln + 3
        elif 'INTERNUCLEAR DISTANCES' in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        crd = (float(line[2])/B_TO_A, float(line[3])/B_TO_A, float(line[4])/B_TO_A)
        crdl.append(crd)

    linecache.clearcache()

    #---------For ESP points---------
    #Both log and esp files use Bohr and Atomic Unit Charge

    espdict = {}

    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'ELECTROSTATIC POTENTIAL' in line:
            bln = ln + 6
        ln = ln + 1
    fp.close()

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
                print j+i*6, k+i*6, fcmatrix[j+i*6][k+i*6]
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


