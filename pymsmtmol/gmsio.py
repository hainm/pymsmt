"This module is for GAMESS"
import linecache

def get_esp_from_log2(logfile, espfile):

    #------------Coordinate List for the Atom and ESP Center--------------
    #Log file uses Angstrom, esp file uses Bohr

    B_TO_A = 0.529177249 #Bohr to Angstrom

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

def get_esp_from_log(logfile, espfile):

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

get_esp_from_log('test2.log', 'test2.esp')
get_esp_from_log2('1OKL_large_mk.log', '1OKL_large_mk.esp')

