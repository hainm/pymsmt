"This module for SQM"
import linecache
from pymsmtmol.mol import gauatm

def get_crdinfo_from_sqm(outfile):

    gauatms = []

    ln = 1  
    fp = open(outfile, 'r')
    for line in fp:
        if " Final Structure" in line:
            bln = ln + 4
        elif "--------- Calculation Completed ----------" in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    for i in range(bln, eln+1):
        line = linecache.getline(outfile, i)
        line = line.strip('\n')
        line = line.split()
        if line[0] == 'QMMM:':
            atm = gauatm(line[3], float(line[4]), float(line[5]), float(line[6]))
            gauatms.append(atm)

    linecache.clearcache()

    return gauatms

