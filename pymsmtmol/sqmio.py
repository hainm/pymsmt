"This module for SQM"
import linecache
from pymsmtmol.mol import gauatm
from chemistry.periodic_table import AtomicNum

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

def write_sqm_optf(siopf, totchg, gatms):

    sqm_scf = open(siopf, 'w')
    print >> sqm_scf, "Run semi-empirical minimization"
    print >> sqm_scf, " &qmmm"
    print >> sqm_scf, " qm_theory='PM6', grms_tol=0.0002,"
    print >> sqm_scf, " tight_p_conv=1, scfconv=1.d-10, qmcharge=%d," %int(totchg)
    print >> sqm_scf, " /"
    for gatmi in gatms:
        nuchg = int(AtomicNum[gatmi.element])
        print >> sqm_scf, "%-2s %5s %10.4f %10.4f %10.4f" \
        %(nuchg, gatmi.element, gatmi.crdx, gatmi.crdy, gatmi.crdz)
    sqm_scf.close()


#def write_sqm_mkf(simkf, totchg, gatms):
#    sqm_lgf = open(simkf, 'w')
#    print >> sqm_lgf, "Run semi-empirical minimization"
#    print >> sqm_lgf, " &qmmm"
#    print >> sqm_lgf, " qm_theory='PM6', grms_tol=0.0002,"
#    print >> sqm_lgf, " tight_p_conv=1, scfconv=1.d-10, qmcharge=%d," %int(totchg)
#    print >> sqm_lgf, " /"
#    sqm_lgf.close()







