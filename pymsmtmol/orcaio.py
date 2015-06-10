"This module is for ORCA"

def write_orca_optf(outf, ofcf, totchg, SpinNum, gatms):

    ##ORCA constant calculation file
    orcafcf = open(ofcf, 'w')
    print >> orcafcf, "%pal"
    print >> orcafcf, "  nprocs 2"
    print >> orcafcf, "end"
    print >> orcafcf, "%MaxCore 1500"
    print >> orcafcf, "! B3LYP 6-31G* TightSCF NRSCF"
    print >> orcafcf, "%method"
    print >> orcafcf, "  Grid 6"
    print >> orcafcf, "end"
    print >> orcafcf, "! Opt"
    print >> orcafcf, "%%base \"%s\"" %outf
    print >> orcafcf, '*xyz', str(int(totchg)), str(SpinNum)
    orcafcf.close()

    #Coordinates
    for gatmi in gatms:
      write_gauatm(gatmi, ofcf)

    #END
    orcafcf = open(ofcf, 'a')
    print >> orcafcf, "*"
    orcafcf.close()


