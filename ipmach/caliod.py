#------------------------------------------------------------------------------
# This module is for IOD and CN calculation
#------------------------------------------------------------------------------
import numpy
import os

def cal_iod(rs, grs, maxind):
    rlbind = maxind - 10
    rubind = maxind + 11
    xl = rs[rlbind:rubind]
    yl = grs[rlbind:rubind]
    a, b, c = numpy.polyfit(xl, yl, 2)
    iod = -b/(2*a)
    return iod

def MD_simulation(exe, md_prmtop, md_inpcrd, md_min_steps, md_nvt_steps, md_npt_steps, md_md_steps, ifc4):

    print "Perform MD calculation..."

    #Normal MD simulation-IOD, CN
    #MIN
    print "Perform md_min %d steps..." %md_min_steps
    md_minf = open('md_min.in', 'w')
    print >> md_minf, "rdf simulation first step min 1000"
    print >> md_minf, "&cntrl"
    print >> md_minf, " imin=1, maxcyc=%d, ncyc=1000," %md_min_steps
    print >> md_minf, " ntb=1, cut=10.0, iwrap=1, ioutfm=1,"
    if ifc4 == 1:
      print >> md_minf, " lj1264=1,"
    print >> md_minf, "/"
    md_minf.close()
    os.system("%s -O -i md_min.in -o md_min.out -p %s -c %s -r md_min.rst -x md_min.netcdf" %(exe, md_prmtop, md_inpcrd))

    #NVT
    print "Perform md_nvt %d steps..." %md_nvt_steps
    md_nvtf = open('md_nvt.in', 'w')
    print >> md_nvtf, "rdf simulation second step NVT"
    print >> md_nvtf, " &cntrl"
    print >> md_nvtf, "  imin=0, irest=0, ntx=1, ig=-1,"
    print >> md_nvtf, "  nstlim=%d, dt=0.001," %md_nvt_steps
    print >> md_nvtf, "  cut=10.0, ntb=1,"
    print >> md_nvtf, "  ntc=2, ntf=2, "
    print >> md_nvtf, "  ntt=3, gamma_ln=5.0, tempi=0.0, temp0=300.0,"
    print >> md_nvtf, "  ntpr=500, ntwx=2000, ntwr=2000, ioutfm=1, ntwv=-1, iwrap=1,"
    if ifc4 == 1:
      print >> md_nvtf, " lj1264=1,"
    print >> md_nvtf, "/"
    md_nvtf.close()
    os.system("%s -O -i md_nvt.in -o md_nvt.out -p %s -c md_min.rst -r md_nvt.rst -x md_nvt.netcdf" %(exe, md_prmtop))

    #NPT
    print "Perform md_npt %d steps..." %md_npt_steps
    md_nptf = open('md_npt.in', 'w')
    print >> md_nptf, "for deltaG calculation"
    print >> md_nptf, " &cntrl"
    print >> md_nptf, "  imin=0, irest=1, ntx=5, ig=-1,"
    print >> md_nptf, "  nstlim=%d, dt=0.001," %md_npt_steps
    print >> md_nptf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
    print >> md_nptf, "  ntc=2, ntf=2, "
    print >> md_nptf, "  tempi=300.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
    print >> md_nptf, "  ntpr=500, ntwr=2000, ntwx=2000, ioutfm=1, iwrap=1, ntwv=-1,"
    if ifc4 == 1:
      print >> md_nptf, " lj1264=1,"
    print >> md_nptf, "/"
    md_nptf.close()
    os.system("%s -O -i md_npt.in -o md_npt.out -p %s -c md_nvt.rst -r md_npt.rst -x md_npt.netcdf" %(exe, md_prmtop))

    #MD
    print "Perform md_md %d steps..." %md_md_steps
    md_mdf = open('md_md.in', 'w')
    print >> md_mdf, "for deltaG calculation"
    print >> md_mdf, " &cntrl"
    print >> md_mdf, "  imin=0, irest=1, ntx=5, ig=-1,"
    print >> md_mdf, "  nstlim=%d, dt=0.001," %md_md_steps
    print >> md_mdf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
    print >> md_mdf, "  ntc=2, ntf=2, "
    print >> md_mdf, "  tempi=300.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
    print >> md_mdf, "  ntpr=500, ntwr=500, ntwx=500, ioutfm=1, iwrap=1, ntwv=-1,"
    if ifc4 == 1:
      print >> md_mdf, " lj1264=1,"
    print >> md_mdf, "/"
    md_mdf.close()
    os.system("%s -O -i md_md.in -o md_md.out -p %s -c md_npt.rst -r md_md.rst -x md_md.netcdf" %(exe, md_prmtop))

    #Cpptraj input file
    cpptrajf = open('cpptraj.in', 'w')
    print >> cpptrajf, "trajin md_md.netcdf 1 100000 1"
    print >> cpptrajf, "radial M_O 0.01 5.0 :1 :WAT@O volume intrdf M_O"
    cpptrajf.close()

    os.system("cpptraj -p %s -i cpptraj.in > cpptraj.out" %md_prmtop)
    #Get the IOD
    rdff = open('M_O', 'r')
    dr = 0.005
    rs = []
    grs = []
    ccns = []
    for line in rdff:
      if line[0] != "#":
        r, gr, ccn = line.split()
        r = float(r)
        gr = float(gr)
        ccn = float(ccn)
        rs.append(r)
        grs.append(gr)
        ccns.append(ccn)
    rdff.close()

    #Find the first peak
    maxind = grs.index(max(grs))
    iod = cal_iod(rs, grs, maxind)
    drs = [abs(i-iod) for i in rs]

    #Find the first peak again
    maxind2 = drs.index(min(drs))
    iod = cal_iod(rs, grs, maxind2)
    iod = round(iod, 2)

    #Find the second peak
    grs2 = grs[maxind+100:]
    maxind3 = grs2.index(max(grs2))
    grs3 = grs[maxind:maxind+maxind3]

    #Find the first minimum
    minind = grs3.index(min(grs3))
    cn = ccns[maxind+minind]
    cn = round(cn, 1)

    return iod, cn

