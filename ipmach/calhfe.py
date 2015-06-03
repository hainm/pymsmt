#------------------------------------------------------------------------------
# For TI calculation
#------------------------------------------------------------------------------
import os

def get_lamadas(windows):
    if windows == 12:
      lamadas = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
                 0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]
    elif windows == 9:
      lamadas = [0.01592, 0.08198, 0.19331, 0.33787, 0.50000, 0.66213,
                 0.80669, 0.91802, 0.98408]
    elif windows == 7:
      lamadas = [0.02544, 0.12923, 0.29707, 0.50000, 0.70292, 0.87076,
                 0.97455]
    elif windows == 5:
      lamadas = [0.04691, 0.23076, 0.50000, 0.76923, 0.95308]
    elif windows == 3:
      lamadas = [0.11270, 0.50000, 0.88729] 
    elif windows == 2:
      lamadas = [0.21132, 0.78867]
    elif windows == 1:
      lamadas = [0.50000]
    return lamadas

def get_weights(windows):
    if windows == 12:
      weights = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457, 0.12457,
                 0.11675, 0.10158, 0.08004, 0.05347, 0.02359]
    elif windows == 9:
      weights = [0.04064, 0.09032, 0.13031, 0.15617, 0.16512, 0.15617, 0.13031,
                 0.09032, 0.04064]
    elif windows == 7:
      weights = [0.06474, 0.13985, 0.19091, 0.20897, 0.19091, 0.13985, 0.06474]
    elif windows == 5:
      weights = [0.11846, 0.23931, 0.28444, 0.23931, 0.11846]
    elif windows == 3:
      weights = [0.27777, 0.44444, 0.27777]
    elif windows == 2:
      weights = [0.5, 0.5]
    elif windows == 1:
      weights = [1.0]
    return weights

#------------------------------------------------------------------------------
#-----One Step TI for pmemd
#------------------------------------------------------------------------------

def OneStep_pTI(ti_windows, ti_window_steps, ti_sample_steps, exe,
                ti_prmtop, ti_inpcrd, ti_min_steps, ti_nvt_steps,
                ti_npt_steps, rev):

    print "Perform TI calculation..."

    #Get lamadas
    lamadas = get_lamadas(ti_windows)

    #1. Total appearing--------------------------------------------------------
    for i in range(1, ti_windows+1):
      lamada = lamadas[i-1]

      tif = open(str(i)+'_fwd.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=1, ntx=5, ig=-1,"
      print >> tif, "  nstlim=%d, dt=0.001," %ti_window_steps
      print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
      print >> tif, "  timask1 = \":1\", scmask1 = \":1\","
      print >> tif, "  timask2 = \":2\", scmask2 = \":2\","
      print >> tif, "/"
      tif.close()

      if i == 1:
        #1.1 Minimization
        print "Perform ti_min %d steps" %ti_min_steps
        ti_minf = open('ti_min.in', 'w')
        print >> ti_minf, "for deltaG calculation"
        print >> ti_minf, " &cntrl"
        print >> ti_minf, "  imin=1, irest=0, ntx=1, ntmin=2,"
        print >> ti_minf, "  maxcyc=%d," %ti_min_steps
        print >> ti_minf, "  ntb=1, cut=10.0,"
        print >> ti_minf, "  ntc=1, ntf=1,"
        print >> ti_minf, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1,"
        print >> ti_minf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_minf, "  timask1 = \":1\", scmask1 = \":1\","
        print >> ti_minf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_minf, "/"
        ti_minf.close()

        os.system("%s -O -i ti_min.in -o ti_min.out -p %s -c %s -r ti_min.rst -x ti_min.netcdf" %(exe, ti_prmtop, ti_inpcrd))

        #1.2 NVT heating
        print "Perform nvt %d steps..." %ti_nvt_steps
        ti_nvtf = open('ti_nvt.in', 'w')
        print >> ti_nvtf, "for deltaG calculation"
        print >> ti_nvtf, " &cntrl"
        print >> ti_nvtf, "  imin=0, irest=0, ntx=1,"
        print >> ti_nvtf, "  nstlim=%d, dt=0.001," %ti_nvt_steps
        print >> ti_nvtf, "  ntb=1, cut=10.0,"
        print >> ti_nvtf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nvtf, "  ntc=1, ntf=1,"
        print >> ti_nvtf, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nvtf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_nvtf, "  timask1 = \":1\", scmask1 = \":1\","
        print >> ti_nvtf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_nvtf, "/"
        ti_nvtf.close()

        os.system("%s -O -i ti_nvt.in -o ti_nvt.out -p %s -c ti_min.rst -r ti_nvt.rst -x ti_nvt.netcdf" %(exe, ti_prmtop))

        #1.3 NPT equil
        print "Perform npt %d steps..." %ti_npt_steps
        ti_nptf = open('ti_npt.in', 'w')
        print >> ti_nptf, "for deltaG calculation"
        print >> ti_nptf, " &cntrl"
        print >> ti_nptf, "  imin=0, irest=0, ntx=1,"
        print >> ti_nptf, "  nstlim=%d, dt=0.001," %ti_npt_steps
        print >> ti_nptf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
        print >> ti_nptf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nptf, "  ntc=1, ntf=1,"
        print >> ti_nptf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nptf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_nptf, "  timask1 = \":1\", scmask1 = \":1\","
        print >> ti_nptf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_nptf, "/"
        ti_nptf.close()

        os.system("%s -O -i ti_npt.in -o ti_npt.out -p %s -c ti_nvt.rst -r ti_npt.rst -x ti_npt.netcdf" %(exe, ti_prmtop))

        #1.4 TI
        print "TI step" + str(i) + " %d steps (forward)..." %ti_window_steps
        os.system("%s -O -i %d_fwd.in -o %d_fwd.out -p %s -c ti_npt.rst -r %d_fwd.rst -x %d_fwd.netcdf" %(exe, i, i, ti_prmtop, i, i))

      else:
        j = i - 1
        print "TI step" + str(i) + " %d steps (forward)..." %ti_window_steps
        os.system("%s -O -i %d_fwd.in -o %d_fwd.out -p %s -c %d_fwd.rst -r %d_fwd.rst -x %d_fwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #2. Total disappearing---------------------------------------------------
    if rev == 1:
      for i in range(1, ti_windows+1):

        lamada = lamadas[i-1]

        tif = open(str(i)+'_bwd.in', 'w')
        print >> tif, "for deltaG calculation"
        print >> tif, " &cntrl"
        print >> tif, "  imin=0, irest=0, ntx=1,"
        print >> tif, "  nstlim=%d, dt=0.001," %ti_window_steps
        print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
        print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> tif, "  ntc=1, ntf=1,"
        print >> tif, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                      "iwrap=1, ntwv=-1,"
        print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> tif, "  timask1 = \":2\", scmask1 = \":2\","
        print >> tif, "  timask2 = \":1\", scmask2 = \":1\","
        print >> tif, "/"
        tif.close()

        if i == 1:
          print "TI step" + str(i) + " %d steps (backward)..." %ti_window_steps
          os.system("%s -O -i %d_bwd.in -o %d_bwd.out -p %s -c %d_fwd.rst -r %d_bwd.rst -x %d_bwd.netcdf" %(exe, i, i, ti_prmtop, ti_windows, i, i))
        else:
          print "TI step" + str(i) + " %d steps (backward)..." %ti_window_steps
          j = i - 1
          os.system("%s -O -i %d_bwd.in -o %d_bwd.out -p %s -c %d_bwd.rst -r %d_bwd.rst -x %d_bwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    weights = get_weights(ti_windows)

    dG1 = 0.0
    for i in range(1, ti_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_fwd.out', "DV\/DL", ti_sample_steps/500)).read()
      dvdl = float(dvdl) * weights[i-1]
      dG1 += dvdl
    dG1 = round(dG1, 2)

    if rev == 1:
      dG2 = 0.0
      for i in range(1, ti_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_bwd.out', "DV\/DL", ti_sample_steps/500)).read()
        dvdl = float(dvdl) * weights[i-1]
        dG2 += dvdl
      dG2 = round(dG2, 2)
      dG = (dG1 - dG2)/2
      dG = round(dG, 2)
    else:
      dG = dG1

    return dG

#------------------------------------------------------------------------------
#-------One step TI for sander
#------------------------------------------------------------------------------

def OneStep_sTI(ti_windows, ti_window_steps, ti_sample_steps, minexe, exe,
                md_prmtop, md_inpcrd, md0_prmtop, md0_inpcrd,
                ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4):

    print "Perform TI calculation..."

    #Get lamadas
    lamadas = get_lamadas(ti_windows)

    #1. Total appearing--------------------------------------------------------
    for i in range(1, ti_windows+1):
      lamada = lamadas[i-1]

      #Input file
      tif = open(str(i)+'_v1.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=1, ntx=5, ig=-1,"
      print >> tif, "  nstlim=%d, dt=0.001," %ti_window_steps
      print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
      print >> tif, "  scmask = \":1\","
      if ifc4 == 1:
        print >> tif, "  lj1264=1,"
      print >> tif, "/"
      tif.close()
      os.system("cp %d_v1.in %d_v2.in" %(i, i))
      
      if i == 1:

        #1.1 Minimization
        print "Perform ti_min %d steps" %ti_min_steps
        ti_minf = open('ti_min1.in', 'w')
        print >> ti_minf, "for deltaG calculation"
        print >> ti_minf, " &cntrl"
        print >> ti_minf, "  imin=1, irest=0, ntx=1, ntmin=2,"
        print >> ti_minf, "  maxcyc=%d," %ti_min_steps
        print >> ti_minf, "  ntb=1, cut=10.0,"
        print >> ti_minf, "  ntc=1, ntf=1,"
        print >> ti_minf, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1,"
        print >> ti_minf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_minf, "  scmask = \":1\","
        if ifc4 == 1:
          print >> ti_minf, "  lj1264=1,"
        print >> ti_minf, "/"
        ti_minf.close()

        #Group file
        os.system("cp ti_min1.in ti_min2.in")
        groupf = open('ti_min.group', 'w')
        print >> groupf, "-O -i ti_min1.in -o ti_min1.out -p %s -c %s -inf ti_min1.info -x ti_min1.netcdf -r ti_min1.rst" %(md0_prmtop, md0_inpcrd)
        print >> groupf, "-O -i ti_min2.in -o ti_min2.out -p %s -c %s -inf ti_min2.info -x ti_min2.netcdf -r ti_min2.rst" %(md_prmtop, md_inpcrd)
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_min.group" %minexe)

        #1.2 NVT heating
        print "Perform nvt %d steps..." %ti_nvt_steps

        ti_nvtf = open('ti_nvt1.in', 'w')
        print >> ti_nvtf, "for deltaG calculation"
        print >> ti_nvtf, " &cntrl"
        print >> ti_nvtf, "  imin=0, irest=0, ntx=1,"
        print >> ti_nvtf, "  nstlim=%d, dt=0.001," %ti_nvt_steps
        print >> ti_nvtf, "  ntb=1, cut=10.0,"
        print >> ti_nvtf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nvtf, "  ntc=1, ntf=1,"
        print >> ti_nvtf, "  ntpr=500, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nvtf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_nvtf, "  scmask = \":1\","
        if ifc4 == 1:
          print >> ti_nvtf, "  lj1264=1,"
        print >> ti_nvtf, "/"
        ti_nvtf.close()

        #Group file
        os.system("cp ti_nvt1.in ti_nvt2.in")
        groupf = open('ti_nvt.group', 'w')
        print >> groupf, "-O -i ti_nvt1.in -o ti_nvt1.out -p %s -c ti_min1.rst -inf ti_nvt1.info -x ti_nvt1.netcdf -r ti_nvt1.rst" %md0_prmtop
        print >> groupf, "-O -i ti_nvt2.in -o ti_nvt2.out -p %s -c ti_min2.rst -inf ti_nvt2.info -x ti_nvt2.netcdf -r ti_nvt2.rst" %md_prmtop
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_nvt.group" %exe)

        #1.3 NPT equil
        print "Perform npt %d steps..." %ti_npt_steps
        ti_nptf = open('ti_npt1.in', 'w')
        print >> ti_nptf, "for deltaG calculation"
        print >> ti_nptf, " &cntrl"
        print >> ti_nptf, "  imin=0, irest=0, ntx=1,"
        print >> ti_nptf, "  nstlim=%d, dt=0.001," %ti_npt_steps
        print >> ti_nptf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
        print >> ti_nptf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nptf, "  ntc=1, ntf=1,"
        print >> ti_nptf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nptf, "  icfe=1, ifsc=1, clambda=%6.5f," %lamada
        print >> ti_nptf, "  scmask = \":1\","
        if ifc4 == 1:
          print >> ti_nptf, "  lj1264=1,"
        print >> ti_nptf, "/"
        ti_nptf.close()

        #Group file
        os.system("cp ti_npt1.in ti_npt2.in")
        groupf = open('ti_npt.group', 'w')
        print >> groupf, "-O -i ti_npt1.in -o ti_npt1.out -p %s -c ti_nvt1.rst -inf ti_npt1.info -x ti_npt1.netcdf -r ti_npt1.rst" %md0_prmtop
        print >> groupf, "-O -i ti_npt2.in -o ti_npt2.out -p %s -c ti_nvt2.rst -inf ti_npt2.info -x ti_npt2.netcdf -r ti_npt2.rst" %md_prmtop
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_npt.group" %exe)

        #1.4 TI
        print "TI step" + str(i) + " %d steps (forward)..." %ti_window_steps
        #Group file
        groupf = open(str(i)+'_fwd.group', 'w')
        print >> groupf, "-O -i %d_v1.in -o %d_fwd1.out -p %s -c ti_npt1.rst -inf %d_fwd1.info -x %d_fwd1.netcdf -r %d_fwd1.rst" %(i, i, md0_prmtop, i, i, i)
        print >> groupf, "-O -i %d_v2.in -o %d_fwd2.out -p %s -c ti_npt2.rst -inf %d_fwd2.info -x %d_fwd2.netcdf -r %d_fwd2.rst" %(i, i, md_prmtop, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_fwd.group" %(exe, i))

      else:
        j = i - 1
        print "TI step" + str(i) + " %d steps (forward)..." %ti_window_steps
        #Group file
        groupf = open(str(i)+'_fwd.group', 'w')
        print >> groupf, "-O -i %d_v1.in -o %d_fwd1.out -p %s -c %d_fwd1.rst -inf %d_fwd1.info -x %d_fwd1.netcdf -r %d_fwd1.rst" %(i, i, md0_prmtop, j, i, i, i)
        print >> groupf, "-O -i %d_v2.in -o %d_fwd2.out -p %s -c %d_fwd2.rst -inf %d_fwd2.info -x %d_fwd2.netcdf -r %d_fwd2.rst" %(i, i, md_prmtop, j, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_fwd.group" %(exe, i))

    #2. Total disappearing---------------------------------------------------
    if rev == 1:
      for i in range(1, ti_windows+1):
        if i == 1:
          print "TI step" + str(i) + " %d steps (backward)..." %ti_window_steps
          #Group file
          groupf = open(str(i)+'_bwd.group', 'w')
          print >> groupf, "-O -i %d_v2.in -o %d_bwd2.out -p %s -c %d_fwd2.rst -inf %d_bwd2.info -x %d_bwd2.netcdf -r %d_bwd2.rst" %(i, i, md_prmtop, ti_windows, i, i, i)
          print >> groupf, "-O -i %d_v1.in -o %d_bwd1.out -p %s -c %d_fwd1.rst -inf %d_bwd1.info -x %d_bwd1.netcdf -r %d_bwd1.rst" %(i, i, md0_prmtop, ti_windows, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_bwd.group" %(exe, i))
        else:
          print "TI step" + str(i) + " %d steps (backward)..." %ti_window_steps
          j = i - 1
          #Group file
          groupf = open(str(i)+'_bwd.group', 'w')
          print >> groupf, "-O -i %d_v2.in -o %d_bwd2.out -p %s -c %d_bwd2.rst -inf %d_bwd2.info -x %d_bwd2.netcdf -r %d_bwd2.rst" %(i, i, md_prmtop, j, i, i, i)
          print >> groupf, "-O -i %d_v1.in -o %d_bwd1.out -p %s -c %d_bwd1.rst -inf %d_bwd1.info -x %d_bwd1.netcdf -r %d_bwd1.rst" %(i, i, md0_prmtop, j, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_bwd.group" %(exe, i))

    weights = get_weights(ti_windows)

    dG1 = 0.0
    for i in range(1, ti_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_fwd1.out', "DV\/DL", ti_sample_steps/500)).read()
      dvdl = float(dvdl) * weights[i-1]
      dG1 += dvdl
    dG1 = round(dG1, 2)

    if rev == 1:
      dG2 = 0.0
      for i in range(1, ti_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_bwd1.out', "DV\/DL", ti_sample_steps/500)).read()
        dvdl = float(dvdl) * weights[i-1]
        dG2 += dvdl
      dG2 = round(dG2, 2)
      dG = (dG1 - dG2)/2
      dG = round(dG, 2)
    else:
      dG = dG1

    return dG

#------------------------------------------------------------------------------
#---------Two Step TI for pmemd
#------------------------------------------------------------------------------
def TwoStep_pTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
               chg_window_steps, vdw_sample_steps, chg_sample_steps,
               exe, vdw_prmtop, ti_prmtop, ti_inpcrd, ti_min_steps,
               ti_nvt_steps, ti_npt_steps, rev):

    print "Perform TI calculation..."

    #Get lamadas
    vdw_lamadas = get_lamadas(ti_vdw_windows)
    chg_lamadas = get_lamadas(ti_chg_windows)

    #1. VDW appearing-------------------------------------------------------
    for i in range(1, ti_vdw_windows+1):

      vdw_lamada = vdw_lamadas[i-1]

      tif = open(str(i)+'_vdw_fwd.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=0, ntx=1,"
      print >> tif, "  nstlim=%d, dt=0.001," %vdw_window_steps
      print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
      print >> tif, "  timask1 = \":1\", scmask1 = \":1\", crgmask = \":1|:2\","
      print >> tif, "  timask2 = \":2\", scmask2 = \":2\","
      print >> tif, "/"
      tif.close()

      if i == 1:
        #1.1 Minimization
        print "Perform min %d steps..." %ti_min_steps
        ti_minf = open('ti_min.in', 'w')
        print >> ti_minf, "for deltaG calculation"
        print >> ti_minf, " &cntrl"
        print >> ti_minf, "  imin=1, ntx=1, ntmin=2,"
        print >> ti_minf, "  maxcyc=%d," %ti_min_steps
        print >> ti_minf, "  ntb=1, cut=10.0,"
        print >> ti_minf, "  ntc=1, ntf=1,"
        print >> ti_minf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1,"
        print >> ti_minf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_minf, "  timask1 = \":1\", scmask1 = \":1\", crgmask = \":1|:2\","
        print >> ti_minf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_minf, "/"
        ti_minf.close()

        os.system("%s -O -i ti_min.in -o ti_min.out -p %s -c %s -r ti_min.rst -x ti_min.netcdf" %(exe, ti_prmtop, ti_inpcrd))

        #1.2 NVT heating
        print "Perform nvt %d steps..." %ti_nvt_steps
        ti_nvtf = open('ti_nvt.in', 'w')
        print >> ti_nvtf, "for deltaG calculation"
        print >> ti_nvtf, " &cntrl"
        print >> ti_nvtf, "  imin=0, irest=0, ntx=1, ig=-1,"
        print >> ti_nvtf, "  nstlim=%d, dt=0.001," %ti_nvt_steps
        print >> ti_nvtf, "  ntb=1, cut=10.0,"
        print >> ti_nvtf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nvtf, "  ntc=1, ntf=1,"
        print >> ti_nvtf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nvtf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_nvtf, "  timask1 = \":1\", scmask1 = \":1\", crgmask = \":1|:2\","
        print >> ti_nvtf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_nvtf, "/"
        ti_nvtf.close()

        os.system("%s -O -i ti_nvt.in -o ti_nvt.out -p %s -c ti_min.rst -r ti_nvt.rst -x ti_nvt.netcdf" %(exe, ti_prmtop))

        #1.3 NPT equil
        print "Perform npt %d steps..." %ti_npt_steps
        ti_nptf = open('ti_npt.in', 'w')
        print >> ti_nptf, "for deltaG calculation"
        print >> ti_nptf, " &cntrl"
        print >> ti_nptf, "  imin=0, irest=1, ntx=7, ig=-1,"
        print >> ti_nptf, "  nstlim=%d, dt=0.001," %ti_npt_steps
        print >> ti_nptf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
        print >> ti_nptf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nptf, "  ntc=1, ntf=1,"
        print >> ti_nptf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nptf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_nptf, "  timask1 = \":1\", scmask1 = \":1\", crgmask = \":1|:2\","
        print >> ti_nptf, "  timask2 = \":2\", scmask2 = \":2\","
        print >> ti_nptf, "/"
        ti_nptf.close()

        os.system("%s -O -i ti_npt.in -o ti_npt.out -p %s -c ti_nvt.rst -r ti_npt.rst -x ti_npt.netcdf" %(exe, ti_prmtop))

        #1.4 TI
        print "VDW TI step" + str(i) + " %d steps (forward)..." %vdw_window_steps
        os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c ti_npt.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe, i, i, ti_prmtop, i, i))
      else:
        j = i - 1
        print "VDW TI step" + str(i) + " %d steps (forward)..." %vdw_window_steps
        os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c %d_vdw_fwd.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #One more step to equ the structure
    tif = open(str(ti_vdw_windows+1)+'_vdw_fwd.in', 'w')
    print >> tif, "for deltaG calculation"
    print >> tif, " &cntrl"
    print >> tif, "  imin=0, irest=0, ntx=1,"
    print >> tif, "  nstlim=%d, dt=0.001," %vdw_window_steps
    print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
    print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
    print >> tif, "  ntc=1, ntf=1,"
    print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                  "iwrap=1, ntwv=-1,"
    print >> tif, "  icfe=1, ifsc=1, clambda=0.98000,"
    print >> tif, "  timask1 = \":1\", scmask1 = \":1\", crgmask = \":1|:2\","
    print >> tif, "  timask2 = \":2\", scmask2 = \":2\","
    print >> tif, "/"
    tif.close()

    print "VDW TI step" + str(ti_vdw_windows+1) + " to equ structure %d steps (forward)..." %vdw_window_steps
    os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c %d_vdw_fwd.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe,
               ti_vdw_windows+1, ti_vdw_windows+1, ti_prmtop, ti_vdw_windows, ti_vdw_windows+1, ti_vdw_windows+1))

    #transfer the rst file
    os.system("awk 'NR<=2' %d_vdw_fwd.rst > %d_vdw_fwd_merge.rst" %(ti_vdw_windows+1, ti_vdw_windows+1))
    os.system("awk 'NR==3' %d_vdw_fwd.rst | awk '{printf( \"%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f\\n\", $1, $2, $3, $1, $2, $3)}' >> %d_vdw_fwd_merge.rst" %(ti_vdw_windows+1, ti_vdw_windows+1))
    os.system("awk 'NR>3' %d_vdw_fwd.rst >> %d_vdw_fwd_merge.rst" %(ti_vdw_windows+1, ti_vdw_windows+1))

    #2. Charge appearing-------------------------------------------------------
    for i in range(1, ti_chg_windows+1):

      chg_lamada = chg_lamadas[i-1]

      tif = open(str(i)+'_chg_fwd.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=1, ntx=7, ig=-1,"
      print >> tif, "  nstlim=%d, dt=0.001," %chg_window_steps
      print >> tif, "  ntb=2, cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, clambda=%6.5f," %chg_lamada
      print >> tif, "  timask1 = \":1\","
      print >> tif, "  timask2 = \":2\","
      print >> tif, "/"
      tif.close()

      if i == 1:
        print "Charge TI step" + str(i) + " %d steps (forward)..." %chg_window_steps
        os.system("%s -O -i %d_chg_fwd.in -o %d_chg_fwd.out -p %s -c %d_vdw_fwd_merge.rst -r %d_chg_fwd.rst -x %d_chg_fwd.netcdf" %(exe, i, i, vdw_prmtop, ti_vdw_windows+1, i, i))
      else:
        print "Charge TI step" + str(i) + " %d steps (forward)..." %chg_window_steps
        j = i - 1
        os.system("%s -O -i %d_chg_fwd.in -o %d_chg_fwd.out -p %s -c %d_chg_fwd.rst -r %d_chg_fwd.rst -x %d_chg_fwd.netcdf" %(exe, i, i, vdw_prmtop, j, i, i))

    if rev == 1:
      #3. Charge disappearing----------------------------------------------------
      for i in range(1, ti_chg_windows+1):
        chg_lamada = chg_lamadas[i-1]
        tif = open(str(i)+'_chg_bwd.in', 'w')
        print >> tif, "for deltaG calculation"
        print >> tif, " &cntrl"
        print >> tif, "  imin=0, irest=1, ntx=7, ig=-1,"
        print >> tif, "  nstlim=%d, dt=0.001," %chg_window_steps
        print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
        print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> tif, "  ntc=1, ntf=1,"
        print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                      "iwrap=1, ntwv=-1,"
        print >> tif, "  icfe=1, clambda=%6.5f," %chg_lamada
        print >> tif, "  timask1 = \":2\","
        print >> tif, "  timask2 = \":1\","
        print >> tif, "/"
        tif.close()

        if i == 1:
          print "Charge TI step" + str(i) + " %d steps (backward)..." %chg_window_steps
          os.system("%s -O -i %d_chg_bwd.in -o %d_chg_bwd.out -p %s -c %d_chg_fwd.rst -r %d_chg_bwd.rst -x %d_chg_bwd.netcdf" %(exe, i, i, vdw_prmtop, ti_chg_windows, i, i))
        else:
          print "Charge TI step" + str(i) + " %d steps (backward)..." %chg_window_steps
          j = i - 1
          os.system("%s -O -i %d_chg_bwd.in -o %d_chg_bwd.out -p %s -c %d_chg_bwd.rst -r %d_chg_bwd.rst -x %d_chg_bwd.netcdf" %(exe, i, i, vdw_prmtop, j, i, i))

      #4. VDW disappearing----------------------------------------------------
      for i in range(1, ti_vdw_windows+1):

        vdw_lamada = vdw_lamadas[i-1]

        tif = open(str(i)+'_vdw_bwd.in', 'w')
        print >> tif, "for deltaG calculation"
        print >> tif, " &cntrl"
        print >> tif, "  imin=0, irest=1, ntx=7, ig=-1,"
        print >> tif, "  nstlim=%d, dt=0.001," %vdw_window_steps
        print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=2,"
        print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> tif, "  ntc=1, ntf=1,"
        print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                      "iwrap=1, ntwv=-1,"
        print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> tif, "  timask1 = \":2\", scmask1 = \":2\", crgmask = \":1|:2\","
        print >> tif, "  timask2 = \":1\", scmask2 = \":1\","
        print >> tif, "/"
        tif.close()

        if i == 1:
          print "VDW TI step" + str(i) + " %d steps (backward)..." %vdw_window_steps
          os.system("%s -O -i %d_vdw_bwd.in -o %d_vdw_bwd.out -p %s -c %d_chg_bwd.rst -r %d_vdw_bwd.rst -x %d_vdw_bwd.netcdf" %(exe, i, i, ti_prmtop, ti_chg_windows, i, i))
        else:
          print "VDW TI step" + str(i) + " %d steps (backward)..." %vdw_window_steps
          j = i - 1
          os.system("%s -O -i %d_vdw_bwd.in -o %d_vdw_bwd.out -p %s -c %d_vdw_bwd.rst -r %d_vdw_bwd.rst -x %d_vdw_bwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #Calcualte the free energies
    vdw_weights = get_weights(ti_vdw_windows)
    chg_weights = get_weights(ti_chg_windows)

    dG1 = 0.0
    for i in range(1, ti_vdw_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_vdw_fwd.out', "DV\/DL", vdw_sample_steps/500)).read()
      dvdl = float(dvdl) * vdw_weights[i-1]
      dG1 += dvdl
    dG1 = round(dG1, 2)

    dG2 = 0.0
    for i in range(1, ti_chg_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_chg_fwd.out', "DV\/DL", chg_sample_steps/500)).read()
      dvdl = float(dvdl) * chg_weights[i-1]
      dG2 += dvdl
    dG2 = round(dG2, 2)

    if rev == 1:
      dG3 = 0.0
      for i in range(1, ti_chg_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_chg_bwd.out', "DV\/DL", chg_sample_steps/500)).read()
        dvdl = float(dvdl) * chg_weights[i-1]
        dG3 += dvdl
      dG3 = round(dG3, 2)

      dG4 = 0.0
      for i in range(1, ti_vdw_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_vdw_bwd.out', "DV\/DL", vdw_sample_steps/500)).read()
        dvdl = float(dvdl) * vdw_weights[i-1]
        dG4 += dvdl
      dG4 = round(dG4, 2)

      dG = (dG1 + dG2 - dG3 - dG4)/2
      dG = round(dG, 2)
    else:
      dG = dG1 + dG2
      dG = round(dG, 2)

    return dG

#------------------------------------------------------------------------------
#---------Two Step TI for sander
#------------------------------------------------------------------------------
def TwoStep_sTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
            chg_window_steps, vdw_sample_steps, chg_sample_steps,
            minexe, exe, md0_prmtop, md0_inpcrd, mdv_prmtop, mdv_inpcrd, md_prmtop,
            ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4):

    print "Perform TI calculation..."

    #Get lamadas
    vdw_lamadas = get_lamadas(ti_vdw_windows)
    chg_lamadas = get_lamadas(ti_chg_windows)

    #1. VDW appearing-------------------------------------------------------
    for i in range(1, ti_vdw_windows+1):

      vdw_lamada = vdw_lamadas[i-1]

      tif = open(str(i)+'_vdw1.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=0, ntx=1,"
      print >> tif, "  nstlim=%d, dt=0.001," %vdw_window_steps
      print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
      print >> tif, "  scmask = \":1\", crgmask = \":1\","
      print >> tif, "/"
      tif.close()
      os.system("cp %d_vdw1.in %d_vdw2.in" %(i, i))

      if i == 1:
        #1.1 Minimization
        print "Perform min %d steps..." %ti_min_steps
        ti_minf = open('ti_min1.in', 'w')
        print >> ti_minf, "for deltaG calculation"
        print >> ti_minf, " &cntrl"
        print >> ti_minf, "  imin=1, ntx=1, ntmin=2,"
        print >> ti_minf, "  maxcyc=%d," %ti_min_steps
        print >> ti_minf, "  ntb=1, cut=10.0,"
        print >> ti_minf, "  ntc=1, ntf=1,"
        print >> ti_minf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1,"
        print >> ti_minf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_minf, "  scmask = \":1\", crgmask = \":1\","
        print >> ti_minf, "/"
        ti_minf.close()
        os.system("cp ti_min1.in ti_min2.in")

        #Group file
        groupf = open('ti_min.group', 'w')
        print >> groupf, "-O -i ti_min1.in -o ti_min1.out -p %s -c %s -inf ti_min1.info -x ti_min1.netcdf -r ti_min1.rst" %(md0_prmtop, md0_inpcrd)
        print >> groupf, "-O -i ti_min2.in -o ti_min2.out -p %s -c %s -inf ti_min2.info -x ti_min2.netcdf -r ti_min2.rst" %(mdv_prmtop, mdv_inpcrd)
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_min.group" %minexe)

        #1.2 NVT heating
        print "Perform nvt %d steps..." %ti_nvt_steps
        ti_nvtf = open('ti_nvt1.in', 'w')
        print >> ti_nvtf, "for deltaG calculation"
        print >> ti_nvtf, " &cntrl"
        print >> ti_nvtf, "  imin=0, irest=0, ntx=1, ig=-1,"
        print >> ti_nvtf, "  nstlim=%d, dt=0.001," %ti_nvt_steps
        print >> ti_nvtf, "  ntb=1, cut=10.0,"
        print >> ti_nvtf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nvtf, "  ntc=1, ntf=1,"
        print >> ti_nvtf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nvtf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_nvtf, "  scmask = \":1\", crgmask = \":1\","
        print >> ti_nvtf, "/"
        ti_nvtf.close()
        os.system("cp ti_nvt1.in ti_nvt2.in")

        #Group file
        groupf = open('ti_nvt.group', 'w')
        print >> groupf, "-O -i ti_nvt1.in -o ti_nvt1.out -p %s -c ti_min1.rst -inf ti_nvt1.info -x ti_nvt1.netcdf -r ti_nvt1.rst" %md0_prmtop
        print >> groupf, "-O -i ti_nvt2.in -o ti_nvt2.out -p %s -c ti_min2.rst -inf ti_nvt2.info -x ti_nvt2.netcdf -r ti_nvt2.rst" %mdv_prmtop
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_nvt.group" %exe)

        #1.3 NPT equil
        print "Perform npt %d steps..." %ti_npt_steps
        ti_nptf = open('ti_npt1.in', 'w')
        print >> ti_nptf, "for deltaG calculation"
        print >> ti_nptf, " &cntrl"
        print >> ti_nptf, "  imin=0, irest=1, ntx=7, ig=-1,"
        print >> ti_nptf, "  nstlim=%d, dt=0.001," %ti_npt_steps
        print >> ti_nptf, "  cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
        print >> ti_nptf, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
        print >> ti_nptf, "  ntc=1, ntf=1,"
        print >> ti_nptf, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, " + \
                       "ioutfm=1, iwrap=1, ntwv=-1,"
        print >> ti_nptf, "  icfe=1, ifsc=1, clambda=%6.5f," %vdw_lamada
        print >> ti_nptf, "  scmask = \":1\", crgmask = \":1\","
        print >> ti_nptf, "/"
        ti_nptf.close()
        os.system("cp ti_npt1.in ti_npt2.in")

        #Group file
        groupf = open('ti_npt.group', 'w')
        print >> groupf, "-O -i ti_npt1.in -o ti_npt1.out -p %s -c ti_nvt1.rst -inf ti_npt1.info -x ti_npt1.netcdf -r ti_npt1.rst" %md0_prmtop
        print >> groupf, "-O -i ti_npt2.in -o ti_npt2.out -p %s -c ti_nvt2.rst -inf ti_npt2.info -x ti_npt2.netcdf -r ti_npt2.rst" %mdv_prmtop
        groupf.close()
        os.system("%s -ng 2 -groupfile ti_npt.group" %exe)

        #1.4 TI
        print "VDW TI step" + str(i) + " %d steps (forward)..." %vdw_window_steps
        #Group file
        groupf = open(str(i)+'_vdw_fwd.group', 'w')
        print >> groupf, "-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c ti_npt1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(i, i, md0_prmtop, i, i, i)
        print >> groupf, "-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c ti_npt2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(i, i, mdv_prmtop, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, i))
      else:
        j = i - 1
        print "VDW TI step" + str(i) + " %d steps (forward)..." %vdw_window_steps
        #Group file
        groupf = open(str(i)+'_vdw_fwd.group', 'w')
        print >> groupf, "-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(i, i, md0_prmtop, j, i, i, i)
        print >> groupf, "-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(i, i, mdv_prmtop, j, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, i))

    #One more step to equ the structure
    k = ti_vdw_windows + 1
    tif = open(str(k)+'_vdw1.in', 'w')
    print >> tif, "for deltaG calculation"
    print >> tif, " &cntrl"
    print >> tif, "  imin=0, irest=0, ntx=1,"
    print >> tif, "  nstlim=%d, dt=0.001," %vdw_window_steps
    print >> tif, "  cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
    print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
    print >> tif, "  ntc=1, ntf=1,"
    print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                  "iwrap=1, ntwv=-1,"
    print >> tif, "  icfe=1, ifsc=1, clambda=0.98000,"
    print >> tif, "  scmask = \":1\", crgmask = \":1\","
    print >> tif, "/"
    tif.close()
    os.system("cp %d_vdw1.in %d_vdw2.in" %(k, k))

    print "VDW TI step" + str(k) + " to equ structure %d steps (forward)..." %vdw_window_steps
    j = k - 1
    #Group file
    groupf = open(str(k)+'_vdw_fwd.group', 'w')
    print >> groupf, "-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(k, k, md0_prmtop, j, k, k, k)
    print >> groupf, "-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(k, k, mdv_prmtop, j, k, k, k)
    groupf.close()
    os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, k))

    #2. Charge appearing-------------------------------------------------------
    for i in range(1, ti_chg_windows+1):

      chg_lamada = chg_lamadas[i-1]

      tif = open(str(i)+'_chg1.in', 'w')
      print >> tif, "for deltaG calculation"
      print >> tif, " &cntrl"
      print >> tif, "  imin=0, irest=1, ntx=7, ig=-1,"
      print >> tif, "  nstlim=%d, dt=0.001," %chg_window_steps
      print >> tif, "  ntb=2, cut=10.0, ntp=1, pres0=1.01325, barostat=1,"
      print >> tif, "  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,"
      print >> tif, "  ntc=1, ntf=1,"
      print >> tif, "  ntpr=100, ntwr=500, ntwx=10000, ntave=500, ioutfm=1, " + \
                    "iwrap=1, ntwv=-1,"
      print >> tif, "  icfe=1, clambda=%6.5f," %chg_lamada
      if ifc4 == 1:
        print >> tif, "  lj1264=1,"
      print >> tif, "/"
      tif.close()
      os.system("cp %d_chg1.in %d_chg2.in" %(i, i))

      if i == 1:
        print "Charge TI step" + str(i) + " %d steps (forward)..." %chg_window_steps
        #Group file
        groupf = open(str(i)+'_chg_fwd.group', 'w')
        print >> groupf, "-O -i %d_chg1.in -o %d_chg_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_chg_fwd1.info -x %d_chg_fwd1.netcdf -r %d_chg_fwd1.rst" %(i, i, mdv_prmtop, k, i, i, i)
        print >> groupf, "-O -i %d_chg2.in -o %d_chg_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_chg_fwd2.info -x %d_chg_fwd2.netcdf -r %d_chg_fwd2.rst" %(i, i, md_prmtop, k, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_chg_fwd.group" %(exe, i))
      else:
        print "Charge TI step" + str(i) + " %d steps (forward)..." %chg_window_steps
        j = i - 1
        #Group file
        groupf = open(str(i)+'_chg_fwd.group', 'w')
        print >> groupf, "-O -i %d_chg1.in -o %d_chg_fwd1.out -p %s -c %d_chg_fwd1.rst -inf %d_chg_fwd1.info -x %d_chg_fwd1.netcdf -r %d_chg_fwd1.rst" %(i, i, mdv_prmtop, j, i, i, i)
        print >> groupf, "-O -i %d_chg2.in -o %d_chg_fwd2.out -p %s -c %d_chg_fwd2.rst -inf %d_chg_fwd2.info -x %d_chg_fwd2.netcdf -r %d_chg_fwd2.rst" %(i, i, md_prmtop, j, i, i, i)
        groupf.close()
        os.system("%s -ng 2 -groupfile %d_chg_fwd.group" %(exe, i))

    if rev == 1:
      #3. Charge disappearing----------------------------------------------------
      for i in range(1, ti_chg_windows+1):
        if i == 1:
          print "Charge TI step" + str(i) + " %d steps (backward)..." %chg_window_steps
          #Group file
          groupf = open(str(i)+'_chg_bwd.group', 'w')
          print >> groupf, "-O -i %d_chg2.in -o %d_chg_bwd2.out -p %s -c %d_chg_fwd2.rst -inf %d_chg_bwd2.info -x %d_chg_bwd2.netcdf -r %d_chg_bwd2.rst" %(i, i, md_prmtop, ti_chg_windows, i, i, i)
          print >> groupf, "-O -i %d_chg1.in -o %d_chg_bwd1.out -p %s -c %d_chg_fwd1.rst -inf %d_chg_bwd1.info -x %d_chg_bwd1.netcdf -r %d_chg_bwd1.rst" %(i, i, mdv_prmtop, ti_chg_windows, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_chg_bwd.group" %(exe, i))
        else:
          print "Charge TI step" + str(i) + " %d steps (backward)..." %chg_window_steps
          j = i - 1
          #Group file
          groupf = open(str(i)+'_chg_bwd.group', 'w')
          print >> groupf, "-O -i %d_chg2.in -o %d_chg_bwd2.out -p %s -c %d_chg_bwd2.rst -inf %d_chg_bwd2.info -x %d_chg_bwd2.netcdf -r %d_chg_bwd2.rst" %(i, i, md_prmtop, j, i, i, i)
          print >> groupf, "-O -i %d_chg1.in -o %d_chg_bwd1.out -p %s -c %d_chg_bwd1.rst -inf %d_chg_bwd1.info -x %d_chg_bwd1.netcdf -r %d_chg_bwd1.rst" %(i, i, mdv_prmtop, j, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_chg_bwd.group" %(exe, i))

      #4. VDW disappearing----------------------------------------------------
      for i in range(1, ti_vdw_windows+1):
        if i == 1:
          print "VDW TI step" + str(i) + " %d steps (backward)..." %vdw_window_steps
          groupf = open(str(i)+'_vdw_bwd.group', 'w')
          print >> groupf, "-O -i %d_vdw2.in -o %d_vdw_bwd2.out -p %s -c %d_chg_bwd2.rst -inf %d_vdw_bwd2.info -x %d_vdw_bwd2.netcdf -r %d_vdw_bwd2.rst" %(i, i, mdv_prmtop, ti_chg_windows, i, i, i)
          print >> groupf, "-O -i %d_vdw1.in -o %d_vdw_bwd1.out -p %s -c %d_chg_bwd1.rst -inf %d_vdw_bwd1.info -x %d_vdw_bwd1.netcdf -r %d_vdw_bwd1.rst" %(i, i, md0_prmtop, ti_chg_windows, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_vdw_bwd.group" %(exe, i))
        else:
          print "VDW TI step" + str(i) + " %d steps (backward)..." %vdw_window_steps
          j = i - 1
          #Group file
          groupf = open(str(i)+'_vdw_bwd.group', 'w')
          print >> groupf, "-O -i %d_vdw2.in -o %d_vdw_bwd2.out -p %s -c %d_vdw_bwd2.rst -inf %d_vdw_bwd2.info -x %d_vdw_bwd2.netcdf -r %d_vdw_bwd2.rst" %(i, i, mdv_prmtop, j, i, i, i)
          print >> groupf, "-O -i %d_vdw1.in -o %d_vdw_bwd1.out -p %s -c %d_vdw_bwd1.rst -inf %d_vdw_bwd1.info -x %d_vdw_bwd1.netcdf -r %d_vdw_bwd1.rst" %(i, i, md0_prmtop, j, i, i, i)
          groupf.close()
          os.system("%s -ng 2 -groupfile %d_vdw_bwd.group" %(exe, i))

    #Calcualte the free energies
    vdw_weights = get_weights(ti_vdw_windows)
    chg_weights = get_weights(ti_chg_windows)

    dG1 = 0.0
    for i in range(1, ti_vdw_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_vdw_fwd1.out', "DV\/DL", vdw_sample_steps/500)).read()
      dvdl = float(dvdl) * vdw_weights[i-1]
      dG1 += dvdl
    dG1 = round(dG1, 2)

    dG2 = 0.0
    for i in range(1, ti_chg_windows+1):
      dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                      "| tail -n %d | "
                      "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                      %("DV\/DL", str(i)+'_chg_fwd1.out', "DV\/DL", chg_sample_steps/500)).read()
      dvdl = float(dvdl) * chg_weights[i-1]
      dG2 += dvdl
    dG2 = round(dG2, 2)

    if rev == 1:
      dG3 = 0.0
      for i in range(1, ti_chg_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_chg_bwd1.out', "DV\/DL", chg_sample_steps/500)).read()
        dvdl = float(dvdl) * chg_weights[i-1]
        dG3 += dvdl
      dG3 = round(dG3, 2)

      dG4 = 0.0
      for i in range(1, ti_vdw_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+'_vdw_bwd1.out', "DV\/DL", vdw_sample_steps/500)).read()
        dvdl = float(dvdl) * vdw_weights[i-1]
        dG4 += dvdl
      dG4 = round(dG4, 2)

      dG = (dG1 + dG2 - dG3 - dG4)/2
      dG = round(dG, 2)
    else:
      dG = dG1 + dG2
      dG = round(dG, 2)

    return dG


