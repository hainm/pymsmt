#!/usr/bin/env python
import os
from optparse import OptionParser
import numpy
import math
from scipy.optimize import fmin #as fmin
from ipmach import calhfe
from ipmach import caliod
from ipmach import tifiles

#------------------------------------------------------------------------------
# Other functions
#------------------------------------------------------------------------------
def get_ep(rmin):
    rmin = round(rmin, 3)
    ep = 10**(-57.36* math.exp(-2.471 * rmin))
    ep = round(ep, 8)
    return rmin, ep

class ION:
    def __init__(self, resname, atname, attype, element, charge, rmin, ep):
        self.resname = resname
        self.atname = atname
        self.attype = attype
        self.element = element
        self.charge = charge
        self.rmin = rmin
        self.ep = ep

def get_dG():
    if prog == 'sander':
      if TISteps == 1:
        dG = calhfe.OneStep_sTI(ti_windows, ti_window_steps, ti_sample_steps, minexe, exe,
                                md_prmtop, md_inpcrd, md0_prmtop, md0_inpcrd,
                                ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4)
      elif TISteps == 2:
        dG = calhfe.TwoStep_sTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
                                chg_window_steps, vdw_sample_steps, chg_sample_steps,
                                minexe, exe, md0_prmtop, md0_inpcrd, mdv_prmtop, mdv_inpcrd, md_prmtop,
                                ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4)
    elif prog == 'pmemd':
      if TISteps == 1:
        dG = calhfe.OneStep_pTI(ti_windows, ti_window_steps,
                                ti_sample_steps, exe, ti_prmtop, ti_inpcrd,
                                ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4)
      elif TISteps == 2:
        dG = calhfe.TwoStep_pTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
                                chg_window_steps, vdw_sample_steps, chg_sample_steps,
                                exe, vdw_prmtop, ti_prmtop, ti_inpcrd, ti_min_steps, ti_nvt_steps,
                                ti_npt_steps, rev, ifc4)
    return dG

def get_IOD():
    iod, cn = caliod.MD_simulation(exe, md_prmtop, md_inpcrd, md_min_steps, md_nvt_steps,
                            md_npt_steps, md_md_steps, ifc4)
    return iod, cn

def gethfeiod(params):

    rmin, c4v = params[0], params[1]

    ion1.rmin, ion1.ep = get_ep(rmin)

    tifiles.gene_topcrd(ion0, ion1, 1, c4v)

    dG = get_dG()
    iod, cn = get_IOD()

    print "This result of this cycle:"
    print "    rmin=%5.3f, ep=%10.8f, c4=%5.0f " %(ion1.rmin, ion1.ep, c4v)
    print "    HFE=%6.1f, IOD=%5.2f, CN=%3.1f" %(dG, iod, cn)
    print "    Exp:HFE=%6.1f, IOD=%5.2f" %(HFE_VAL, IOD_VAL)

    HFEerr = abs(dG - HFE_VAL)
    IODerr = abs(iod - IOD_VAL)

    TOTerr = HFEerr + IODerr * 100.0

    if HFEerr <= hfetol and IODerr <= iodtol:
      print "Find the parameters!"
      print "    rmin=%5.3f, ep=%10.8f" %(ion1.rmin, ion1.ep)
      print "    HFE=%6.1f, IOD=%5.2f, CN=%3.1f" %(dG, iod, cn)
      print "    Exp:HFE=%6.1f, IOD=%5.2f" %(HFE_VAL, IOD_VAL)
      quit()

    return TOTerr

def gethfe(rmin):
    ion1.rmin, ion1.ep = get_ep(rmin)
    tifiles.gene_topcrd(ion0, ion1)
    dG = get_dG()

    print "This result of this cycle:"
    print "    rmin=%5.3f, ep=%10.8f" %(ion1.rmin, ion1.ep)
    print "    HFE=%6.1f" %dG

    HFEerr = abs(dG - HFE_VAL)

    if HFEerr <= hfetol:
      print "Find the parameters!"
      iod, cn = get_IOD()
      print "    rmin=%5.3f, ep=%10.8f" %(ion1.rmin, ion1.ep)
      print "    HFE=%6.1f, IOD=%5.2f, CN=%3.1f" %(dG, iod, cn)
      print "    Exp:HFE=%6.1f, IOD=%5.2f" %(HFE_VAL, IOD_VAL)
      quit()

    return HFEerr

def getiod(rmin):

    ion1.rmin, ion1.ep = get_ep(rmin)
    tifiles.gene_topcrd(ion0, ion1)
    iod, cn = get_IOD()

    print "This result of this cycle:"
    print "    rmin=%5.3f, ep=%10.8f" %(ion1.rmin, ion1.ep)
    print "    IOD=%5.2f, CN=%3.1f" %(iod, cn)

    IODerr = abs(iod - IOD_VAL)

    if IODerr <= iodtol:
      print "Find the parameters!"
      dG = get_dG()
      print "    rmin=%5.3f, ep=%10.8f" %(ion1.rmin, ion1.ep)
      print "    HFE=%6.1f, IOD=%5.2f, CN=%3.1f" %(dG, iod, cn)
      print "    Exp:HFE=%6.1f, IOD=%5.2f" %(HFE_VAL, IOD_VAL)
      quit()

    return IODerr

#----------------------------------------------------------------------------#
#                            Main Program                                    #
#----------------------------------------------------------------------------#
parser = OptionParser("usage: -i inputfile")
parser.add_option("-i", dest="inputf", type='string',
                  help="Input file name")
(options, args) = parser.parse_args()

#---------------------------Default values------------------------------------
cpus = 4
gpus = 0
TIsteps = 2
ti_windows = 7
ti_vdw_windows = 3
ti_chg_windows = 7
rev = 1
Cal = 'HFE'
max_iternum = 100
mode = 'normal'
c4v = 0.0
prog = 'pmemd'
hfetol = 1.0
iodtol = 0.01

#----------------------------Read the input files------------------------------
rinput = open(options.inputf, 'r')
for line in rinput:
  line = line.split()
  if line[0].lower() == "resname":
    resname = line[1]
  elif line[0].lower() == "atname":
    atname = line[1]
  elif line[0].lower() == "element":
    element = line[1]
  elif line[0].lower() == "attype":
    attype = line[1]
  elif line[0].lower() == "charge":
    charge = round(float(line[1]), 1)
  elif line[0].lower() == "program":
    prog = line[1].lower()
  elif line[0].lower() == "cpus":
    cpus = int(line[1])
  elif line[0].lower() == "gpus":
    gpus = int(line[1])
  elif line[0].lower() == "tisteps":
    TISteps = int(line[1])
  elif line[0].lower() == "ti_windows":
    ti_windows = int(line[1])
  elif line[0].lower() == "vdw_windows":
    ti_vdw_windows = int(line[1])
  elif line[0].lower() == "chg_windows":
    ti_chg_windows = int(line[1])
  elif line[0].lower() == "rev":
    rev = int(line[1])
  elif line[0].lower() == "rmin":
    rmin = round(float(line[1]), 3)
  elif line[0].lower() == "hfe":
    HFE_VAL = round(float(line[1]), 1)
  elif line[0].lower() == "iod":
    IOD_VAL = round(float(line[1]), 2)
  elif line[0].lower() == "set":
    Cal = line[1].upper()
  elif line[0].lower() == "maxiter":
    max_iternum = int(line[1])
  elif line[0].lower() == "mode":
    mode = line[1].lower()
  elif line[0].lower() == "hfetol":
    hfetol = float(line[1])
  elif line[0].lower() == "iodtol":
    iodtol = float(line[1])
  elif line[0].lower() == "c4":
    ifc4 = 1
    c4v = float(line[1])
    c4v = round(c4v, 1)
rinput.close()

#-------------------Setting for the program to use---------------------------
if gpus == 1:
  exe = 'pmemd.cuda'
else:
  if cpus == 1:
    exe = prog
    if exe == 'sander':
      raise ValueError('Could not perform sander TI calculation with one cpu')
  elif cpus > 1:
    exe = 'mpirun -np %d %s.MPI' %(cpus, prog)
    if prog == 'sander':
      minexe = 'mpirun -np 2 %s.MPI' %(prog)

if Cal == '1264' and exe == 'pmemd':
  raise ValueError('Could not perform 12-6-4 TI calculation with pmemd!')

#--------------------------------Setting for TI-------------------------------
ion0 = ION('M0', 'M0', 'M0', element, 0.0, 0.0, 0.0)
rmin, ep = get_ep(rmin)
ion1 = ION(resname, atname, attype, element, charge, rmin, ep)

ti_prmtop = element + "_wat_pti.prmtop"
ti_inpcrd = element + "_wat_pti.inpcrd"
vdw_prmtop = element + "_wat_pvdw.prmtop"

#For TI min-nvt-npt preparation
if mode == "test":
  ti_min_steps = 500
  ti_nvt_steps = 500
  ti_npt_steps = 500
elif mode == "normal":
  ti_min_steps = 2000
  ti_nvt_steps = 500000
  ti_npt_steps = 500000

if TISteps == 1:
  if mode == "test":
    ti_window_steps = 2000
    ti_sample_steps = ti_window_steps * 3/4
  elif mode == "normal":
    ti_window_steps = 200000
    ti_sample_steps = ti_window_steps * 3/4
elif TISteps == 2:
  if mode == "test":
    vdw_window_steps = 2000
    vdw_sample_steps = vdw_window_steps * 3/4
    chg_window_steps = 2000
    chg_sample_steps = chg_window_steps * 3/4
  elif mode == "normal":
    vdw_window_steps = 200000
    vdw_sample_steps = vdw_window_steps * 3/4
    chg_window_steps = 200000
    chg_sample_steps = chg_window_steps * 3/4

#--------------------------------Setting for MD-------------------------------
md0_prmtop = ion0.atname + "_wat_md.prmtop"
md0_inpcrd = ion0.atname + "_wat_md.inpcrd"

mdv_prmtop = ion0.atname + "_wat_vdw.prmtop"
mdv_inpcrd = ion0.atname + "_wat_vdw.inpcrd"

md_prmtop = element + "_wat_md.prmtop"
md_inpcrd = element + "_wat_md.inpcrd"

if mode == "test":
  md_min_steps = 500 
  md_nvt_steps = 500 
  md_npt_steps = 500 
  md_md_steps = 2000 
elif mode == "normal":
  md_min_steps = 2000
  md_nvt_steps = 500000
  md_npt_steps = 500000
  md_md_steps = 2000000

#--------------------------------Doing the optimization------------------------
if Cal == "1264":
  initparams = [rmin, c4v]
  rminc4opt = fmin(gethfeiod, initparams)#, epsilon=0.100)
elif Cal == "HFE":
  rminopt = fmin(gethfe, rmin)#, epsilon=0.100)
elif Cal == "IOD":
  rminopt = fmin(getiod, rmin)#, epsilon=0.100)


