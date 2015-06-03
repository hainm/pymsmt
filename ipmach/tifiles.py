#------------------------------------------------------------------------------
#-------- Generate the TI files
#------------------------------------------------------------------------------
import os

#Mass
Mass = {'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.06,
        'P': 30.97, 'F': 19.00, 'Cl': 35.45, 'Br': 79.90, 'I': 126.9,
        'Li': 6.94, 'Na': 22.99, 'K': 39.10, 'Rb': 85.47, 'Cs': 132.91,
        'Be': 9.01,'Cu': 63.55,'Ni': 58.69,'Pt': 195.08,'Zn': 65.4,
        'Co': 58.93,'Pd': 106.42,'Ag': 107.87,'Cr': 52.00,'Fe': 55.85,
        'Mg': 24.305,'V':  50.94,'Mn': 54.94,'Hg': 200.59,'Cd': 112.41,
        'Yb': 173.05,'Ca': 40.08,'Sn': 118.71,'Pb': 207.2,'Eu': 151.96,
        'Sr': 87.62,'Sm': 150.36,'Ba': 137.33,'Ra': 226.03,
        'Al': 26.98,'Fe': 55.85,'Cr': 52.00,'In': 114.82,'Tl': 204.38,
        'Y': 88.91,'La': 138.91,'Ce': 140.12,'Pr': 140.91,'Nd': 144.24,
        'Sm': 150.36,'Eu': 151.96,'Gd': 157.25,'Tb': 158.93,'Dy': 162.5,
        'Er': 167.26,'Tm': 168.93,'Lu': 174.97, 'As': 74.92, 'Ru': 101.07,
        'Hf': 178.49,'Zr': 91.22,'Ce': 140.12,'U': 238.03,'Pu': 244.06,
        'Th': 232.04, 'Mo': 95.96
        }

def write_frcmod(ion, fname):
    mass = Mass[ion.element]
    frcmodf = open(fname, 'w')
    print >> frcmodf, "# remark goes here"
    print >> frcmodf, "MASS"
    print >> frcmodf, "%s  %5.2f" %(ion.attype, mass)
    print >> frcmodf, " "
    print >> frcmodf, "BOND"
    print >> frcmodf, " "
    print >> frcmodf, "ANGL"
    print >> frcmodf, " "
    print >> frcmodf, "DIHE"
    print >> frcmodf, " "
    print >> frcmodf, "NONB"
    print >> frcmodf, "%s    %5.3f    %10.8f" %(ion.attype, ion.rmin, ion.ep)
    print >> frcmodf, " "
    frcmodf.close()

def write_cmd(ion):
    cmdf = open(ion.resname+'.cmd', 'w')
    print >> cmdf, "i = createAtom %s %s %3.1f" %(ion.atname, ion.attype, ion.charge)
    print >> cmdf, "set i element \"%s\"" %ion.element
    print >> cmdf, "set i position { 0 0 0 }"
    print >> cmdf, "r = createResidue %s" %ion.resname
    print >> cmdf, "add r i"
    print >> cmdf, "%s = createUnit %s" %(ion.resname, ion.resname)
    print >> cmdf, "add %s r" %ion.resname
    print >> cmdf, "saveOff %s ./%s.lib" %(ion.resname, ion.resname)
    print >> cmdf, "quit"
    cmdf.close()

def write_leapin(ion0, ion1):
    leapf = open('tleap.in', 'w')
    print >> leapf, "source leaprc.ff14SB"
    print >> leapf, "loadoff %s.lib" %ion0.resname
    print >> leapf, "ion = loadpdb ION.pdb"
    print >> leapf, "solvatebox ion OPCBOX 13"
    print >> leapf, "loadamberparams frcmod.opc"
    print >> leapf, "loadamberparams %s.frcmod" %ion0.attype  #VDW=0
    print >> leapf, "loadamberparams %s.frcmod" %ion1.attype  #VDW!=0
    print >> leapf, "ion2 = copy ion"
    print >> leapf, "remove ion2 ion2.2"
    print >> leapf, "saveamberparm ion2 %s_wat_md.prmtop %s_wat_md.inpcrd" %(ion0.atname, ion0.atname) #For sander: with M0 without VDW
    print >> leapf, "saveamberparm ion %s_wat_ti.prmtop %s_wat_ti.inpcrd" %(ion0.element, ion0.element) #For pmemd: with Na has VDW and M0 doesn't have VDW
    print >> leapf, "loadamberparams %s_vdw.frcmod" %ion0.attype #VDW != 0
    print >> leapf, "saveamberparm ion %s_wat_vdw.prmtop %s_wat_vdw.inpcrd" %(ion0.element, ion0.element) #For pmemd: with Na and M0 all have VDW
    print >> leapf, "remove ion ion.1"
    print >> leapf, "saveamberparm ion %s_wat_md.prmtop %s_wat_md.inpcrd" %(ion0.element, ion0.element) #For md and sander with Na+ has charge and VDW
    print >> leapf, "saveamberparm ion2 %s_wat_vdw.prmtop %s_wat_vdw.inpcrd" %(ion0.atname, ion0.atname) #For sander, with M0 has VDW
    print >> leapf, "quit"
    leapf.close()

def addc4(ion0, ion1, c4v):

    #Without charge and VDW
    c4f0 = open(ion0.atname + '_c4.txt', 'w')
    print >> c4f0, "%s %f" %(ion0.element+str(int(ion0.charge)), 0.0)
    c4f0.close()

    parmf = open('md0_addc4.in', 'w')
    print >> parmf, "loadRestrt %s_wat_md.inpcrd" %ion0.atname
    print >> parmf, "setOverwrite True"
    print >> parmf, "add12_6_4 :1@%s c4file %s_c4.txt" %(ion0.atname, ion0.atname)
    print >> parmf, "outparm %s_wat_md.prmtop %s_wat_md.inpcrd" %(ion0.atname, ion0.atname)
    parmf.close()
    os.system("parmed.py -i md0_addc4.in -p %s_wat_md.prmtop | tail -n 3" %ion0.atname)

    #With VDW and witout charge
    parmf = open('mdv_addc4.in', 'w')
    print >> parmf, "loadRestrt %s_wat_vdw.inpcrd" %ion0.atname
    print >> parmf, "setOverwrite True"
    print >> parmf, "add12_6_4 :1@%s c4file %s_c4.txt" %(ion0.atname, ion0.atname)
    print >> parmf, "outparm %s_wat_vdw.prmtop %s_wat_vdw.inpcrd" %(ion0.atname, ion0.atname)
    parmf.close()
    os.system("parmed.py -i mdv_addc4.in -p %s_wat_vdw.prmtop | tail -n 3" %ion0.atname)

    #With Charge and VDW
    c4f1 = open(ion1.element + '_c4.txt', 'w')
    print >> c4f1, "%s %f" %(ion1.element+str(int(ion1.charge)), c4v)
    c4f1.close()

    parmf = open('md_addc4.in', 'w')
    print >> parmf, "loadRestrt %s_wat_md.inpcrd" %ion1.element
    print >> parmf, "setOverwrite True"
    print >> parmf, "add12_6_4 :1@%s c4file %s_c4.txt" %(ion1.atname, ion1.element)
    print >> parmf, "outparm %s_wat_md.prmtop %s_wat_md.inpcrd" %(ion1.element, ion1.element)
    parmf.close()
    os.system("parmed.py -i md_addc4.in -p %s_wat_md.prmtop | tail -n 3" %ion1.element)

def gene_topcrd(ion0, ion1, ifc4=0, c4v=0.0):

    #os.system("module load AmberTools/15")

    print "Perform the system..."
    #lib file
    write_cmd(ion0)
    write_cmd(ion1)
    os.system("tleap -s -f %s.cmd > %s.out" %(ion0.resname, ion0.resname))
    os.system("tleap -s -f %s.cmd > %s.out" %(ion1.resname, ion1.resname))
    #pdb file
    pdbf = open('ION.pdb', 'w')
    print >> pdbf, "HETATM 2032 %2s    %2s A   1       0.000   0.000   0.000  1.00  0.00" %(ion0.atname, ion0.resname)
    print >> pdbf, "TER"
    print >> pdbf, "HETATM 2032 %2s    %2s A   1       0.000   0.000   0.000  1.00  0.00" %(ion1.atname, ion1.resname)
    pdbf.close()
    #frcmod file
    write_frcmod(ion0, ion0.attype + '.frcmod')
    write_frcmod(ion1, ion1.attype + '.frcmod')
    ion0.rmin = ion1.rmin
    ion0.ep = ion1.ep
    write_frcmod(ion0, ion0.attype + '_vdw.frcmod')
    #tleap file
    write_leapin(ion0, ion1)
    os.system("tleap -s -f tleap.in > tleap.out")
    #add c4 parameters
    if ifc4 == 1:
      addc4(ion0, ion1, c4v)

    #os.system("module load Amber/14")

