import urllib
from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtlib.lib import get_lib_dict
from pymsmtmol.cal import calc_bond
import os

def downloadpdb(pdbid):
    httpadd = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=' + pdbid
    pdbfile = pdbid + '.pdb'
    urllib.urlretrieve(httpadd, pdbfile)

def get_names(fname):
    fp = open(fname, 'r')
    for line in fp:
        line = line.strip('\n')
        names = line.split(',')
        finalnames = [name.strip(' ') for name in names]
    return finalnames

def get_distance(atompairs, pdbids, ffchoice, mol2fs):
    #Lib dictionary
    libdict, chargedict = get_lib_dict(ffchoice)
    for mol2f in mol2fs:
      libdict1, chargedict1 = get_lib_dict(mol2f)
      libdict.update(libdict1)
      chargedict.update(chargedict1)

    #Distance dictionary
    disdict = {}

    for pdbid in pdbids:
        pdbid = pdbid.upper()

        #Download pdb file
        downloadpdb(pdbid)
        chains = os.popen("awk \'$1==\"ATOM\" || $1==\"HETATM\"\' %s.pdb | awk \'{print substr($0, 22, 1)}\' | awk \'!a[$0]++\'" %(pdbid)).read()
        chains = chains.strip('\n')
        chains = chains.split('\n')

        for chainid in chains:
            os.system("awk \'$1==\"ATOM\" || $1==\"HETATM\"\' %s.pdb | awk \'$5==\"%s\"\' > %s_%s.pdb" %(pdbid, chainid, pdbid, chainid))
 
            #Readpdb
            pdbname = pdbid + '_' + chainid + '.pdb'
            mol, atids, resids = get_atominfo_fpdb(pdbname)
 
            for i in range(0, len(atids)):
                crdi = mol.atoms[atids[i]].crd
                resnamei = mol.atoms[atids[i]].resname
                atnamei =  mol.atoms[atids[i]].atname
  
                if atnamei == 'OXT':
                    atnamei = 'O'
  
                if resnamei in chargedict.keys():
                    atypi = libdict[resnamei + '-' + atnamei][0]
                    chgi = libdict[resnamei + '-' + atnamei][1]
  
                    for j in range(i+1, len(atids)):
                        crdj = mol.atoms[atids[j]].crd
                        resnamej = mol.atoms[atids[j]].resname
                        atnamej = mol.atoms[atids[j]].atname
  
                        if atnamej == 'OXT':
                            atnamej = 'O'
  
                        if resnamej in chargedict.keys():
                            atypj = libdict[resnamej + '-' + atnamej][0]
                            chgj = libdict[resnamej + '-' + atnamej][1]
  
                            atypij = (atypi, atypj)
  
                            if atypij in atompairs or atypij[::-1] in atompairs:
                                disij = calc_bond(crdi, crdj)
  
                                if (atypij not in disdict.keys()) and (atypij[::-1] not in disdict.keys()):
                                    disdict[atypij] = [disij]
                                elif (atypij in disdict.keys()):
                                    disl = disdict[atypij]
                                    disl.append(disij)
                                    disdict[atypij] = disl
                                elif (atypij[::-1] in disdict.keys()):
                                    disl = disdict[atypij[::-1]]
                                    disl.append(disij)
                                    disdict[atypij[::-1]] = disl
  

    return disdict

def print_dict(dicti, fname):
    fp = open(fname, 'w')
    for i in dicti.keys():
        print >> fp, "##", i
        for j in dicti[i]:
            print >> fp, j
    fp.close()


