from data_analysis import get_distance, get_names, print_dict

atompairs = [('C0', 'OW')]
pdbids = ['4xw4']
pdbids = get_list('pdbnames.txt')
ffchoice = 'ff14SB'
mol2fs = ['HOH.mol2', 'Ca.mol2']
disdict = get_distance(atompairs, pdbids, ffchoice, mol2fs)
print_dict(disdict, 'dis.txt')

