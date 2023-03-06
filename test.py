from analysis import *

#I have to streamline this process 
# step 1. is to convert the molcas log file into orbitals.txt file
# step 2. is to convert the orbitals.txt file into ras_orbitals.txt or the various scf files
# step 3, is to run the code below

manifold = MolecularManifold('ras_orbitals.txt')
l_dict = manifold.mos[200].l_composition

print(l_dict)


