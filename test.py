from analysis import *

manifold = MolecularManifold('ras_orbitals.txt')
indices = [192, 193, 194, 195, 196]
comp = manifold.active_space_composition(indices)
print(comp)


