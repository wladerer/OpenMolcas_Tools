from analysis import *

manifold = MolecularManifold('natural_orbitals.txt')
print(manifold.homo_index)
print(manifold.active_space_indices('f'))
