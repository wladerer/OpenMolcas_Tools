from analysis import *

manifold = MolecularManifold('natural_orbitals.txt')

print(manifold.to_dataframe())
