from analysis import *

manifold = MolecularManifold('natural_orbitals.txt')
df = manifold.to_dataframe()
df.to_csv('natural_orbitals.csv')