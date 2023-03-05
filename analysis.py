import pandas as pd
import re
from dictionaries import l_map, l_pam, orbital_counts
import numpy as np


def open_file(filename) -> list[str]:
    """Open the file and return the data as a list of strings."""
    with open(filename, 'r') as f:
        data = f.readlines()
    return data


def get_orbitals(filename: str):
    ''' Createw a orbitals.txt file which contains alpha, beta, and natural orbitals'''
    import os
    os.system(
        f'cat {filename} | sed -n \'/Index Energy  Occupation Coefficients/,/--/p\' >> orbitals.txt')
    print('Orbitals saved to orbitals.txt')


def split_file(filename):
    '''Split orbitals data into multiple files'''
    data = open_file(filename)
    # join the data into a single string
    data = ''.join(data)
    # the delimiter is --
    data = data.split('--')

    # remove data[-1]
    data = data[:-1]

    titles = ['uhf_alpha.txt', 'uhf_beta.txt', 'natural_orbitals.txt']
    # save the files
    for i, section in enumerate(data):
        with open(f'{titles[i]}', 'w') as f:
            f.write(section)

    # remove the first line of the second and third files
    for i in range(1, 3):
        with open(f'{titles[i]}', 'r') as f:
            data = f.readlines()
        data = data[1:]
        with open(f'{titles[i]}', 'w') as f:
            f.writelines(data)

def split_ras_file(filename: str):
    '''Split ras orbitals data into multiple files'''
    data = open_file(filename)
    # join the data into a single string
    data = ''.join(data)
    # the delimiter is --
    data = data.split('--')

    #keep only data[-2]
    data = data[-2]

    #remove the first line
    data = data[1:]

    #write the data to a file
    with open('ras_orbitals.txt', 'w') as f:
        f.write(data)

    # replace the first instance of '1-(?=\d)' with '1 -'
    with open('ras_orbitals.txt', 'r') as f:
        data = f.read()
    data = re.sub(r'1-(?=\d)', '1 -', data)
    with open('ras_orbitals.txt', 'w') as f:
        f.write(data)


        

def separate_scf_mos(data: list[str]):
    '''Seperate the data into a list of strings'''

    # data is currently a list of strings, but we want it to be a single string
    data = ''.join(data)

    # the data is seperated by the regex ^\s*$ which is a blank line
    mos = re.split(r'^\s*$', data, flags=re.MULTILINE)

    # remove the first and last elements of the list
    mos = mos[1:-1]

    return mos

def seperate_ras_mos(data: list[str]):
    '''Seperate the data into a list of strings'''

    data = ''.join(data)

    mos = re.split(r'(^\s*\d+\s+[-\d]+\.\d+\s+\d+\.\d+)', data, flags=re.MULTILINE)
    
    # remove the first element of the list
    mos = mos[1:]

    #fuse n and n+1 elements of the list
    mos = [f'{mos[i]} {mos[i+1]}' for i in range(0, len(mos)-1, 2)]

    return mos

def format_mo(mo):
    '''Clean a single mo and return a dictionary of the data'''

    # if there is a substring 1- , seperate the 1 and the - by a space, but only if it is followed by a digit
    mo = re.sub(r'1-(?=\d)', '1 -', mo)

    # split the mo by whitespace
    mo = mo.split()

    # remove parentheses
    mo = [x.replace('(', '') for x in mo]

    # remove commas and empty strings
    mo = [x for x in mo if x != ',' and x != '']

    # remove comma substrings from each element
    mo = [x.replace(',', '') for x in mo]

    # remove ) from each element in the list
    mo = [x.replace(')', '') for x in mo]

    return mo


def get_mos(file: str) -> list:
    data = open_file(file)

    #there are two types of mos, ras and scf and there is a different function for each
    #The header of the ras mos is INDEX  ENERGY  OCCUPATION COEFFICIENTS ...
    #The header of the scf mos is Index Energy  Occupation Coefficients ...
    #check the header of the file to determine which function to use
    if 'INDEX' in data[0]:
        mos = seperate_ras_mos(data)
    elif 'Index' in data[0]:
        mos = separate_scf_mos(data)
    else:
        raise Exception('ERROR: File not formatted correctly')
    # format the mos
    mos = [format_mo(mo) for mo in mos]

    # remove empty lists
    mos = [mo for mo in mos if mo != []]

    return mos



class Coefficient:

    # creat a dictionary of the atomic orbital types
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}

    def __init__(self, coeff: tuple):

        self.ao_index = coeff[0]
        self.atom_type = coeff[1]
        self.ao_type = coeff[2]
        self.value = float(coeff[3])
        self.n = self.get_n()
        self.l_int, self.l_sym = self.get_l()

    def __str__(self):
        return f'Atomic Orbital {self.ao_index} of type {self.ao_type} from atom {self.atom_type} with coefficient {self.value}'

    def __repr__(self):
        return f'Atomic Orbital {self.ao_index} of type {self.ao_type} from atom {self.atom_type} with coefficient {self.value}'

    def get_n(self):

        # the n value is the first character of the ao_type
        try:
            n = int(self.ao_type[0])
        except ValueError:
            print(f'Diagnostics: {self.ao_type} is {type(self.ao_type)}')
            print(f'Diagnostics: ao index {self.ao_index}')
            print('ERROR: ao_type is not formatted correctly')
        return n

    def get_l(self):

        # the l value is the second character of the ao_type
        l_sym = str(self.ao_type[1])

        # convert the l value to an int
        l_int = int(self.l_map[l_sym])

        return l_int, l_sym


class molecularOrbital:

    def __init__(self, mo):
        self.mo = mo
        self.index = mo[0]
        self.energy = float(mo[1])
        self.occupation = float(mo[2])
        self.coefficients = self.get_coefficents(mo)
        self.composition = self.get_composition(by_l=True)
        self.elemental_composition = self.get_elemental_composition()
        self.active_space = None

    def get_coefficents(self, mo):

        coefficients = mo[3:]

        # separated in groups of 4 -> len(coefficients) % 4 == 0 and len(coefficients) // 4 == number of groups
        # each group is a tuple of 4 elements
        coefficients = [tuple(coefficients[i:i+4])
                        for i in range(0, len(coefficients), 4)]
        coefficients = [Coefficient(coeff) for coeff in coefficients]

        return coefficients

    def get_composition(self, by_l: bool = False):
        '''Calculates the percentage of orbital composition by n '''

        if by_l:

            # create a tuple of coeff.n , coeff.l_int and coeff.value^2
            composition = [(coeff.n, coeff.l_int, coeff.value**2)
                           for coeff in self.coefficients]

            # sort the tuple by n and l_int
            composition = sorted(composition, key=lambda x: (x[0], x[1]))

            # sum the squared values of matching n and l_int
            composition = {n: {l: sum([x[2] for x in composition if x[0] == n and x[1] == l]) for l in range(
                0, max([x[1] for x in composition]) + 1)} for n in range(1, max([x[0] for x in composition]) + 1)}

            return composition

        else:
            # create a tuple of coeff.n and coeff.value^2
            composition = [(coeff.n, coeff.value**2)
                           for coeff in self.coefficients]

            # sort the tuple by n
            composition = sorted(composition, key=lambda x: x[0])

            # create a dictionary of n and the sum of the squared values
            composition = {n: sum([x[1] for x in composition if x[0] == n]) for n in range(
                1, max([x[0] for x in composition]) + 1)}

            return composition

    def get_elemental_composition(self):

        # create a list of atoms from the coefficients
        atoms = [coeff.atom_type for coeff in self.coefficients]

        # remove digits from the atoms
        atoms = [re.sub(r'\d', '', atom) for atom in atoms]

        # if the element has more than one letter, make the second letter lowercase
        atoms = [atom[0] + atom[1:].lower() for atom in atoms]

        # get the percentage of each atom
        elemental_composition = {atom: atoms.count(
            atom) / len(atoms) for atom in atoms}

        return elemental_composition

    def __str__(self):

        return f'Orbital {self.index} with energy {self.energy} and occupation {self.occupation}'
    
    def set_active_space(self, active_space: list[int]):
        self.active_space = active_space


class MolecularManifold:

    def __init__(self, file):

        self.mos = get_mos(file)
        self.mos = [molecularOrbital(mo) for mo in self.mos]
        self.energies = [float(mo.energy) for mo in self.mos]
        self.coefficients = [mo.coefficients for mo in self.mos]
        self.homo, self.homo_index = self.get_homo()
        self.lumo, self.lumo_index = self.mos[self.homo_index + 1], self.homo_index + 1
        # self.occupation is the total of the occupation of the orbitals
        self.occupation = sum([float(mo.occupation) for mo in self.mos])

    def get_homo(self):

        # get the index of the molecular orbital with the highest index and occupation > 0
        homo_index = max([i for i, mo in enumerate(self.mos) if mo.occupation > 0])

        return self.mos[homo_index], homo_index

    def __str__(self):
        return f'Molecular Manifold with {len(self.mos)} molecular orbitals'

    def to_dataframe(self):
        ''' Returns a pandas dataframe with the molecular orbitals '''

        # create a dataframe with the energies and occupations
        df = pd.DataFrame({'Energy': self.energies, 'Occupation': [
                          mo.occupation for mo in self.mos]})

        
        for mo in self.mos: 
            orbital_character_dict = {'s': 0, 'p': 0, 'd': 0, 'f': 0, 'g': 0, 'h': 0}
            n_max = max(mo.composition.keys())
            for l in mo.composition[n_max].keys():
                orbital_character_dict[l_pam[l]] = mo.composition[n_max][l]

        for i in orbital_character_dict.keys():
            df[i] = 0 
        
        for i, mo in enumerate(self.mos):
            n_max = max(mo.composition.keys())
            for l in mo.composition[n_max].keys():
                df.loc[i, l_pam[l]] = mo.composition[n_max][l]

        df_subset = df[['s', 'p', 'd', 'f', 'g', 'h']]
        df_subset = df_subset.div(df_subset.sum(axis=1), axis=0)
        #set precision to 3rd decimal
        df_subset = df_subset.round(3)
        df[['s', 'p', 'd', 'f', 'g', 'h']] = df_subset
        
        return df

    def active_space_indices(self, active_space_orbitals: list[str]) -> list[int]:
        ''' Returns the indices of the active space orbitals '''

        #each element in the array is the n and l of the orbital
        #regex for a single letter is 
        active_space_orbitals = [re.findall(r'[a-z]', orbital)[0] for orbital in active_space_orbitals]
        
        #map the letter to the number of orbitals
        n_orbitals = sum([orbital_counts[orbital] for orbital in active_space_orbitals])

        #return a list that goes from homo_index - n_orbitals to homo_index, incremented by 1
        indices = np.arange(self.homo_index + 1 - n_orbitals, self.homo_index + 1, 1)

        return indices
    
    def active_space_composition(self, active_space_indices: list[int]):
        '''Returns the orbital composition of the active space indices'''
        
        dataframe = self.to_dataframe()
        active_space_orbitals = dataframe.iloc[active_space_indices]
        #sum the composition of the orbitals and report it as a percentage of s, p, d, f, g, h
        active_space_composition = active_space_orbitals[['s', 'p', 'd', 'f', 'g', 'h']].sum()
        active_space_composition = active_space_composition/active_space_composition.sum()

        return active_space_composition




