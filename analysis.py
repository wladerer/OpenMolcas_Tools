import pandas as pd
import re


def open_file(filename) -> list:
    """Open the file and return the data as a list of strings."""
    with open(filename, 'r') as f:
        data = f.readlines()
    return data

def split_file(filename):
    '''Split orbitals data into multiple files'''
    data = open_file(filename)
    #join the data into a single string
    data = ''.join(data)
    #the delimiter is -- 
    data = data.split('--')

    #remove data[-1]
    data = data[:-1]

    titles = ['uhf_alpha.txt', 'uhf_beta.txt', 'natural_orbitals.txt']
    #save the files 
    for i, section in enumerate(data):
        with open(f'{titles[i]}', 'w') as f:
            f.write(section)

    #remove the first line of the second and third files
    for i in range(1, 3):
        with open(f'{titles[i]}', 'r') as f:
            data = f.readlines()
        data = data[1:]
        with open(f'{titles[i]}', 'w') as f:
            f.writelines(data)


    


def separate_mos(data):
    '''Seperate the data into a list of strings'''

    # data is currently a list of strings, but we want it to be a single string
    data = ''.join(data)

    # the data is seperated by the regex ^\s*$ which is a blank line
    mos = re.split(r'^\s*$', data, flags=re.MULTILINE)

    # remove the first and last elements of the list
    mos = mos[1:-1]


    return mos

def format_mo(mo):
    '''Clean a single mo and return a dictionary of the data'''

    #if there is a substring 1- , seperate the 1 and the - by a space, but only if it is followed by a digit
    mo = re.sub(r'1-(?=\d)', '1 -', mo)

    # split the mo by whitespace
    mo = mo.split()

    # remove parentheses
    mo = [x.replace('(', '') for x in mo]

    # remove commas and empty strings
    mo = [x for x in mo if x != ',' and x != '']

    # remove comma substrings from each element
    mo = [x.replace(',', '') for x in mo]

    #remove ) from each element in the list
    mo = [x.replace(')', '') for x in mo]

    return mo

def get_mos(file: str) -> list:
    data = open_file(file)
    mos = separate_mos(data)
    mos = [format_mo(mo) for mo in mos]

    #remove empty lists
    mos = [mo for mo in mos if mo != []]

    return mos

class Coefficient:

    # creat a dictionary of the atomic orbital types
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g':4, 'h':5}

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

        #the n value is the first character of the ao_type
        n = int(self.ao_type[0])

        return n
    
    def get_l(self):

        #the l value is the second character of the ao_type
        l_sym = str(self.ao_type[1])

        #convert the l value to an int
        l_int = int(self.l_map[l_sym])

        return l_int, l_sym


class molecularOrbital:

    def __init__(self, mo):
        self.mo = mo
        self.index = mo[0]
        self.energy = mo[1]
        self.occupation = mo[2]
        self.coefficients = self.get_coefficents(mo)
        self.composition = self.get_composition(by_l=True)
        self.elemental_composition = self.get_elemental_composition()

    def get_coefficents(self, mo):

        coefficients = mo[3:]

        #separated in groups of 4 -> len(coefficients) % 4 == 0 and len(coefficients) // 4 == number of groups
        #each group is a tuple of 4 elements
        coefficients = [tuple(coefficients[i:i+4]) for i in range(0, len(coefficients), 4)]
        coefficients = [Coefficient(coeff) for coeff in coefficients]

        return coefficients
    
    def get_composition(self, by_l : bool = False):
        '''Calculates the percentage of orbital composition by n '''

        if by_l:

            #create a tuple of coeff.n , coeff.l_int and coeff.value^2
            composition = [(coeff.n, coeff.l_int, coeff.value**2) for coeff in self.coefficients]

            #sort the tuple by n and l_int
            composition = sorted(composition, key=lambda x: (x[0], x[1]))

            #sum the squared values of matching n and l_int
            composition = {n: {l: sum([x[2] for x in composition if x[0] == n and x[1] == l]) for l in range(0, max([x[1] for x in composition])+ 1)} for n in range(1, max([x[0] for x in composition])+ 1)}

            return composition

        else:
            #create a tuple of coeff.n and coeff.value^2
            composition = [(coeff.n, coeff.value**2) for coeff in self.coefficients]

            #sort the tuple by n
            composition = sorted(composition, key=lambda x: x[0])

            #create a dictionary of n and the sum of the squared values
            composition = {n: sum([x[1] for x in composition if x[0] == n]) for n in range(1, max([x[0] for x in composition])+ 1)}

            return composition
        
    def get_elemental_composition(self):

        # create a list of atoms from the coefficients
        atoms = [coeff.atom_type for coeff in self.coefficients]
        
        # remove digits from the atoms
        atoms = [re.sub(r'\d', '', atom) for atom in atoms]

        # if the element has more than one letter, make the second letter lowercase
        atoms = [atom[0] + atom[1:].lower() for atom in atoms]

        # get the percentage of each atom
        elemental_composition = {atom: atoms.count(atom) / len(atoms) for atom in atoms}

        return elemental_composition

        
    def __str__(self):

        return f'Orbital {self.index} with energy {self.energy} and occupation {self.occupation}'
    

class MolecularManifold:

    def __init__(self, file):
        
        self.mos = get_mos(file)
        self.mos = [molecularOrbital(mo) for mo in self.mos]
        self.energies = [float(mo.energy) for mo in self.mos]
        self.coefficients = [mo.coefficients for mo in self.mos]
        self.homo, self.homo_index = self.get_homo()
        self.lumo, self.lumo_index = self.get_lumo()
        #self.occupation is the total of the occupation of the orbitals
        self.occupation = sum([float(mo.occupation) for mo in self.mos])

    def get_homo(self):

        # the homo is the orbital that is the least negative in energy, but still negative
        # find the least negative energy
        homo_index = self.energies.index(max([energy for energy in self.energies if energy < 0]))

        return self.mos[homo_index], homo_index
    
    def get_lumo(self):

        # the lumo is the orbital that is the most positive in energy, but still positive
        # find the most positive energy
        lumo_index = self.energies.index(min([energy for energy in self.energies if energy > 0]))

        return self.mos[lumo_index], lumo_index

    def __str__(self):
        return f'Molecular Manifold with {len(self.mos)} molecular orbitals'
    

manifold = MolecularManifold('natural_orbitals.txt')
print(manifold.occupation)




