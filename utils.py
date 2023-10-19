# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""
from math import sqrt, isclose
import random
import json


def is_valid_struct(struct : list) -> bool:
    '''
    Check if the structure inserted is valid: is SAW (self avoid walk) and distances between consecutive elements are 1.

    Parameters
    ----------
    struct : list
        Structur of the protein containing x and y coordinate in a list.

    Returns
    -------
    bool
        True if the structure is valid, False if is not.
    '''
    unique_struct = [] # counter of the monomer positions
    n = len(struct) # length of the sequence
    
    for i in range(n):
        if struct[i] in unique_struct:
            return False
        else:
            unique_struct.append(struct[i])
        
        if (i<n-1): # check distance between the monomer i and the following one
            if not isclose(get_dist(struct[i], struct[i+1]),1):
                return False
            
    return True


def is_valid_sequence(seq : str) -> bool:
    '''
    Check if the protein sequence contains only H and/or P and its lengths is almost 3.

    Parameters
    ----------
    seq : str
        The protein sequence (caps sensitive).

    Returns
    -------
    bool
        True if the sequence is valid, False if is not.
    '''
    if len(seq) < 3:
        print('The sequence is too short. It must be at least 3.')
        return False
    
    unique_c = set(seq) # set of all the letters present in the sequence
    
    if len(unique_c) == 2:
        if 'H' in unique_c and 'P' in unique_c:
            return True
    elif len(unique_c) == 1:
        if 'H' in unique_c or 'P' in unique_c:
            return True   
        
    return False # if True is not returned yet, return Fasle since the sequence must be invalid 
    

def linear_struct(seq : str) -> list:
    '''
    Create the linear structure of the length of the input sequence

    Parameters
    ----------
    seq : str
        Sequence which need a linear structure
    
    Returns
    -------
    list : 
        Linear structure
    '''
    struct = []

    for i in range(len(seq)):
        struct.append([i,0])

    print('\033[43mLinear initial structure assumed \033[0;0m')
    
    return struct


def get_dist(coord1 : list, coord2 : list) -> float:
    '''
    Compute the distance of two points in the lattice

    Parameters
    ----------
    coord1 : list
        a and y coordinate in the lattice of the first monomer.
    coord2 : list
        a and y coordinate in the lattice of the second monomer.

    Returns
    -------
    float
        Euclidean distance as a float.
    '''
    dist = sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2)
    return dist


def diagonal_move(struct : list, previous : list) -> list:
    '''
    Move the first monomer along a diagonal looking at the previous and following monomers in the sequence. \n
    It is assumed that the protein structure starts in [0,0], so the coordinates of the surrounding monomers must have only zeros and ones. \n
    This function is called only if the structure can accept a diagonal move (surrounding monomers not aligned).

    Parameters
    ----------
    struct : list
        Protein structure starting in [0,0].
    previous : list
        x and y coordinates of the previous monomer shifted such that the first monomer of struct is [0,0]

    Returns
    -------
    list
        The structure with the first monomer moved.
    '''
    case = random.randint(0, 1) # random movement respect the second monomer
    x_prev, y_prev = previous # previous monomer coords
    x_foll, y_foll = struct[1] # following monomer coordinates

    struct[0] = [x_prev+x_foll,y_prev+y_foll]
    
    return struct


def tail_fold(struct : list, method : int, previous : list) -> list:
    '''
    Apply a rotation/inversion of symmetry at the sequence inserted.\n
    7 methods are present:
        1: 90° clockwise rotation
        2: 90° anticlockwise rotation
        3: 180° rotaion
        4: x-axis refletion
        5: y-axis reflection
        6: 1 and 3 quadrant bisector symmetry
        7: 2 and 4 quadrant bisector symmetry
        8: movement on a digaonal of a random monomer

    Parameters
    ----------
    struct : list
        Structure of the sequence for which each element is the x and y coordinates of the monomer.
    method : int
        The method to apply to the structure.
    previous : list
        x and y coordinates of the previous monomer shifted such that the first monomer of struct is [0,0]

    Returns
    -------
    list
        The structure transformed.
    ''' 
    new_tail = []
    
    if method == 1: # 90 rotation clockwise 
        for x,y in struct:
            new_tail.append([y,-x])
    if method == 2: # 90 rotation anticlockwise
        for x,y in struct:
            new_tail.append([-y,x])
    if method == 3: # 180 rotation
        for x,y in struct:
            new_tail.append([-x,-y])
    if method == 4: # x-axis refletion
        for x,y in struct:
            new_tail.append([x,-y])
    if method == 5: # y-axis reflection
        for x,y in struct:
            new_tail.append([-x,y])
    if method == 6: # 1 and 3 quadrant bisector symmetry
        for x,y in struct:
            new_tail.append([-y,-x])
    if method == 7: # 2 and 4 quadrant bisector symmetry
        for x,y in struct:
            new_tail.append([y,x])
    if method == 8: # movement on the digonal
        new_tail = diagonal_move(struct, previous)
                
    return new_tail


def hp_sequence_transform(seq : str) -> str :
    '''
    Transform a compleate sequence of 20 amino-acids into the HP sequence used in the code as model.

    Parameters
    ----------
    seq : list
        Sequence containing the 20 different amino-acids (RNDQEHKSTACGILMFPWYV), they must be upper case letters.

    Returns
    -------
    str
        The sequence converted into only H/P.  
    '''

    polar = 'RNDQEHKST' # polar amino acids 
    hydr = 'ACGILMFPWYV' # hydrophobic amino acids

    hp_seq ='' # string where to save the hp sequence

    for amin in seq:
        if amin in polar:
            hp_seq += 'P'
        elif amin in hydr:
            hp_seq += 'H'
        else:
            raise ValueError(f'Amino acids {amin} not recognized')
    
    return hp_seq


def progress_bar(progress : int, total : int) -> None:
    '''
    Print a progress bar on terminal when used inside a loop.

    Parameters
    ----------
    progress : int
        Enumeration during the loop to take counts of the progress during the process.
    total : int
        Total number of steps for the evolution

    Returns
    -------
    None
        It only print on terminal the progress bar
    '''
    percentage = progress / float(total) * 100 
    left = int(percentage/10) # number of  '#' to print (up to 10)
    rigth = 10 - left # number of blank spaces to be included
    bar = '[' + '#' * left + ' ' * rigth + ']' 
    print(f'\r{bar} {percentage:.2f}%',end='')
    if progress == total:
        print(f'\r{bar} {percentage:.2f}%') # to have a newline after the process is finished


class Configuration():
    '''
    Class to save and ordinate the parameters given in the input file.

    Parameters
    ----------
    config : ConfigParser
        ConfigParser of the desired input file.

    Returns
    -------
    None.
    '''
    def __init__(self, config) -> None:
        self.seq = config['SEQUENCE']['sequence'] # selected sequence
        self.folds = config['PROCESS'].getint('folding_steps') # number of folds
        self.use_struct = config['optional'].getboolean('use_structure') # if use the structure present in config file or use linear structure
        self.annealing = config['optional'].getboolean('annealing') # if use annealing or not
        self.T = config['optional'].getfloat('T') # starting temperature
        if self.use_struct:
            struct = config['optional']['structure'] # structure if TRUE in input file
            self.struct = json.loads(struct)
        self.gif = config['optional'].getboolean('create_gif')
        self.seed = config['random_seed']['seed'] # get the random seed
        if self.seed == 'None': # generate a random seed if None
            self.seed = random.randint(0,10000)
        self.seed = int(self.seed) # convers the seed to int in any case