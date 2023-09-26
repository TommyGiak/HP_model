# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:14:08 2023

@author: Tommaso Giacometti
"""
from math import sqrt, isclose
import random


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
    else:
        return False
    
    
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


def diagonal_move(struct : list) -> list:
    '''
    Move the first monomer along a diagonal looking for the second monomer in the sequence. \n
    The information of the following monomer helps to the optimization of the algorithm. \n
    It is assumed that the protein structure starts in [0,0], so the second monomer must have a coordinate eqaul
    to +/- 1 and the other equal to zero.

    Parameters
    ----------
    struct : list
        Protein structure starting in [0,0].

    Returns
    -------
    list
        The structure with the first monomer moved.
    '''
    case = random.randint(0, 1) # random movement respect the second monomer
    x,y = struct[1] # second monomer coordinates

    # movement settings
    if x == 1:
        if case == 0:
            struct[0] = [1,1]
        if case == 1:
            struct[0] = [1,-1]
    elif x == -1:
        if case == 0:
            struct[0] = [-1,1]
        if case == 1:
            struct[0] = [-1,-1]
    elif y == 1:
        if case == 0:
            struct[0] = [1,1]
        if case == 1:
            struct[0] = [-1,1]
    elif y == -1:
        if case == 0:
            struct[0] = [1,-1]
        if case == 1:
            struct[0] = [-1,-1]
    else:
        raise RuntimeError('No compatible information about the second monomer')
    
    return struct

