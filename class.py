# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 15:54:27 2023

@author: Tommaso Giacometti
"""
from math import sqrt, isclose

class Protein():
    '''
    Protein class that contains: the sequence of the protein (Protein.seq), the length of the sequence 
    (Protein.n) and the structure of the protein (Protein.struc). \n
    When initialized it automatically check if the sequence and teh structure are valid. The structure is optional.

    Parameters
    ----------
    seq : str
        Is a string of H and P of any length, if the sequence contains other letters, an error is raised. It is caps sensitive! 
        So H and P must be upper case only.
    struct : list, optional
        List containing the x and y coordinate of every monomer of the protein as INTEGER, 
        the number of element must be the same of the length of the sequence.\n
        The protein is assumed to move in a 2D lattice.\n
        EXAMPLE of length 6: [[0,0],[0,1],[1,1],[1,2],[1,3],[2,3]].\n
        By default (if nothing is inserted) a linear structure is assumed.

    Raises
    ------
    AssertionError
        If the sequence or structure are not valid or are not of compatible lengths.
    '''
    def __init__(self, seq : str, struct : list = None):
        if is_valid_sequence(seq): # check that the sequence is valid
            self.seq = seq
        else:
            raise AssertionError('The protein sequence is not valid')
            
        self.n = len(seq) # length of the sequence
        
        if struct is None: # linear structur if struct is not specified
            self.struct = []
            for i in range(self.n):
                self.struct.append([i,0])
            print('Linear initial structure assumed')
        else:
            try: # check that sequence has the right length
                assert len(struct) == self.n
            except:
                raise AssertionError('The lengths of the sequence and the structure are not the same')
            self.struct = struct
        
        if not is_valid_struct(self.struct): # check that the sequence is valid
            raise AssertionError('The structure is not a self avoid walk (SAW) or the distances between consecutive points is different from 1')
            
        pass




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
    Check if the protein sequence contains only H and/or P.

    Parameters
    ----------
    seq : str
        The protein sequence (caps sensitive).

    Returns
    -------
    bool
        True if the sequence is valid, False if is not.
    '''
    
    unique_c = set(seq) # set of all the letters present in the sequence
    if len(unique_c) == 2:
        if 'H' in unique_c and 'P' in unique_c:
            return True
    elif len(unique_c) == 1:
        if 'H' in unique_c or 'P' in unique_c:
            return True
    else:
        return False


def get_neig(coord : list) -> list:
    raise NotImplementedError()
    
    
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
    
    
    
    
    
    
    