# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 15:54:27 2023

@author: Tommaso Giacometti
"""
from math import sqrt, isclose
import matplotlib.pyplot as plt
import random

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
    
    
    def view(self):
        '''
        Function to plot the protein structure with matplotlib.
        '''
        assert is_valid_struct(self.struct)
        
        x = [] # x coordinates of the monomers (ordered)
        y = [] # y coordinates of the monomers (ordered)
        
        fig, ax = plt.subplots()
        for i in range(self.n):
            x.append(self.struct[i][0])
            y.append(self.struct[i][1])    
        ax.plot(x,y, alpha = 0.5)
        for i, coord in enumerate(self.struct):
            ax.scatter(x[i], y[i], marker='$'+self.seq[i]+'$', s=25, color = 'red')
        ax.set_xlim(-self.n,self.n)
        ax.set_ylim(-self.n,self.n)
        plt.show()
        
        pass
        
    
    def energy(self, e = 1.) -> float:
        '''
        Function to compute the energy of the protein srtucture. The binding energy can be inserted.

        Parameters
        ----------
        e : float, optional
            e represent the binding energy.\n
            The default is 1.

        Returns
        -------
        float
            The energy of the protein structure.
        '''
        count_h = 0 # counter of H-H neighbor pairs (exluding protein's backbone bonds)
        
        for i,mon in enumerate(self.struct):
            if self.seq[i] == 'H':
                neig = self.get_neig_of(i) # string of neighbors of i-th monomer (exluding protein's backbone bonds)
                count_h += neig.count('H') 
        
        tot_en = -e*count_h/2 # total energy of the prot struct (/2 because each bond is counted twice)
        return tot_en
                
    
    def get_neig_of(self, i : int) -> str:
        '''
        Function to see which are the neighbors of the i-th monomer of the protein sequence.
        The bounded monomer are not considered neighbors.\n
        The function return a string with the type of neighbor monomers: H/P.

        Parameters
        ----------
        i : int
            Monomer position in the sequence for which we want the neighbors.

        Returns
        -------
        str
            Neighbors type H/P.
        '''
        neig = '' # string to save the neighbors
        x,y = self.struct[i] # coordinates of the monomer
        
        if [x+1,y] in self.struct: # check each neighbor exist
            ind = self.struct.index([x+1,y]) # get the position on the sequence of the neighbor
            if ind != i+1 and ind != i-1: # exluding the backbone of the protein from neighbors
                neig += self.seq[ind] # get the neighbor type H/P
        
        if [x-1,y] in self.struct: # repeat the previus sequence for each neighbor in the sequence
            ind = self.struct.index([x-1,y])
            if ind != i+1 and ind != i-1:
                neig += self.seq[ind]
        
        if [x,y+1] in self.struct:
            ind = self.struct.index([x,y+1])
            if ind != i+1 and ind != i-1:
                neig += self.seq[ind]
        
        if [x,y-1] in self.struct:
            ind = self.struct.index([x,y-1])
            if ind != i+1 and ind != i-1:
                neig += self.seq[ind]
        
        return neig
    
    
    def random_fold(self) -> list:
        '''
        Randomly choos a monomer in the protein (exluding the first and the last) and fold the protein with a
        random method using the tail_fold function. If the structure generated is not valid
        the process is repited until a valid structure is found.

        Returns
        -------
        list
            The new rotein streucture randomly folded (valid).
        '''
        while True: # cycle valid until a valid structure is found
            index = random.randint(1, self.n-2) # select a random monomer where start the folding
            x, y = self.struct[index] 
            tail = self.struct[index:] # tail of the structure that wil be folded
            
            for i,mon in enumerate(tail): # shifting the tail with start on zero for the movement
                tail[i] = [mon[0]-x, mon[1]-y]
                
            tail = tail_fold(tail) # fold the tail with a random method
            
            for i,mon in enumerate(tail): # shifting the folded tail in the correct position
                tail[i] = [mon[0]+x, mon[1]+y]
            
            new_struct = self.struct[:index] # construction of the new structure generated
            for mon in tail: # pasting the new tail
                new_struct.append(mon)
                
            if is_valid_struct(new_struct): # if the structure is valid and the cycle 
                break
            
        return new_struct


def tail_fold(struct : list, method : int = None) -> list:
    '''
    Apply a rotation/inversion of symmetry at the sequence inserted.\n
    If no method is specified a random one is choosed.\n
    Are present 7 methods:
        1: 90° clockwise rotation
        2: 90° anticlockwise rotation
        3: 180° rotaion
        4: x-axis refletion
        5: y-axis reflection
        6: 1 and 3 quadrant bisector symmetry
        7: 2 and 4 quadrant bisector symmetry

    Parameters
    ----------
    struct : list
        Structure of the sequence for which each element is the x and y coordinates of the monomer.
    method : int, optional
        The method to apply to the structure. If no method is specified a random one is choosed.

    Returns
    -------
    list
        The structure transformed.
    '''
    if method == None: # if a specific method is not specified a random one is choosed
        method = random.randint(1, 7)
        
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
    return new_tail


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
    
    
    
if __name__ == '__main__':
    random.seed(0)
    a = [[0, 0],[0, 1],[1, 1],[1, 2],[1, 3],[2, 3],[2, 2],
         [2, 1],[2, 0],[2, -1],[1, -1],[0, -1],[-1, -1],
         [-1, -2],[-1, -3],[-1, -4],[-1, -5]]
    seq = 'HPPHHPHPHPHHPHHHH'
    prot = Protein(seq,a)
    prot.view()
    
    # for i in range(100):
    #     prot.struct = prot.random_fold()
    #     prot.view()
        
        
        
        