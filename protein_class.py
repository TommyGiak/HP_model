# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 15:54:27 2023

@author: Tommaso Giacometti
"""
import matplotlib.pyplot as plt
import random
import utils
import math
import time

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
        
        if utils.is_valid_sequence(seq): # check that the sequence is valid
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
        
        if not utils.is_valid_struct(self.struct): # check that the sequence is valid
            raise AssertionError('The structure is not a self avoid walk (SAW) or the distances between consecutive points is different from 1')
        
        self.min_en_struct = self.struct # current min energy structure
        self.en_evo = [self.energy()] # to keep track of the energy evolution
        
    
    def view(self):
        '''
        Function to plot the protein structure with matplotlib.
        '''
        assert utils.is_valid_struct(self.struct)
        
        x = [] # x coordinates of the monomers (ordered)
        y = [] # y coordinates of the monomers (ordered)
        
        fig, ax = plt.subplots()
        for i in range(self.n):
            x.append(self.struct[i][0])
            y.append(self.struct[i][1])    
        ax.plot(x,y, alpha = 0.5)
        for i, coord in enumerate(self.struct):
            ax.scatter(x[i], y[i], marker='$'+self.seq[i]+'$', s=20, color = 'red')
        ax.set_xlim(min(x)-6,max(x)+6)
        ax.set_ylim(min(y)-6,max(y)+6)
        en = self.energy()
        string = f'Energy: {en}'
        ax.text(0.01,0.99,string, ha='left', va='top', transform=ax.transAxes)
        plt.show()

        
    def evolution(self, T : float = 0.5, steps : int = 10000):
        '''
        Let the system evolving for a certain number of steps of protein folds. 
        The new structures are accepted following the Metropolis algorithm.\n
        The temperature can be controlled.\n
        The energy evolution values and min energy structure conformation are saved.

        Parameters
        ----------
        T : float, optional
            Temperature of the enviroment. The default is 0.5.
        steps : int, optional
            Number protein folds. The default is 10000.

        Returns
        -------
        None.
        '''
        for i in range(steps):
            en = self.energy() # current protein energy
            init_str = self.struct # current protein structure
            self.struct = self.random_fold() # the new structure already replace the previus one
            new_en = self.energy() # the new energy is computed
            
            if new_en > en: # if the new energy higher the previus, the old structure remain following the Metropolis alg.
                d_en = new_en-en
                r = random.uniform(0, 1)
                p = math.exp(-d_en/T)
                if r > p:
                    self.struct = init_str
                    
            if new_en < min(self.en_evo): # to save the min enrergy and structure
                self.min_en_struct = self.struct
            self.en_evo.append(new_en) # record the energy evolution
    
    
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
                
            if utils.is_valid_struct(new_struct): # if the structure is valid and the cycle 
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
        8: movement on a digaonal of a random monomer

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
        method = random.randint(1, 8)
        
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
        new_tail = utils.diagonal_move(struct)
                
    return new_tail
    
    
if __name__ == '__main__':
    random.seed(0)
    a = [[0, 0],[0, 1],[1, 1],[1, 2],[1, 3],[2, 3],[2, 2],
         [2, 1],[2, 0],[2, -1],[1, -1],[0, -1],[-1, -1],
         [-1, -2],[-1, -3],[-1, -4],[-1, -5]]
    seq = 'PHPPPHPPHHPHPPHPPPHPHPHPPPHPPPHHHPHPHPH'
    prot = Protein(seq)
    prot.view()
    
    start = time.time()
    
    prot.evolution(T=0.5,steps=10000)
    prot.view()

    print(f'It took {time.time()-start:.3f} seconds')
        
        
        