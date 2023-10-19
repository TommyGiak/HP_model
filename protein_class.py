# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import random
import utils
import math
import numpy as np


class Protein():
    '''
    Protein class that contains all the information on the protein and the function to makes the system evolve.\n
    It takes as input the Configuration class present in the utils.py file.\n
    The protein sequence can be also coded in the 20 different amino acids (RNDQEHKSTACGILMFPWYV), in this case it will
    be automatically converted into the HP sequence considering the polar and hydrophobic amino acids,
    but it must be setted in the configuration file given as input.

    Parameters
    ----------
    config : utils.Configuration
        Configuration class already assigned using the selected input file.
    '''

    def __init__(self, config : utils.Configuration) -> None:
        
        if utils.is_valid_sequence(config.seq): # check that the sequence is valid (contains only HP)
            self.seq = config.seq
        else:
            self.seq = utils.hp_sequence_transform(config.seq) # if the sequence include the 20 different amino acids it will be coded in HP only
            print('\033[42mThe sequence was converted into binary configuration -> ', self.seq, '\033[0;0m')

        self.n = len(config.seq) # length of the sequence
        
        if not config.use_struct: # linear structur assumed if struct is not specified as input
            self.struct = utils.linear_struct(self.seq)
        else:
            try: # check that sequence has the right length
                assert len(config.struct) == self.n
            except:
                raise AssertionError('The lengths of the sequence and the structure are not the same')
            self.struct = config.struct
        
        if not utils.is_valid_struct(self.struct): # check that the sequence is valid
            raise AssertionError('The structure is not a self avoid walk (SAW) or the distances between consecutive points are different from 1')
        
        # parameters setting
        self.annealing = config.annealing
        self.T_in = config.T
        self.steps = config.folds
        self.gif = config.gif
        self.gif_struct = []

        self.min_en_struct = self.struct # variable to record the min energy structure (for now is the only structure)
        self.en_evo = [self.energy()] # list to keep track of the energy evolution
        self.T = [] # list to keep track of the temperature evolution
        self.counter = [] # counter of number of folding per step
        self.comp_evo = [self.compactness()] # list to keep track of the compactness evolution
        self.max_comp_struct = self.struct # variable to record the max compact structure (for now is the only structure)

        
    def evolution(self):
        '''
        Let the system evolving for a certain number of steps. 
        New structures are accepted following the Metropolis algorithm (this function basically apply the Metropolis alg).\n
        All the parameters are taken from the initial configuration.\n
        The energy evolution values, min energy structure and compactness conformations are saved.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        '''
        T = self.T_in
        self.T.append(T) # initial temperature
        m = -T/self.steps # angolar coefficient for the annealing

        for i in range(self.steps):
            utils.progress_bar(i+1,self.steps) # print the progress bar of the evolution

            if self.annealing and T > 0.002 : T = m*(i - self.steps) # temperature decrease linearly w.r.t. the steps, if annealing is True
            en = self.energy() # current protein energy
            init_str = self.struct # current protein structure
            self.struct = self.random_fold() # new structure is generated
            new_en = self.energy() # the energy of the new structure is computed
            
            if new_en > en: # if the new energy is higher to the previus one, the new structure is accepted following the Metropolis alg
                d_en = new_en - en # energy difference of the two states
                r = random.uniform(0, 1)
                p = math.exp(-d_en/T) # probability to accept the new structure
                if r > p:
                    self.struct = init_str # the new structure is not accepted (overwrite the initial structure)
                    
            if new_en < min(self.en_evo): # to save the min enrergy and structure
                self.min_en_struct = self.struct
            self.en_evo.append(new_en) # record the energy evolution

            self.comp_evo.append(self.compactness()) # save the compactness
            if self.comp_evo[-1] > max(self.comp_evo[:-1]):
                self.max_comp_struct = self.struct

            self.T.append(T) # record the T evolution

            if self.gif:
                if i%(int(self.steps/100)) == 0:
                    self.gif_struct.append(self.struct)
    
    
    def energy(self, e = 1.) -> float:
        '''
        Function to compute the energy of the protein structure. The binding energy can be changed.

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
        
        for i,seq in enumerate(self.seq):
            if seq == 'H':
                neig = self.get_neig_of(i) # string of neighbors of i-th monomer (exluding protein's backbone bonds)
                count_h += neig.count('H') 
        
        tot_en = -e*count_h/2 # total energy of the prot struct (/2 because each bond is counted twice)
        return tot_en
    
    
    def compactness(self) -> int:
        '''
        Function to compute the compactness of the structure.
        The compactness is the total number of neighbours of each monomer (backbone excluded).\n
        The return is twice the number of neighbours, but since the compactness is then normalized by its maximum value,
        is not necessary to divide by two

        Parameters
        ----------
        None

        Returns
        -------
        int :
            The total number of neighbours counted (doubled)
        '''
        count_neig = 0 
        
        for i,seq in enumerate(self.seq):
            count_neig += len(self.get_neig_of(i)) 
        
        return count_neig
                
    
    def get_neig_of(self, i : int) -> str:
        '''
        Function to see which are the neighbors of the i-th monomer of the protein sequence.
        The bounded monomer are not considered neighbors.\n
        The function return a string with the type of neighbour monomers: H/P.

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

        possible_neig = [[x-1,y],[x,y-1],[x+1,y],[x,y+1]]

        if i > 0: # if we are not in the starting monomer we remove the prevous monomer from the list
            possible_neig.remove(self.struct[i-1])

        if i < (self.n - 1): # if we are not in the ending monomer we remove the following monomer from the list
            possible_neig.remove(self.struct[i+1])

        for monomer in possible_neig: # check if the remaning points in the list contains some monomers (some neighbours)
            if monomer in self.struct:
                ind = self.struct.index(monomer) # get the position on the sequence of the neighbor
                neig += self.seq[ind] # get the H/P monomer

        return neig
        
    
    def random_fold(self) -> list:
        '''
        Randomly choose a monomer in the protein (exluding the first and the last) and fold the protein with a
        random method using the tail_fold function. If the structure generated is not valid
        the process is repited until a valid structure is found.

        Returns
        -------
        list
            The new rotein streucture randomly folded (valid).
        '''
        c = 0 # counter of the number of folding until a valid sequence is founded
        
        while True: # cycle valid until a valid structure is found
            index = random.randint(1, self.n-2) # select a random monomer where start the folding
            x, y = self.struct[index] 
            tail = self.struct[index:] # tail of the structure that will be folded

            # Excluding diagonal move is the sequence cannot support it (the previous and following monomer are aligned)
            distance_sur = utils.get_dist(self.struct[index-1],self.struct[index+1])
            diag_move = True if math.isclose(distance_sur,math.sqrt(2)) else False
            
            for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
                tail[i] = [mon[0]-x, mon[1]-y]
            previous = self.struct[index-1] # recording the prev monomer
            previous = [previous[0]-x, previous[1]-y]
                
            # choose a random method for the protein folding
            method = random.randint(1, 8) if diag_move else random.randint(1, 7) # exclude diagonal move if the conditions don't match
            tail = utils.tail_fold(struct=tail, method=method, previous=previous) # fold the tail with a random method
            
            for i,mon in enumerate(tail): # shifting the folded tail in the correct position
                tail[i] = [mon[0]+x, mon[1]+y]
            
            new_struct = self.struct[:index] # construction of the new structure generated
            for mon in tail: # pasting the new tail
                new_struct.append(mon)
                
            c += 1
                
            if utils.is_valid_struct(new_struct): # if the structure is valid and the cycle 
                break
            
        self.counter.append(c) # counter of the number of foldings

        return new_struct
