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
import numpy as np

class Protein():
    '''
    Protein class that contains: the sequence of the protein (Protein.seq), the length of the sequence 
    (Protein.n) and the structure of the protein (Protein.struc). \n
    When initialized it automatically check if the sequence and the structure are valid. The structure is optional. \n
    The protein sequence can be also coded in the 20 different amino acids (RNDQEHKSTACGILMFPWYV), in this case it will
    be automatically converted into the HP sequence considering the polar and hydrophobic amino acids.

    Parameters
    ----------
    seq : str
        Is a string of H and P of any length, if the sequence contains the 20 different amino-acids, it will be converted to HP sequences. It is caps sensitive! 
        So use upper case letters.
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
        
        if utils.is_valid_sequence(seq): # check that the sequence is valid (contains only HP)
            self.seq = seq
        else:
            self.seq = utils.hp_sequence_transform(seq) # if the sequence include the 20 different amino acids it will be coded in HP only
            print('\033[42mThe sequence was converted into binary configuration -> ', self.seq, '\033[0;0m')

        self.n = len(seq) # length of the sequence
        
        if struct is None: # linear structur assumed if struct is not specified as input
            self.struct = []
            for i in range(self.n):
                self.struct.append([i,0])
            print('\033[43mLinear initial structure assumed \033[0;0m')
        else:
            try: # check that sequence has the right length
                assert len(struct) == self.n
            except:
                raise AssertionError('The lengths of the sequence and the structure are not the same')
            self.struct = struct
        
        if not utils.is_valid_struct(self.struct): # check that the sequence is valid
            raise AssertionError('The structure is not a self avoid walk (SAW) or the distances between consecutive points are different from 1')
        
        self.min_en_struct = self.struct # variable to record the min energy structure (for now is the only structure)
        self.en_evo = [self.energy()] # list to keep track of the energy evolution
        self.T = [] # list to keep track of the temperature evolution
        self.counter = [] # counter of number of folding per step
        self.comp_evo = [self.compactness()] # list to keep track of the compactness evolution
        self.max_comp_struct = self.struct # variable to record the max compact structure (for now is the only structure)


    def view(self, save = True, tit = None):
        '''
        Function to plot the protein structure with matplotlib.
        Title can be optionally inserted.
        If save == True the plot will be also saved as pdf.
        '''
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
        ax.grid(alpha=0.2)
        if tit is not None:
            ax.set_title(tit)
        en = self.energy()
        comp = self.compactness()
        string = f'Energy: {en}'
        string_comp = f'Compactness: {comp/(max(self.comp_evo)+10e-15):.2f}'
        ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
        plt.show(block=False)
        if save:
            plt.savefig("data/prot_view.pdf", format="pdf", bbox_inches="tight")

        
    def evolution(self, annealing : bool = True, T : float = 1., steps : int = 10000):
        '''
        Let the system evolving for a certain number of steps. 
        New structures are accepted following the Metropolis algorithm (this function basically apply the Metropolis alg).\n
        The temperature can be controlled.\n
        The energy evolution values, min energy structure and compactness conformations are saved.

        Parameters
        ----------
        annealing : bool, optional
            If True, the temperature will slowly decrease as the steps increases.
            If False the temperature will remain constant during all the evolution.
            The default is True.
        T : float, optional
            Temperature of the enviroment. The default is 0.5.
        steps : int, optional
            Number protein folds. The default is 10000.

        Returns
        -------
        None.
        '''
        self.T.append(T) # initial temperature
        m = -1/steps # angolar coefficient for the annealing
        print('Evolution started')

        for i in range(steps):
            utils.progress_bar(i+1,steps) # print the progress bar of the evolution

            if annealing and T > 0.002 : T = m*(i - steps) # temperature decrease linearly w.r.t. the steps, if annealing is True
            en = self.energy() # current protein energy
            init_str = self.struct # current protein structure
            self.struct = self.random_fold() # new structure is generated
            new_en = self.energy() # the energy of the new structure is computed
            
            if new_en > en: # if the new energy is higher to the previus one, the new structure is accepted following the Metropolis alg
                d_en = new_en - en # energy difference of the two states
                r = random.uniform(0, 1)
                p = math.exp(-d_en/T) # probability to accept the new structure
                if r > p:
                    self.struct = init_str
                    
            if new_en < min(self.en_evo): # to save the min enrergy and structure
                self.min_en_struct = self.struct
            self.en_evo.append(new_en) # record the energy evolution

            self.comp_evo.append(self.compactness()) # save the compactness
            if self.comp_evo[-1] > max(self.comp_evo[:-1]):
                self.max_comp_struct = self.struct

            self.T.append(T) # record the T evolution
        print('Evolution ended')
    
    
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
        The compactness is the total number of neighbours of each monomer (backbone excluded).

        Parameters
        ----------
        None

        Returns
        -------
        None
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
            
            for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
                tail[i] = [mon[0]-x, mon[1]-y]
                
            tail = utils.tail_fold(tail) # fold the tail with a random method
            
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
    
    
    def view_min_en(self, save = True):
        '''
        Function to plot the protein structure founded whit less energy with matplotlib.
        The plot can be saved with save = True as pdf
        '''
        x = [] # x coordinates of the monomers (ordered)
        y = [] # y coordinates of the monomers (ordered)
        
        fig, ax = plt.subplots()
        for i in range(self.n):
            x.append(self.min_en_struct[i][0])
            y.append(self.min_en_struct[i][1])    
        ax.plot(x,y, alpha = 0.5)
        for i, coord in enumerate(self.struct):
            ax.scatter(x[i], y[i], marker='$'+self.seq[i]+'$', s=20, color = 'red')
        ax.set_xlim(min(x)-6,max(x)+6)
        ax.set_ylim(min(y)-6,max(y)+6)
        ax.grid(alpha=0.2)
        en = min(self.en_evo)
        comp = Protein(seq=self.seq,struct=self.min_en_struct).compactness()
        string = f'Energy: {en}'
        string_comp = f'Compactness: {comp/(max(self.comp_evo)+10e-15):.2f}'
        ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
        ax.set_title('Min energy structure')
        plt.show(block=False)
        if save:
            plt.savefig("data/min_energy_view.pdf", format="pdf", bbox_inches="tight")

             
    def view_max_comp(self, save = True):
        '''
        Function to plot the protein structure founded whit the max compactness.
        The plot can be saved with save = True as pdf

        '''
        x = [] # x coordinates of the monomers (ordered)
        y = [] # y coordinates of the monomers (ordered)
        
        fig, ax = plt.subplots()
        for i in range(self.n):
            x.append(self.max_comp_struct[i][0])
            y.append(self.max_comp_struct[i][1])    
        ax.plot(x,y, alpha = 0.5)
        for i, coord in enumerate(self.struct):
            ax.scatter(x[i], y[i], marker='$'+self.seq[i]+'$', s=20, color = 'red')
        ax.set_xlim(min(x)-6,max(x)+6)
        ax.set_ylim(min(y)-6,max(y)+6)
        ax.grid(alpha=0.2)
        ax.set_title('Max compactness structure')
        en = Protein(seq=self.seq,struct=self.max_comp_struct).energy()
        comp = max(self.comp_evo)
        string = f'Energy: {en}'
        string_comp = f'Compactness: {comp/(max(self.comp_evo)+10e-15):.2f}'
        ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
        plt.show(block=False)
        if save:
            plt.savefig("data/max_compactness_view.pdf", format="pdf", bbox_inches="tight")

        
    
    def plot_energy(self, avg : int = 1, save = True) -> None:
        '''
        plot the energy evolution of the system
        The plot can be saved with save = True as pdf

        Parameters
        ----------
        avg : int, optional
            The energy will be averaged every avg steps. The default is 1.

        Returns
        -------
        Plot
        '''
        en_evo = np.array(self.en_evo[1:].copy())
        T = np.array(self.T[1:])
        x = np.arange(0, len(en_evo), avg)
        try:
            en_evo = en_evo.reshape(-1,avg).mean(axis=1)
            T = T.reshape(-1,avg).mean(axis=1)
        except:
            print(f'Mean procedure skipped since the number of time steps is not a multiple of the avarage required: {avg}')
        fig, ax = plt.subplots()
        ax.set_title(f'Energy evolution of the system averaged by {avg} time steps')
        ax.set_xlabel('Time step')
        ax.set_ylabel('Energy', color = 'b')
        ax.tick_params(axis='y', labelcolor='b')
        ax.plot(x, en_evo, color ='b')
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color = 'r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T, color = 'r')
        fig.tight_layout()
        plt.show(block=False)
        if save:
            plt.savefig("data/energy_evolution.pdf", format="pdf", bbox_inches="tight")
        
        
    def plot_compactness(self, avg : int = 1, save = True) -> None:
        '''
        plot the compactness evolution of the system
        The plot can be saved with save = True as pdf

        Parameters
        ----------
        avg : int, optional
            The compactness will be averaged every avg steps. The default is 1.

        Returns
        -------
        Plot
        '''
        comp = np.array(self.comp_evo[1:].copy())
        comp = comp/max(comp)
        T = np.array(self.T[1:])
        x = np.arange(0, len(comp), avg)
        try:
            comp = comp.reshape(-1,avg).mean(axis=1)
            T = T.reshape(-1,avg).mean(axis=1)
        except:
            print(f'Mean procedure skipped since the number of time steps is not a multiple of the avarage required: {avg}')
        fig, ax = plt.subplots()
        ax.set_title(f'Compactness evolution of the system averaged by {avg} time steps')
        ax.set_xlabel('Time step')
        ax.set_ylabel('Compactness', color = 'b')
        ax.tick_params(axis='y', labelcolor='b')
        ax.plot(x, comp, color ='b')
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color = 'r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T, color = 'r')
        fig.tight_layout()
        plt.show(block=False)
        if save:
            plt.savefig("data/compactness_evolution.pdf", format="pdf", bbox_inches="tight")
