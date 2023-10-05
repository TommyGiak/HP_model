# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""

import protein_class as p
import utils
import configparser
from math import isclose, sqrt
import hypothesis

configuration = configparser.ConfigParser()
configuration.read('config.txt')

config = utils.Configuration(configuration)


correct_structures = [[[0, 0],[0, 1],[1, 1],[1, 2],[1, 3],[2, 3],[2, 2],
                      [2, 1],[2, 0],[2, -1],[1, -1],[0, -1],[-1, -1]],                    
                      [[0, 0],[0, 1],[1, 1],[1, 2],[1, 3],[2, 3],[2, 2],
                      [2, 1],[2, 0],[2, -1],[2, -2],[2, -3],[2, -4]],
                      [[0,0],[0,1],[1,1],[1,2],[1,3],[2,3]],
                      [[0,0],[1,0],[1,1],[1,2],[1,3],[2,3]],
                      [[0,0],[0,-1],[1,-1],[1,-2],[1,-3],[2,-3]],
                      [[0,0],[1,0],[2,0],[2,1],[2,2],[2,3]],
                      [[0,0],[-1,0],[-1,1],[-1,2],[-1,3],[-2,3]],
                      [[0,0],[0,1],[1,1],[1,2],[1,3],[2,3]],
                      [[0,0],[1,0],[1,1],[1,2],[1,3]],
                      [[0,0],[0,-1],[1,-1],[1,-2],[1,-3]],
                      [[0,0],[1,0],[2,0],[2,1],[2,2]],
                      [[0,0],[-1,0],[-1,1]],
                      [[0,0],[1,0],[1,1],[1,2]],
                      [[0,0],[0,-1],[1,-1],[1,-2]],
                      [[0,0],[1,0],[2,0],[2,1]]]

seq = 'HPPHHPHPHPHHP'
seq1 = 'HHHPHPHPPPPPHPHPHPHHPHHPHPHHPHPPH'
seq2 = 'VFCNKASIRIPWTKLKTHPICLSLDKVIMEMSTCEEPRSPFAEK'
seq3 = 'VVEGISVSVNSIVIRIGAKAFNASFELSQLRIYSVNAHWEHGDLRFTRIQDPQRGEV'
seq4 = 'DLMSVVVFKITGVNGEIDIRGEDTEICLQVNQVTPDQLGNISLRHYLCNRPVGSDQKAVATVMPMKIQVSNTKINLKDDSPRSSTVSLEPAPVTVHIDHLVVERSDDGSFHIRDSHMLNTGNDLKENVKSDSV'
seq5 = 'LTSGKYDLKKQRSVTQATQTSPGVPWPSQSANFPEFSFDFTREQLMEENESLKQELAKAKMALAEAHLEKDALLHHIKKMTVE'
seq_invalid = 'ASDHLKGFDKJHDCVNB'


wrong_structures = [[[0,0],[0,0],[0,1],[0,2],[0,3],[1,3]],
                      [[0,0],[1,1],[1,2],[1,3],[2,3]],
                      [[0,0],[1,-1],[1,-2],[1,-3],[2,-3]],
                      [[0,0],[1,0],[3,0],[2,1],[2,2],[2,3]],
                      [[0,0],[-1,0],[-1,1],[-1,2],[-1,1],[0,1]],
                      [[0,0],[1,-1],[1,-2],[1,-3]],
                      [[0,0],[1,0],[3,0],[2,1],[2,3]],
                      [[0,0],[-1,0],[-1,1],[-1,2],[-1,1],[0,1]]]


def test_is_valid_struct_when_correct(structures = correct_structures):
    '''
    Test the is_valid_struct when a correct structure is given, a list of correct structures are given in 
    the first part of the file.
    
    GIVEN: a correct structure\n
    WHEN: I want to verify if the structure is actually see as true with is_valid_struct\n
    THEN: I expect a True response from the function
    '''
    
    for struct in structures:
        assert utils.is_valid_struct(struct)
    
    
def test_is_valid_struct_when_wrong(structures = wrong_structures):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    
    for struct in structures:
        assert not utils.is_valid_struct(struct)
    

def test_is_valid_sequence_when_correct():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of H/P\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    
    assert utils.is_valid_sequence('HHHHHH')
    assert utils.is_valid_sequence('PPPPPPPPP')
    assert utils.is_valid_sequence('PHPHPPPPHHHPHPHPH')
    assert utils.is_valid_sequence('HPHPHHHP')

        
def test_is_valid_sequence_when_wrong():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    
    assert not utils.is_valid_sequence('HHHhHH')
    assert not utils.is_valid_sequence('PPPpPPPPP')
    assert not utils.is_valid_sequence('PHPAPPPPHHHPHPHPH')
    assert not utils.is_valid_sequence('HPHPLHHP')
    assert not utils.is_valid_sequence('HA')
    assert not utils.is_valid_sequence('HP')
    assert not utils.is_valid_sequence('P')
    assert not utils.is_valid_sequence('PPPHHHl')        
        
        
def test_get_dist():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((1,1), (1,1)), 0)
    assert isclose(utils.get_dist((1,1), (2,1)), 1)
    assert isclose(utils.get_dist((1,1), (3,1)), 2)
    assert isclose(utils.get_dist((0,0), (1,1)), sqrt(2))
    assert isclose(utils.get_dist((-1,1), (1,3)), 2*sqrt(2))  
    
        
def test_energy_computation():
    '''
    Test the energy computation of the protein structures for two structures defined above.
    '''
    prot1, prot2 = p.Protein(config), p.Protein(config)
    prot1.seq = seq
    prot1.struct = correct_structures[0]
    prot2.seq = seq
    prot2.struct = correct_structures[1]
    prot1.n, prot2.n = len(seq), len(seq)
    assert isclose(prot1.energy(), -2.)
    assert isclose(prot2.energy(), -1.)
    
    
def test_get_neig():
    '''
    Test the get neighbors function using a hand made structure and a linear structure (which should not has neighbors)
    '''
    prot1 = p.Protein(config)
    prot1.seq = seq
    prot1.struct = correct_structures[0]
    prot1.n = len(seq)
    neig = ['H','','P','H','','','H','P','','','','H','']
    for i in range(prot1.n):
        assert prot1.get_neig_of(i) == neig[i]
    prot2 = p.Protein(config)
    prot2.seq = 'HPHPHPHPPPPHHHHPPP'
    prot2.struct = utils.linear_struct(prot2.seq)
    prot2.n = len(prot2.seq)
    for i in range(prot2.n):
        assert prot2.get_neig_of(i) == ''
        
    
def test_random_fold_valid_struc():
    '''
    Test that the random fald of the protein gives a valid structure

    GIVEN: a valid protein structure
    WHEN: I want to randomly fold the protein
    THEN: I expect a valid protein structure
    '''
    n = 1000
    
    prot1 = p.Protein(config)
    prot1.seq = seq
    prot1.struct = correct_structures[0]
    prot1.n = len(seq)
    for i in range(n):
        prot1.struct = prot1.random_fold()
        assert utils.is_valid_struct(prot1.struct)
        
    prot2 = p.Protein(config)
    prot2.seq = 'HPHPHPHPHPHHHHHPPHPHPHPPHHPPPPHHPP'
    prot2.struct = utils.linear_struct(prot2.seq)
    prot2.n = len(prot2.seq)
    for i in range(n):
        prot2.struct = prot2.random_fold()
        assert utils.is_valid_struct(prot2.struct)
    
    
def test_tail_fold_valid_struct():
    '''
    Test tail_fold gives valid sequences for each method

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to randomly fold the protein
    THEN: I expect a valid protein structure
    '''
    for struct in correct_structures:
        for i in range(7):
            assert utils.tail_fold(struct,i+1)
            
            
def test_tail_fold_correct_length():
    '''
    Test that tail_fold does not change the sequence length

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to randomly fold the protein
    THEN: I expect the structure length unchanged
    '''
    for struct in correct_structures:
        for i in range(7):
            l = len(struct)
            l_new = len(utils.tail_fold(struct,i+1))
            assert l == l_new
  
            
def test_diagonal_move_length():
    '''
    Test the constant length of the protein when diagonal_move is applied

    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: the structure length should not change
    '''
    for struct in correct_structures:
        l = len(struct)
        
        assert l == len(utils.diagonal_move(struct))
      

def test_diagonal_move_equal_struct():
    '''
    Test that diagolan_move let unchanged the structure from the second monomer.
    
    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: the structure starting from the second monomer should not change
    '''
    for struct in correct_structures:
        new_struct = utils.diagonal_move(struct)
        assert new_struct[1:] == struct[1:]
        
    
def test_diagonal_move_first_mon_move():
    '''
    Test that diagolan_move moves the first monomer near the second one.
    
    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: I expect the distance between first and second monomer equal to one
    '''
    for struct in correct_structures:
        struct = utils.diagonal_move(struct)
        d = utils.get_dist(struct[0], struct[1])
        assert isclose(d, 1)
    
    
def test_hp_sequence_transform_letters():
    '''
    Test that hp_sequence_transform return a str with only H and P, if an invalid sequence is passed is should fail.
    
    GIVEN: a list of random sequences and an invalid sequence
    WHEN: I want to convert a complete amino-acid sequence into a sequence with only H and P
    THEN: I expect a string as return containing only H and P
    '''
    assert set(utils.hp_sequence_transform(seq1)) == {'H', 'P'}
    assert set(utils.hp_sequence_transform(seq2)) == {'H', 'P'}
    assert set(utils.hp_sequence_transform(seq3)) == {'H', 'P'}
    assert set(utils.hp_sequence_transform(seq4)) == {'H', 'P'}
    assert set(utils.hp_sequence_transform(seq5)) == {'H', 'P'}
    try:
        s = utils.hp_sequence_transform(seq_invalid)
        raise ValueError('An invalid sequence is passed')
    except:
        pass


def test_hp_sequence_transform_lenght():
    '''
    Test that hp_sequence_transform conserve the length of the sequence.
    
    GIVEN: a list of random sequences
    WHEN: I want to convert a complete amino-acid sequence into a sequence with only H and P
    THEN: I expect that the converted string has the same lenght than before
    '''
    assert len(utils.hp_sequence_transform(seq1)) == len(seq1)
    assert len(utils.hp_sequence_transform(seq2)) == len(seq2)
    assert len(utils.hp_sequence_transform(seq3)) == len(seq3)
    assert len(utils.hp_sequence_transform(seq4)) == len(seq4)
    assert len(utils.hp_sequence_transform(seq5)) == len(seq5)


def test_evolution_minimize_energy():
    '''
    Test that the evolution of the protein takes to a minimization of the energy. 
    The energy should never be grater than zero
    
    GIVEN: a list of random sequences
    WHEN: I evolve the system for a certain number of steps
    THEN: I expect that the energy of the protein decrease (or remain equal to zero)
    '''
    # definition of the protein to test
    prot1, prot2, prot3 = p.Protein(config), p.Protein(config), p.Protein(config)
    prot1.seq, prot2.seq, prot3.seq = seq, seq1, seq2
    prot1.struct = utils.linear_struct(prot1.seq)
    prot2.struct = utils.linear_struct(prot2.seq)
    prot3.struct = utils.linear_struct(prot3.seq)
    prot1.n, prot2.n, prot3.n = len(seq), len(seq1), len(seq2)

    prot1.steps, prot2.steps, prot3.steps = 500, 500, 500

    en1, en2, en3 = prot1.energy(), prot2.energy(), prot3.energy()

    # energy shoud be zero for the linear structure
    assert isclose(en1,0.)
    assert isclose(en2,0.)
    assert isclose(en3,0.)

    # evolitions
    prot1.evolution()
    prot2.evolution()
    prot3.evolution()

    # asserts for the energy minimization after the evolution (energy shoul not be grater than zero)
    assert prot1.energy() <= en1
    assert prot2.energy() <= en2
    assert prot3.energy() <= en3


def test_evolution_maximize_compactness():
    '''
    Test that the evolution of the protein takes does not takes the compactness below zero. 
    
    GIVEN: a list of random sequences
    WHEN: I evolve the system for a certain number of steps
    THEN: I expect that the compactness of the protein never goes below zero
    '''
    # definition of the protein to test
    prot1, prot2, prot3 = p.Protein(config), p.Protein(config), p.Protein(config)
    prot1.seq, prot2.seq, prot3.seq = seq, seq1, seq2
    prot1.struct = utils.linear_struct(prot1.seq)
    prot2.struct = utils.linear_struct(prot2.seq)
    prot3.struct = utils.linear_struct(prot3.seq)
    prot1.n, prot2.n, prot3.n = len(seq), len(seq1), len(seq2)

    prot1.steps, prot2.steps, prot3.steps = 500, 500, 500

    comp1, comp2, comp3 = prot1.compactness(), prot2.compactness(), prot3.compactness()

    # energy shoud be zero for the linear structure
    assert isclose(comp1,0.)
    assert isclose(comp2,0.)
    assert isclose(comp3,0.)

    # evolitions
    prot1.evolution()
    prot2.evolution()
    prot3.evolution()

    # asserts for the energy minimization after the evolution (energy shoul not be grater than zero)
    assert prot1.compactness() >= comp1
    assert prot2.compactness() >= comp2
    assert prot3.compactness() >= comp3
