# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 18:00:40 2023

@author: Tommaso Giacometti
"""

import protein_class as p
from math import isclose, sqrt


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
                      [[0,0],[1,0],[2,0],[2,1]],
                      [[1,0],[0,0],[0,1],[0,2]]]

seq = 'HPPHHPHPHPHHP'

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
        assert p.is_valid_struct(struct)
    
    
def test_is_valid_struct_when_wrong(structures = wrong_structures):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    
    for struct in structures:
        assert not p.is_valid_struct(struct)
    

def test_is_valid_sequence_when_correct():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of H/P\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    
    assert p.is_valid_sequence('HHHHHH')
    assert p.is_valid_sequence('PPPPPPPPP')
    assert p.is_valid_sequence('PHPHPPPPHHHPHPHPH')
    assert p.is_valid_sequence('HPHPHHHP')

        
def test_is_valid_sequence_when_wrong():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    
    assert not p.is_valid_sequence('HHHhHH')
    assert not p.is_valid_sequence('PPPpPPPPP')
    assert not p.is_valid_sequence('PHPAPPPPHHHPHPHPH')
    assert not p.is_valid_sequence('HPHPLHHP')
    assert not p.is_valid_sequence('HA')
    assert not p.is_valid_sequence('HP')
    assert not p.is_valid_sequence('P')
    assert not p.is_valid_sequence('PPPHHHl')        
        
        
def test_get_dist():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(p.get_dist((1,1), (1,1)), 0)
    assert isclose(p.get_dist((1,1), (2,1)), 1)
    assert isclose(p.get_dist((1,1), (3,1)), 2)
    assert isclose(p.get_dist((0,0), (1,1)), sqrt(2))
    assert isclose(p.get_dist((-1,1), (1,3)), 2*sqrt(2))  
    
        
def test_energy_computation():
    '''
    Test the energy computation of the protein structures for two structures defined above.
    '''
    prot1 = p.Protein(seq, correct_structures[0])
    prot2 = p.Protein(seq, correct_structures[1])
    assert isclose(prot1.energy(), -2.)
    assert isclose(prot2.energy(), -1.)
    
    
def test_get_neig():
    '''
    Test the get neighbors function using a hand made structure and a linear structure (which should not has neighbors)
    '''
    prot1 = p.Protein(seq, correct_structures[0])
    neig = ['H','','P','H','','','H','P','','','','H','']
    for i in range(prot1.n):
        assert prot1.get_neig_of(i) == neig[i]
    prot2 = p.Protein('HPHPHPHPPPPHHHHPPP')
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
    
    prot1 = p.Protein(seq,correct_structures[0])
    for i in range(n):
        prot1.struct = prot1.random_fold()
        assert p.is_valid_struct(prot1.struct)
        
    prot2 = p.Protein('HPHPHPHPHPHHHHHPPHPHPHPPHHPPPPHHPP')
    for i in range(n):
        prot2.struct = prot2.random_fold()
        assert p.is_valid_struct(prot2.struct)
    
    
def test_tail_fold_valid_struct():
    '''
    Test tail_fold gives valid sequences for each method

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to randomly fold the protein
    THEN: I expect a valid protein structure
    '''
    for struct in correct_structures:
        for i in range(7):
            assert p.tail_fold(struct,i+1)
            
            
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
            l_new = len(p.tail_fold(struct,i+1))
            assert l == l_new
            
        
      
    
    
    
    
        