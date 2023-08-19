# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 18:00:40 2023

@author: Tommaso Giacometti
"""

from protein_class import is_valid_sequence, is_valid_struct, get_dist
from math import isclose, sqrt


correct_structures = [[[0,0],[0,1],[1,1],[1,2],[1,3],[2,3]],
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
        assert is_valid_struct(struct)
    
    
def test_is_valid_struct_when_wrong(structures = wrong_structures):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    
    for struct in structures:
        assert not is_valid_struct(struct)
    

def test_is_valid_sequence_when_correct():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of H/P\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    
    assert is_valid_sequence('HHHHHH')
    assert is_valid_sequence('PPPPPPPPP')
    assert is_valid_sequence('PHPHPPPPHHHPHPHPH')
    assert is_valid_sequence('HPHPHHHP')
    assert is_valid_sequence('H')
    assert is_valid_sequence('P')

        
def test_is_valid_sequence_when_wrong():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    
    assert not is_valid_sequence('HHHhHH')
    assert not is_valid_sequence('PPPpPPPPP')
    assert not is_valid_sequence('PHPAPPPPHHHPHPHPH')
    assert not is_valid_sequence('HPHPLHHP')
    assert not is_valid_sequence('HA')
    assert not is_valid_sequence('PPPHHHl')        
        
        
def test_get_dist():
    '''
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(get_dist((1,1), (1,1)), 0)
    assert isclose(get_dist((1,1), (2,1)), 1)
    assert isclose(get_dist((1,1), (3,1)), 2)
    assert isclose(get_dist((0,0), (1,1)), sqrt(2))
    assert isclose(get_dist((-1,1), (1,3)), 2*sqrt(2))  
        
        