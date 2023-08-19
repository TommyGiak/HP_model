# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 18:00:40 2023

@author: Tommaso Giacometti
"""

from protein_class import is_valid_sequence, is_valid_struct, get_dist

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


def test_is_valid_struct_when_correct(correct_structures = correct_structures):
    '''
    Test the is_valid_struct when a correct structure is given, a list of correct structures are given in 
    the first part of the file.
    
    GIVEN: a correct structure
    WHEN: I want to verify if the structure is actually true with is_valid_struct
    THEN: I expect a True response from the function
    '''
    
    for struct in correct_structures:
        assert is_valid_struct(struct)
    
