# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""

import protein_class as p
import utils
import configparser
from math import isclose, sqrt
import random
import hypothesis

configuration = configparser.ConfigParser()
configuration.read('config_test.txt')

config = utils.Configuration(configuration)


correct_structure = [[0,0],[0,1],[1,1],[1,2],[1,3],[2,3],[2,2],[2,1],[2,0],[2,-1],[1,-1],[0,-1],[-1,-1]]

seq = 'HPPHHPHPHPHHP'
seq1 = 'HHHPHPHPPPPPHPHPHPHHPHHPHPHHPHPPH'
seq2 = 'VFCNKASIRIPWTKLKTHPICLSLDKVIMEMSTCEEPRSPFAEK'
seq3 = 'VVEGISVSVNSIVIRIGAKAFNASFELSQLRIYSVNAHWEHGDLRFTRIQDPQRGEV'
seq4 = 'DLMSVVVFKITGVNGEIDIRGEDTEICLQVNQVTPDQLGNISLRHYLCNRPVGSDQKAVATVMPMKIQVSNTKINLKDDSPRSSTVSLEPAPVTVHIDHLVVERSDDGSFHIRDSHMLNTGNDLKENVKSDSV'
seq5 = 'LTSGKYDLKKQRSVTQATQTSPGVPWPSQSANFPEFSFDFTREQLMEENESLKQELAKAKMALAEAHLEKDALLHHIKKMTVE'

seq_invalid = 'ASDHLKGFDKJHDCVNB' # contains letters which don't represent any amino-acids

wrong_str_double_point = [[0,0],[0,0],[0,1],[0,2],[0,3],[1,3]]
wrong_str_skip_step_square = [[0,0],[1,1],[1,2],[1,3],[2,3]]
wrong_str_skip_step_linear = [[0,0],[1,0],[3,0],[2,1],[2,2],[2,3]]
wrong_str_return_point = [[0,0],[-1,0],[-1,1],[-1,2],[-1,1],[0,1]]


def test_is_valid_struct_when_correct(structure = correct_structure):
    '''
    Test the is_valid_struct when a correct structure is given, a list of correct structures are given in 
    the first part of the file.
    
    GIVEN: a correct structure\n
    WHEN: I want to verify if the structure is actually see as true with is_valid_struct\n
    THEN: I expect a True response from the function
    '''
    
    assert utils.is_valid_struct(structure)
    
    
def test_is_valid_struct_when_wrong_double_point(structure = wrong_str_double_point):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure with a same point repited twice\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    assert not utils.is_valid_struct(structure)
    
def test_is_valid_struct_when_wrong_step_square(structure = wrong_str_skip_step_square):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure with a point is skipped -> the distance between two consecutive points is sqrt(2)\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    assert not utils.is_valid_struct(structure)
    
def test_is_valid_struct_when_wrong_step_linear(structure = wrong_str_skip_step_linear):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure with a point is skipped -> the distance between two consecutive points is 2\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    assert not utils.is_valid_struct(structure)

def test_is_valid_struct_when_wrong_return_point(structure = wrong_str_return_point):
    '''
    Test the is_valid_struct when a wrong structure is given, a list of wrong structures are given in 
    the first part of the file.
    
    GIVEN: a wrong structure with sequence that return on the same point (overlap)\n
    WHEN: I want to verify if the structure is actually see as wrong with is_valid_struct\n
    THEN: I expect a False response from the function
    '''
    assert not utils.is_valid_struct(structure)


def test_is_valid_sequence_when_correct_only_H():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of H\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    
    assert utils.is_valid_sequence('HHHHHHHHHHHHH')

def test_is_valid_sequence_when_correct_only_P():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of P\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    assert utils.is_valid_sequence('PPPPPPPPP')

def test_is_valid_sequence_when_correct_both_HP():
    '''
    Test the is_valid_sequence when the sequence is correct.
    A bounch of test sequence are given, including cases with only H or only P.

    GIVEN: a correct protein sequence of H/P\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return True
    '''
    assert utils.is_valid_sequence('PHPHPPPPHHHPHPHPH')

        
def test_is_valid_sequence_when_wrong_small_case_h():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a smaller case letter h\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('HHHhHH')

def test_is_valid_sequence_when_wrong_small_case_p():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a smaller case letter p\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('PPPpPPPPP')
  
def test_is_valid_sequence_when_wrong_not_HP_upper():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a letter different from H/P upper case\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('PHPAPPPPHHHPHPHPH')
    assert not utils.is_valid_sequence('HPHPLHHP')
     
def test_is_valid_sequence_when_wrong_not_HP_upper():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a letter different from H/P lower case\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('PHPaPPPPHHHPHPHPH')
    assert not utils.is_valid_sequence('HPHPlHHP')

def test_is_valid_sequence_when_wrong_too_short_2():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a too short sequence (2)\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('HP')

def test_is_valid_sequence_when_wrong_too_short_1():
    '''
    Test the is_valid_sequence when the sequence is wrong.
    A bounch of test sequence are given, including cases with lowercase h and p.

    GIVEN: a wrong protein sequence with a too short sequence (1)\n
    WHEN: is_valid_sequence is apply\n
    THEN: I expect that the function return False
    '''
    assert not utils.is_valid_sequence('P')


def test_get_dist_0():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice distance zero\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((1,1), (1,1)), 0)
    
def test_get_dist_1():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice distance one\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((1,1), (2,1)), 1)   
    
def test_get_dist_2():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice distance 2\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((1,1), (3,1)), 2)

def test_get_dist_sqrt_2():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice distance sqrt(2)\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((0,0), (1,1)), sqrt(2))

def test_get_dist_2_sqrt_2():
    ''' 
    Test get_dist that computes correctly the euclidean distance
    
    GIVEN: different points in the lattice distance 2*sqrt(2)\n
    WHEN: I want to compute the euclidean distance\n
    THEN: I expect the correct euclidean distance
    '''
    assert isclose(utils.get_dist((-1,1), (1,3)), 2*sqrt(2))  


def test_energy_computation():
    '''
    Test the energy computation of the protein structures for two structures defined above.

    GIVEN: a specific porotein structure\n
    WHEN: I want to compute the energy\n
    THEN: I expect that tha energy is exactly -2
    '''
    prot1 = p.Protein(config)
    prot1.seq = seq
    prot1.struct = correct_structure
    prot1.n = len(seq)
    assert isclose(prot1.energy(), -2.)

    
def test_get_neig_linear():
    '''
    Test the get neighbors function using a hand made structure and a linear structure (which should not has neighbors)

    GIVEN: a linear structure for a protein\n
    WHEN: I want to check the neighbours of each monomer\n
    THEN: I expect no neighbours for the monomers
    '''
    prot = p.Protein(config)
    prot.seq = 'HPHPHPHPPPPHHHHPPP'
    prot.struct = utils.linear_struct(prot.seq)
    prot.n = len(prot.seq)
    for i in range(prot.n):
        assert prot.get_neig_of(i) == ''

def test_get_neig_composite():
    '''
    Test the get neighbors function using a hand made structure and a linear structure (which should not has neighbors)

    GIVEN: a composite structure for a protein\n
    WHEN: I want to check the neighbours of each monomer\n
    THEN: I expect specific neighbours for each monomer
    '''
    prot = p.Protein(config)
    prot.seq = seq
    prot.struct = correct_structure
    prot.n = len(seq)
    neig = ['H','','P','H','','','H','P','','','','H','']
    for i in range(prot.n):
        assert prot.get_neig_of(i) == neig[i]
        
    
def test_random_fold_valid_struc_linear():
    '''
    Test that the random fald of the protein gives a valid structure

    GIVEN: a valid linear protein structure
    WHEN: I want to randomly fold the protein
    THEN: I expect a valid protein structure
    '''
    random.seed(4326748)
    n = 1000
    prot = p.Protein(config)
    prot.seq = 'HPHPHPHPHPHHHHHPPHPHPHPPHHPPPPHHPP'
    prot.struct = utils.linear_struct(prot.seq)
    prot.n = len(prot.seq)
    for i in range(n):
        prot.struct = prot.random_fold()
        assert utils.is_valid_struct(prot.struct)

def test_random_fold_valid_struc_composite():
    '''
    Test that the random fald of the protein gives a valid structure

    GIVEN: a valid composite protein structure
    WHEN: I want to randomly fold the protein
    THEN: I expect a valid protein structure
    '''
    random.seed(7694)
    n = 1000
    prot = p.Protein(config)
    prot.seq = seq
    prot.struct = correct_structure
    prot.n = len(seq)
    for i in range(n):
        prot.struct = prot.random_fold()
        assert utils.is_valid_struct(prot.struct)

    
def test_tail_fold_valid_struct_1():
    '''
    Test tail_fold gives valid sequences for 90 deg rotation clock

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 1, previous))

def test_tail_fold_valid_struct_2():
    '''
    Test tail_fold gives valid sequences for 90 deg rotation anticlock

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
        
    assert utils.is_valid_struct(utils.tail_fold(tail, 2, previous))

def test_tail_fold_valid_struct_3():
    '''
    Test tail_fold gives valid sequences for 180 deg rotation

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 3, previous))
            
def test_tail_fold_valid_struct_4():
    '''
    Test tail_fold gives valid sequences for x-reflection

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 4, previous))

def test_tail_fold_valid_struct_5():
    '''
    Test tail_fold gives valid sequences for y-reflection

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 5, previous))

def test_tail_fold_valid_struct_6():
    '''
    Test tail_fold gives valid sequences for bisec symmetry

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 6, previous))

def test_tail_fold_valid_struct_7():
    '''
    Test tail_fold gives valid sequences for -bisec symmetry

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect a valid protein structure
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]
    
    assert utils.is_valid_struct(utils.tail_fold(tail, 7, previous))


def test_tail_fold_correct_length_1():
    '''
    Test that tail_fold does not change the sequence length for 90 deg rotation clock

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,1,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_2():
    '''
    Test that tail_fold does not change the sequence length for 90 deg rotation anticlock

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,2,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_3():
    '''
    Test that tail_fold does not change the sequence length for 180 deg rotation

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,3,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_4():
    '''
    Test that tail_fold does not change the sequence length for x-reflection

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,4,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_5():
    '''
    Test that tail_fold does not change the sequence length for y-reflection

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,5,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_6():
    '''
    Test that tail_fold does not change the sequence length for bisec symmetry

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,6,[0,0]))
    assert l == l_new

def test_tail_fold_correct_length_7():
    '''
    Test that tail_fold does not change the sequence length for -bisec symmetry

    GIVEN: a valid protein structure and a specific method to use
    WHEN: I want to fold the protein
    THEN: I expect the structure length unchanged
    '''
    l = len(correct_structure)
    l_new = len(utils.tail_fold(correct_structure,7,[0,0]))
    assert l == l_new
  
            
def test_diagonal_move_length():
    '''
    Test the constant length of the protein when diagonal_move is applied

    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: the structure length should not change
    '''
    tail = correct_structure[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = correct_structure[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]

    l = len(tail)
    assert l == len(utils.diagonal_move(tail,previous))
    

def test_diagonal_move_equal_struct():
    '''
    Test that diagolan_move lets unchanged the structure from the second monomer.
    
    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: the structure starting from the second monomer should not change
    '''
    new_struct = utils.diagonal_move(correct_structure,[0,0])
    assert new_struct[1:] == correct_structure[1:]
        
    
def test_diagonal_move_first_mon_move():
    '''
    Test that diagolan_move moves the first monomer near the second one.
    
    GIVEN: a structure starting in [0,0]
    WHEN: I want to move the first monomer on a diagonal to fold the protein
    THEN: I expect the distance between first and second monomer equal to one
    '''
    struct = [[0,0],[1,0],[1,1],[2,1],[2,2],[3,2],[3,3]]
    tail = struct[1:]
    x,y = tail[0]
    for i,mon in enumerate(tail): # shifting the tail start in [0,0] for the folding
            tail[i] = [mon[0]-x, mon[1]-y]
    previous = struct[0] # recording the prev monomer
    previous = [previous[0]-x, previous[1]-y]

    tail = utils.diagonal_move(tail,previous)
    d = utils.get_dist(tail[0], tail[1])
    assert isclose(d, 1)
    
    
def test_hp_sequence_transform_letters_correct():
    '''
    Test that hp_sequence_transform return a str with only H and P.
    
    GIVEN: a list of random sequences and an invalid sequence
    WHEN: I want to convert a complete amino-acid sequence into a sequence with only H and P
    THEN: I expect a string as return containing only H and P
    '''
    assert set(utils.hp_sequence_transform(seq2)) == {'H', 'P'}

def test_hp_sequence_transform_letters_wrong():
    '''
    Test that hp_sequence_transform fails with a invalid sequence
    
    GIVEN: an invalid sequence
    WHEN: I want to convert a complete amino-acid sequence into a sequence with only H and P
    THEN: I and error arised in the function
    '''
    try:
        s = utils.hp_sequence_transform(seq_invalid)
        raise ValueError('An invalid sequence is passed')
    except:
        pass


def test_hp_sequence_transform_lenght_correct():
    '''
    Test that hp_sequence_transform conserve the length of the sequence.
    
    GIVEN: a list of random sequences
    WHEN: I want to convert a complete amino-acid sequence into a sequence with only H and P
    THEN: I expect that the converted string has the same lenght than before
    '''
    assert len(utils.hp_sequence_transform(seq2)) == len(seq2)


def test_evolution_minimize_energy():
    '''
    Test that the evolution of the protein takes to a minimization of the energy. 
    The energy should never be grater than zero
    
    GIVEN: a list of random sequences
    WHEN: I evolve the system for a certain number of steps
    THEN: I expect that the energy of the protein decrease (or remain equal to zero)
    '''
    # definition of the protein to test
    prot1 = p.Protein(config)
    prot1.seq = seq
    prot1.struct = utils.linear_struct(prot1.seq)
    prot1.n = len(seq)
    prot1.steps = 500
    en1 = prot1.energy()
    # energy shoud be zero for the linear structure
    assert isclose(en1,0.)
    # evolition
    prot1.evolution()
    # asserts for the energy minimization after the evolution (energy shoul not be grater than zero)
    assert prot1.energy() <= en1


def test_evolution_maximize_compactness():
    '''
    Test that the evolution of the protein takes does not takes the compactness below zero. 
    
    GIVEN: a list of random sequences
    WHEN: I evolve the system for a certain number of steps
    THEN: I expect that the compactness of the protein never goes below zero
    '''
    # definition of the protein to test
    prot1 = p.Protein(config)
    prot1.seq = seq
    prot1.struct = utils.linear_struct(prot1.seq)
    prot1.n = len(seq)
    prot1.steps = 500
    comp1 = prot1.compactness()
    # energy shoud be zero for the linear structure
    assert isclose(comp1,0.)
    # evolitions
    prot1.evolution()
    # asserts for the energy minimization after the evolution (energy shoul not be grater than zero)
    assert prot1.compactness() >= comp1
