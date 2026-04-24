"""
@author: Tommaso Giacometti
"""

import copy
import os
import random
import sys

from math import isclose, sqrt

# Use insert(0, ...) so src/ is resolved before any other path entry,
# preventing standard-library or installed packages from shadowing local
# modules with the same name.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import protein as p

from config import Configuration
from geometry import get_distance, generate_linear_fold
from validation import is_valid_fold, is_valid_sequence
from transforms import tail_fold, diagonal_move
from utils import convert_to_hp

from fold_sampler import FoldSampler
from simulation import Simulation

config_path = os.path.join(os.path.dirname(__file__), 'config_test.yaml')
config = Configuration(config_path)

correct_structure = [[0, 0], [0, 1], [1, 1], [1, 2], [1, 3], [2, 3], [2, 2], [2, 1], [2, 0], [2, -1], [1, -1], [0, -1],
                     [-1, -1]]

seq = 'HPPHHPHPHPHHP'
seq1 = 'HHHPHPHPPPPPHPHPHPHHPHHPHPHHPHPPH'
seq2 = 'VFCNKASIRIPWTKLKTHPICLSLDKVIMEMSTCEEPRSPFAEK'
seq3 = 'VVEGISVSVNSIVIRIGAKAFNASFELSQLRIYSVNAHWEHGDLRFTRIQDPQRGEV'
seq4 = 'DLMSVVVFKITGVNGEIDIRGEDTEICLQVNQVTPDQLGNISLRHYLCNRPVGSDQKAVATVMPMKIQVSNTKINLKDDSPRSSTVSLEPAPVTVHIDHLVVERSDDGSFHIRDSHMLNTGNDLKENVKSDSV'
seq5 = 'LTSGKYDLKKQRSVTQATQTSPGVPWPSQSANFPEFSFDFTREQLMEENESLKQELAKAKMALAEAHLEKDALLHHIKKMTVE'

seq_invalid = 'ASDHLKGFDKJHDCVNB'

wrong_str_double_point = [[0, 0], [0, 0], [0, 1], [0, 2], [0, 3], [1, 3]]
wrong_str_skip_step_square = [[0, 0], [1, 1], [1, 2], [1, 3], [2, 3]]
wrong_str_skip_step_linear = [[0, 0], [1, 0], [3, 0], [2, 1], [2, 2], [2, 3]]
wrong_str_return_point = [[0, 0], [-1, 0], [-1, 1], [-1, 2], [-1, 1], [0, 1]]


def test_is_valid_struct_when_correct():
    assert is_valid_fold(correct_structure)


def test_is_valid_struct_when_wrong_double_point():
    assert not is_valid_fold(wrong_str_double_point)


def test_is_valid_struct_when_wrong_step_square():
    assert not is_valid_fold(wrong_str_skip_step_square)


def test_is_valid_struct_when_wrong_step_linear():
    assert not is_valid_fold(wrong_str_skip_step_linear)


def test_is_valid_struct_when_wrong_return_point():
    assert not is_valid_fold(wrong_str_return_point)


def test_is_valid_sequence_when_correct_only_H():
    assert is_valid_sequence('HHHHHHHHHHHHH')


def test_is_valid_sequence_when_correct_only_P():
    assert is_valid_sequence('PPPPPPPPP')


def test_is_valid_sequence_when_correct_both_HP():
    assert is_valid_sequence('PHPHPPPPHHHPHPHPH')


def test_is_valid_sequence_when_wrong_small_case_h():
    assert not is_valid_sequence('HHHhHH')


def test_is_valid_sequence_when_wrong_small_case_p():
    assert not is_valid_sequence('PPPpPPPPP')


def test_is_valid_sequence_when_wrong_not_HP_upper():
    assert not is_valid_sequence('PHPAPPPPHHHPHPHPH')
    assert not is_valid_sequence('HPHPLHHP')


def test_is_valid_sequence_when_wrong_not_HP_lower():
    assert not is_valid_sequence('PHPaPPPPHHHPHPHPH')
    assert not is_valid_sequence('HPHPlHHP')


def test_is_valid_sequence_when_wrong_too_short_2():
    assert not is_valid_sequence('HP')


def test_is_valid_sequence_when_wrong_too_short_1():
    assert not is_valid_sequence('P')


def test_get_dist_0():
    assert isclose(get_distance((1, 1), (1, 1)), 0)


def test_get_dist_1():
    assert isclose(get_distance((1, 1), (2, 1)), 1)


def test_get_dist_2():
    assert isclose(get_distance((1, 1), (3, 1)), 2)


def test_get_dist_sqrt_2():
    assert isclose(get_distance((0, 0), (1, 1)), sqrt(2))


def test_get_dist_2_sqrt_2():
    assert isclose(get_distance((-1, 1), (1, 3)), 2 * sqrt(2))


def test_energy_computation():
    prot1 = p.Protein(config)
    prot1.sequence = seq
    prot1.fold = correct_structure
    prot1.sequence_length = len(seq)
    assert isclose(prot1.get_energy(), -2.)


def test_get_neig_linear():
    prot = p.Protein(config)
    prot.sequence = 'HPHPHPHPPPPHHHHPPP'
    prot.fold = generate_linear_fold(prot.sequence)
    prot.sequence_length = len(prot.sequence)

    for i in range(prot.sequence_length):
        assert prot.get_neighbors(i) == ''


def test_get_neig_composite():
    prot = p.Protein(config)
    prot.sequence = seq
    prot.fold = correct_structure
    prot.sequence_length = len(seq)

    neig = ['H', '', 'P', 'H', '', '', 'H', 'P', '', '', '', 'H', '']

    for i in range(prot.sequence_length):
        assert prot.get_neighbors(i) == neig[i]


def test_random_fold_valid_struc_linear():
    random.seed(4326748)

    prot = p.Protein(config)
    prot.sequence = 'HPHPHPHPHPHHHHHPPHPHPHPPHHPPPPHHPP'
    prot.fold = generate_linear_fold(prot.sequence)
    prot.sequence_length = len(prot.sequence)

    sampler = FoldSampler()

    for _ in range(1000):
        prot.fold = sampler.sample(prot)
        assert is_valid_fold(prot.fold)


def test_random_fold_valid_struc_composite():
    random.seed(7694)

    prot = p.Protein(config)
    prot.sequence = seq
    prot.fold = correct_structure
    prot.sequence_length = len(seq)

    sampler = FoldSampler()

    for _ in range(1000):
        prot.fold = sampler.sample(prot)
        assert is_valid_fold(prot.fold)


def test_tail_fold_valid():
    tail = correct_structure[1:]
    x, y = tail[0]

    for i, mon in enumerate(tail):
        tail[i] = [mon[0] - x, mon[1] - y]

    previous = correct_structure[0]
    previous = [previous[0] - x, previous[1] - y]

    for method in range(1, 8):
        assert is_valid_fold(tail_fold(tail, method, previous))


def test_tail_fold_correct_length():
    l = len(correct_structure)

    for method in range(1, 8):
        l_new = len(tail_fold(correct_structure, method, [0, 0]))
        assert l == l_new


def test_diagonal_move_length():
    tail = correct_structure[1:]
    previous = correct_structure[0]

    assert len(tail) == len(diagonal_move(tail, previous))


def test_diagonal_move_equal_struct():
    new_struct = diagonal_move(correct_structure, [0, 0])
    assert new_struct[1:] == correct_structure[1:]


def test_hp_sequence_transform_letters_correct():
    assert set(convert_to_hp(seq2)) == {'H', 'P'}


def test_hp_sequence_transform_lenght_correct():
    assert len(convert_to_hp(seq2)) == len(seq2)


def test_evolution_minimize_energy():
    random.seed(config.seed)
    prot = p.Protein(config)
    prot.sequence = seq
    prot.fold = generate_linear_fold(seq)
    prot.sequence_length = len(seq)

    en1 = prot.get_energy()

    sim_config = copy.copy(config)
    sim_config.n_steps = 500

    sim = Simulation(prot, sim_config)
    sim.run()

    assert prot.get_energy() <= en1


def test_evolution_maximize_compactness():
    random.seed(config.seed)
    prot = p.Protein(config)
    prot.sequence = seq
    prot.fold = generate_linear_fold(seq)
    prot.sequence_length = len(seq)

    comp = prot.get_compactness()

    sim_config = copy.copy(config)
    sim_config.n_steps = 500

    sim = Simulation(prot, sim_config)
    sim.run()

    assert prot.get_compactness() >= comp
