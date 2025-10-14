"""
@author: Tommaso Giacometti
"""
import json

from math import isclose


def is_valid_struct(struct: list[list[float]]) -> bool:
    """
    Check if the given protein structure is valid:
    - Self-avoiding walk (no overlapping monomers)
    - Distance between consecutive monomers is 1

    Parameters
    ----------
    struct : list[list[float]]
        Structure of the protein containing x and y coordinates for each monomer.

    Returns
    -------
    bool
        True if the structure is valid, False otherwise.
    """
    seen = set()  # positions already visited

    for i, pos in enumerate(struct):
        pos_tuple = tuple(pos)  # convert to tuple to store in set
        if pos_tuple in seen:
            return False
        seen.add(pos_tuple)

        # check distance to next monomer
        if i < len(struct) - 1:
            if not isclose(get_dist(pos, struct[i + 1]), 1):
                return False

    return True


def is_valid_sequence(sequence: str) -> bool:
    """
    Check if the protein sequence contains only 'H' and/or 'P' and has a minimum length of 3.

    Parameters
    ----------
    sequence : str
        The protein sequence (case-sensitive).

    Returns
    -------
    bool
        True if the sequence is valid, False otherwise.
    """
    if len(sequence) < 3:
        print('The sequence is too short. It must be at least 3.')
        return False

    allowed_chars = {'H', 'P'}
    sequence_chars = set(sequence)

    # Sequence is valid if all characters are in allowed set
    return sequence_chars.issubset(allowed_chars)


def linear_struct(seq: str) -> list[list[int]]:
    """
    Create a linear structure corresponding to the input sequence.

    Parameters
    ----------
    seq : str
        Sequence for which to create a linear structure.

    Returns
    -------
    list[list[int]]
        Linear structure, where each monomer is positioned on the x-axis.
    """
    struct = [[i, 0] for i in range(len(seq))]
    print("Linear initial structure assumed")
    return struct


from math import sqrt


def get_dist(coord1: list[float], coord2: list[float]) -> float:
    """
    Compute the Euclidean distance between two points in the lattice.

    Parameters
    ----------
    coord1 : list[float]
        [x, y] coordinates of the first monomer.
    coord2 : list[float]
        [x, y] coordinates of the second monomer.

    Returns
    -------
    float
        Euclidean distance.
    """
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    return sqrt(dx ** 2 + dy ** 2)


import random


def diagonal_move(struct: list[list[int]], previous: list[int]) -> list[list[int]]:
    """
    Move the first monomer along a diagonal based on the previous and following monomers.
    Assumes the protein structure starts at [0,0], so coordinates of surrounding monomers are 0 or 1.
    Called only if a diagonal move is possible (surrounding monomers not aligned).

    Parameters
    ----------
    struct : list[list[int]]
        Protein structure starting at [0,0].
    previous : list[int]
        Coordinates [x, y] of the previous monomer, shifted such that the first monomer of struct is [0,0].

    Returns
    -------
    list[list[int]]
        Updated structure with the first monomer moved.
    """
    # Randomly choose one of the two possible diagonal directions
    case = random.randint(0, 1)

    x_prev, y_prev = previous
    x_next, y_next = struct[1]

    # Move the first monomer along the diagonal
    struct[0] = [x_prev + x_next, y_prev + y_next]

    return struct


def tail_fold(struct: list[list[int]], method: int, previous: list[int]) -> list[list[int]]:
    """
    Apply a rotation/reflection or diagonal movement to the protein tail.

    Methods:
        1: 90° clockwise rotation
        2: 90° anticlockwise rotation
        3: 180° rotation
        4: x-axis reflection
        5: y-axis reflection
        6: 1-3 quadrant bisector symmetry
        7: 2-4 quadrant bisector symmetry
        8: diagonal move of the first monomer

    Parameters
    ----------
    struct : list[list[int]]
        Structure of the sequence, each element is [x, y] of a monomer.
    method : int
        The transformation method to apply.
    previous : list[int]
        Coordinates [x, y] of the previous monomer (used for diagonal move).

    Returns
    -------
    list[list[int]]
        Transformed structure.
    """
    if method == 8:
        return diagonal_move(struct, previous)

    # Dictionary mapping methods to lambda transformations
    transforms = {
        1: lambda x, y: [y, -x],
        2: lambda x, y: [-y, x],
        3: lambda x, y: [-x, -y],
        4: lambda x, y: [x, -y],
        5: lambda x, y: [-x, y],
        6: lambda x, y: [-y, -x],
        7: lambda x, y: [y, x]
    }

    transform = transforms.get(method)
    if transform is None:
        raise ValueError(f"Invalid method {method}. Must be an integer between 1 and 8.")

    new_tail = [transform(x, y) for x, y in struct]
    return new_tail


def hp_sequence_transform(seq: str) -> str:
    """
    Transform a protein sequence of standard amino acids into the HP model sequence.

    Parameters
    ----------
    seq : str
        Sequence containing the 20 standard amino acids (RNDQEHKSTACGILMFPWYV),
        all must be uppercase letters.

    Returns
    -------
    str
        The sequence converted into H (hydrophobic) and P (polar).

    Raises
    ------
    ValueError
        If any amino acid in the sequence is not recognized.
    """
    polar = set('RNDQEHKST')  # polar amino acids
    hydrophobic = set('ACGILMFPWYV')  # hydrophobic amino acids

    hp_seq = []

    for aa in seq:
        if aa in polar:
            hp_seq.append('P')
        elif aa in hydrophobic:
            hp_seq.append('H')
        else:
            raise ValueError(f"Amino acid '{aa}' not recognized")

    return ''.join(hp_seq)


def progress_bar(progress: int, total: int) -> None:
    """
    Print a progress bar on the terminal when used inside a loop.

    Parameters
    ----------
    progress : int
        Current progress in the loop.
    total : int
        Total number of steps for the evolution.

    Returns
    -------
    None
        Only prints the progress bar to the terminal.
    """
    percentage = progress / float(total) * 100
    filled = int(percentage / 10)  # number of '#' to print (max 10)
    empty = 10 - filled
    bar = f"[{'#' * filled}{' ' * empty}]"
    print(f"\r{bar} {percentage:.2f}%", end="")
    if progress >= total:
        print()  # newline after completion


class Configuration:
    """
    Class to parse and store parameters from an input configuration file.

    Parameters
    ----------
    config : ConfigParser
        ConfigParser object of the input file.

    Attributes
    ----------
    seq : str
        Protein sequence.
    folds : int
        Number of folding steps.
    use_struct : bool
        Whether to use the structure provided in the config file.
    annealing : bool
        Whether to use annealing.
    temperature : float
        Starting temperature.
    struct : list, optional
        Initial structure, if use_struct is True.
    gif : bool
        Whether to create a GIF of the folding.
    seed : int
        Random seed.
    """

    def __init__(self, config) -> None:
        self.seq = config['SEQUENCE']['sequence']
        self.folds = config['PROCESS'].getint('folding_steps')
        self.use_struct = config['optional'].getboolean('use_structure')
        self.annealing = config['optional'].getboolean('annealing')
        self.temperature = config['optional'].getfloat('T')
        self.struct = None

        if self.use_struct:
            struct_str = config['optional']['structure']
            self.struct = json.loads(struct_str)

        self.gif = config['optional'].getboolean('create_gif')

        seed = config['random_seed']['seed']
        if seed == 'None':
            seed = random.randint(0, 10000)
        self.seed = int(seed)
