from numpy import sqrt


def generate_linear_fold(seq: str) -> list[list[int]]:
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
    return struct


def get_distance(coord1: list[float], coord2: list[float]) -> float:
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
