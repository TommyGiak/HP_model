from numpy import isclose

from geometry import get_distance


def is_valid_fold(fold: list[list[int]]) -> bool:
    """
    Check if the given protein fold is valid:
    - Self-avoiding walk (SAW, no overlapping monomers)
    - Distance between consecutive monomers is 1

    Parameters
    ----------
    fold : list[list[int]]
        Fold of the protein containing x and y coordinates for each monomer.

    Returns
    -------
    bool
        True if the fold is valid, False otherwise.
    """
    seen = set()  # positions already visited

    for i, pos in enumerate(fold):
        pos_tuple = tuple(pos)  # convert to tuple to store in set
        if pos_tuple in seen:
            return False
        seen.add(pos_tuple)

        # check distance to next monomer using squared distance for efficiency
        if i < len(fold) - 1:
            dx = fold[i][0] - fold[i+1][0]
            dy = fold[i][1] - fold[i+1][1]
            if (dx**2 + dy**2) != 1:
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
