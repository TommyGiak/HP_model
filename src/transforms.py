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
