"""
FoldSampler: generates valid random folds for a protein on a 2D lattice.
@author: Tommaso Giacometti, Alessandro Quirile
"""
import math
import random
from typing import TYPE_CHECKING, List, Optional

import geometry
import transforms
import validation

if TYPE_CHECKING:
    from protein import Protein


def _attempt_fold(protein: "Protein") -> Optional[List[List[int]]]:
    """
    Try a single random fold perturbation.

    Returns the new fold if valid (SAW), or None if it self-intersects.
    """
    index = random.randint(1, protein.sequence_length - 2)
    pivot = protein.fold[index]
    x, y = pivot

    # Tail includes the pivot and all subsequent monomers
    tail = [list(coord) for coord in protein.fold[index:]]

    # Diagonal move is allowed when previous and next monomers are not aligned
    # Previous: fold[index-1], Next: fold[index+1]
    dx = protein.fold[index - 1][0] - protein.fold[index + 1][0]
    dy = protein.fold[index - 1][1] - protein.fold[index + 1][1]
    diag_move = (dx**2 + dy**2) == 2

    # Shift tail to origin (pivot becomes [0,0])
    tail = [[cx - x, cy - y] for cx, cy in tail]
    # 'previous' monomer relative to the pivot [0,0]
    previous = [
        protein.fold[index - 1][0] - x,
        protein.fold[index - 1][1] - y,
    ]

    method = random.randint(1, 8) if diag_move else random.randint(1, 7)
    tail = transforms.tail_fold(tail, method, previous)

    # Shift tail back to correct lattice position
    tail = [[cx + x, cy + y] for cx, cy in tail]
    new_fold = protein.fold[:index] + tail

    return new_fold if validation.is_valid_fold(new_fold) else None


class FoldSampler:
    """
    Generates valid random folds (self-avoiding walks, SAW) for a Protein.

    Fold-generation geometry: rotations, reflections, diagonal moves, and SAW validation.
    No knowledge of simulation or tracking.

    New move strategies can be added by extending this class
    without modifying existing simulation code.
    """

    MAX_TRIES: int = 10 ** 9

    def sample(self, protein: "Protein") -> List[List[int]]:
        """
        Generate a new valid random fold based on the protein's current fold.

        Parameters
        ----------
        protein : Protein
            The protein whose fold will be perturbed.

        Returns
        -------
        List[List[int]]
            A new valid fold (self-avoiding walk).

        Raises
        ------
        ValueError
            If the sequence is too short to fold (< 3 monomers).
        RuntimeError
            If no valid fold is found within MAX_TRIES attempts.
        """
        if protein.sequence_length < 3:
            raise ValueError("Sequence too short for folding (minimum length: 3)")

        for _ in range(self.MAX_TRIES):
            candidate = _attempt_fold(protein)
            if candidate is not None:
                return candidate

        raise RuntimeError(
            f"Could not find a valid fold after {self.MAX_TRIES} attempts"
        )
