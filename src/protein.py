"""
Protein domain model: sequence, structure, and energy calculations.
@author: Tommaso Giacometti, Alessandro Quirile
"""

import config
import geometry
import utils
import validation


class Protein:
    """
    A protein in the HP lattice model.

    Encapsulates the biological/physical properties of the protein:
    sequence, 2D fold, and energy/compactness calculations.

    Parameters
    ----------
    config : config.Configuration
        Configuration object with sequence and optional initial structure.
    """

    def __init__(self, config: config.Configuration) -> None:
        self._init_sequence(config)
        self._init_structure(config)

    def _init_sequence(self, config: config.Configuration) -> None:
        """Parse and validate the protein sequence, converting to HP if needed."""
        if validation.is_valid_sequence(config.sequence):
            self.sequence = config.sequence
        else:
            self.sequence = utils.convert_to_hp(config.sequence)
        self.sequence_length = len(self.sequence)

    def _init_structure(self, config: config.Configuration) -> None:
        """Initialize the 2D lattice fold, either linear or user-provided."""

        if config.use_struct:
            if len(config.fold) != self.sequence_length:
                raise AssertionError("Sequence and structure lengths do not match")
            self.fold = config.fold
        else:
            self.fold = geometry.generate_linear_fold(self.sequence)

        if not validation.is_valid_fold(self.fold):
            raise AssertionError(
                "Invalid structure: not a self-avoiding walk (SAW) "
                "or distances between consecutive points differ from 1"
            )

    def get_energy(self, e: float = 1.0) -> float:
        """
        Calculate total energy via the HP model.

        Only non-consecutive H-H (hydrophobic) contacts contribute,
        each contact counted once (divided by 2 to avoid double-counting).

        Parameters
        ----------
        e : float
            Energy scale factor (default 1.0).

        Returns
        -------
        float
            Total energy (≤ 0).
        """
        count_hh = 0
        for i in range(self.sequence_length):
            if self.sequence[i] == 'H':
                neighbors = self.get_neighbors(i)
                count_hh += neighbors.count('H')
        energy = -e * (count_hh * 0.5)
        return energy

    def get_compactness(self) -> int:
        """
        Calculate compactness as the number of unique non-covalent topological contacts.

        Returns
        -------
        int
            Total contact count.
        """
        total = 0
        for i in range(self.sequence_length):
            neighbors = self.get_neighbors(i)
            total += len(neighbors)
        return total // 2

    def get_neighbors(self, i: int) -> str:
        """
        Return the sequence characters (H/P) of lattice neighbors of monomer i,
        excluding backbone (i-1, i+1) neighbors.

        Parameters
        ----------
        i : int
            Index of the monomer.

        Returns
        -------
        str
            Concatenated sequence characters of non-bonded neighbors.
        """
        x, y = self.fold[i]
        candidates = [[x - 1, y], [x + 1, y], [x, y - 1], [x, y + 1]]
        if i > 0 and self.fold[i - 1] in candidates:
            candidates.remove(self.fold[i - 1])
        if i < self.sequence_length - 1 and self.fold[i + 1] in candidates:
            candidates.remove(self.fold[i + 1])
        neighbor_indices = [j for j, pos in enumerate(self.fold) if pos in candidates]
        return ''.join(self.sequence[j] for j in neighbor_indices)
