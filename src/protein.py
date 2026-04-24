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
        self._sequence: str = ""
        self._fold: list[list[int]] = []
        self._pos_map: dict[tuple[int, int], int] = {}
        
        self._init_sequence(config)
        self._init_structure(config)

    @property
    def fold(self) -> list[list[int]]:
        return self._fold

    @fold.setter
    def fold(self, new_fold: list[list[int]]) -> None:
        self._fold = new_fold
        self._pos_map = {tuple(pos): i for i, pos in enumerate(new_fold)}

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
            new_fold = config.fold
        else:
            new_fold = geometry.generate_linear_fold(self.sequence)

        if not validation.is_valid_fold(new_fold):
            raise AssertionError(
                "Invalid structure: not a self-avoiding walk (SAW) "
                "or distances between consecutive points differ from 1"
            )
        
        self.fold = new_fold

    def get_energy(self, epsilon: float = 1.0) -> float:
        """
        Calculate total energy via the HP model.

        Only non-consecutive H-H (hydrophobic) contacts contribute,
        each contact counted once (divided by 2 to avoid double-counting).

        Parameters
        ----------
        epsilon : float
            Energy scale factor (default 1.0).

        Returns
        -------
        float
            Total energy (≤ 0).
        """
        count_hh = 0
        for i in range(self.sequence_length):
            if self.sequence[i] == 'H':
                neighbors = self.get_neighbors(i, self._pos_map)
                count_hh += neighbors.count('H')
        energy = -epsilon * (count_hh * 0.5)
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
            neighbors = self.get_neighbors(i, self._pos_map)
            total += len(neighbors)
        return total // 2

    def get_neighbors(self, i: int, pos_map: dict = None) -> str:
        """
        Return the sequence characters (H/P) of lattice neighbors of monomer i,
        excluding backbone (i-1, i+1) neighbors.

        Parameters
        ----------
        i : int
            Index of the monomer.
        pos_map : dict, optional
            A dictionary mapping (x, y) tuples to monomer indices.
            If None, it is constructed on the fly (less efficient).

        Returns
        -------
        str
            Concatenated sequence characters of non-bonded neighbors.
        """
        if pos_map is None:
            pos_map = {tuple(pos): i for i, pos in enumerate(self.fold)}

        x, y = self.fold[i]
        candidates = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
        
        neighbor_chars = []
        for cand in candidates:
            j = pos_map.get(cand)
            if j is not None:
                # Exclude backbone neighbors
                if j != i - 1 and j != i + 1:
                    neighbor_chars.append(self.sequence[j])
        
        return ''.join(neighbor_chars)
