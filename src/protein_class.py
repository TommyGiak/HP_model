"""
Protein folding simulation class
Refactored for readability and modularity
@author: Tommaso Giacometti
"""
import random
from typing import List

import math

import utils


class Protein:
    """
    Protein class containing all the information on the protein and methods to evolve the system.

    Parameters
    ----------
    config : utils.Configuration
        Configuration class with parameters and initial structure.
    """

    def __init__(self, config: utils.Configuration) -> None:
        self._initialize_sequence(config)
        self._initialize_structure(config)
        self._initialize_parameters(config)
        self._initialize_tracking()

    def print_config(self) -> None:
        """Print configuration options selected by the user."""
        print(f"Sequence: {self.sequence}")
        print(f"Folding steps: {self.steps}")
        print(f"Annealing enabled: {self.annealing}")
        print(f"Starting temperature: {self.starting_temperature}")

    def _initialize_sequence(self, config: utils.Configuration) -> None:
        """Initialize protein sequence, converting to HP if needed."""
        print("Sequence:", config.seq)
        if utils.is_valid_sequence(config.seq):
            self.sequence = config.seq
        else:
            self.sequence = utils.hp_sequence_transform(config.seq)
            # print("HP representation:", self.sequence)
        self.sequence_length = len(self.sequence)

    def _initialize_structure(self, config: utils.Configuration) -> None:
        """Initialize protein structure."""
        if not config.use_struct:
            self.struct = utils.linear_struct(self.sequence)
        else:
            if len(config.struct) != self.sequence_length:
                raise AssertionError("Sequence and structure lengths do not match")
            self.struct = config.struct

        if not utils.is_valid_struct(self.struct):
            raise AssertionError(
                "Invalid structure: not a self-avoiding walk (SAW) "
                "or distances between consecutive points differ from 1"
            )

    def _initialize_parameters(self, config: utils.Configuration) -> None:
        """Load simulation parameters from configuration."""
        self.annealing = config.annealing
        self.starting_temperature = config.temperature
        self.steps = config.folds
        self.gif = config.gif

    def _initialize_tracking(self) -> None:
        """Initialize tracking lists for energy, temperature, compactness, and GIF frames."""
        self.gif_struct: List[List[List[int]]] = []
        self.min_energy_structure = self.struct.copy()
        self.max_comp_struct = self.struct.copy()
        self.energy_evolution: List[float] = [self.energy()]
        self.temperature_evolution: List[float] = [self.starting_temperature]
        self.compactness_evolution: List[int] = [self.compactness()]
        self.n_foldings: List[int] = []

    def evolution(self) -> None:
        """Let the system evolve for `self.steps` steps using the Metropolis algorithm."""
        self.print_config()
        temp = self.starting_temperature
        annealing_coef = -temp / self.steps

        for step in range(self.steps):
            utils.progress_bar(step + 1, self.steps)

            # Annealing temperature decrease
            if self.annealing and temp > 0.002:
                temp = annealing_coef * (step - self.steps)

            self._step_metropolis(temp)

            # Save temperature and GIF frame
            self.temperature_evolution.append(temp)
            if self.gif and (step % max(1, self.steps // 100) == 0):
                self.gif_struct.append([coord.copy() for coord in self.struct])

    def _step_metropolis(self, temperature: float) -> None:
        """Perform a single step of the Metropolis evolution."""
        current_energy = self.energy()
        current_struct = [coord.copy() for coord in self.struct]

        self.struct = self.random_fold()
        new_energy = self.energy()

        # Accept new structure probabilistically if energy increases
        if new_energy > current_energy:
            if not self._accept_higher_energy(current_energy, new_energy, temperature):
                self.struct = current_struct

        # Update min energy and max compact structures
        if new_energy < min(self.energy_evolution):
            self.min_energy_structure = [coord.copy() for coord in self.struct]

        self.energy_evolution.append(new_energy)

        comp = self.compactness()
        self.compactness_evolution.append(comp)
        if comp > max(self.compactness_evolution[:-1]):
            self.max_comp_struct = [coord.copy() for coord in self.struct]

    def _accept_higher_energy(self, current: float, new: float, temperature: float) -> bool:
        """Return True if the new structure is accepted by the Metropolis criterion."""
        delta_e = new - current
        prob = math.exp(-delta_e / temperature)
        return prob >= random.uniform(0, 1)

    def energy(self, e: float = 1.0) -> float:
        """Compute the total energy of the protein structure."""
        count_hh = sum(self.get_neighbors(i).count('H') for i in range(self.sequence_length))
        return -e * count_hh / 2

    def compactness(self) -> int:
        """Compute the compactness of the protein structure."""
        return sum(len(self.get_neighbors(i)) for i in range(self.sequence_length))

    def get_neighbors(self, i: int) -> str:
        """Return string of H/P neighbors for the i-th monomer (excluding backbone)."""
        x, y = self.struct[i]
        candidates = [[x - 1, y], [x + 1, y], [x, y - 1], [x, y + 1]]
        if i > 0: candidates.remove(self.struct[i - 1])
        if i < self.sequence_length - 1: candidates.remove(self.struct[i + 1])
        return ''.join(self.sequence[self.struct.index(pos)] for pos in candidates if pos in self.struct)

    def random_fold(self) -> List[List[int]]:
        """Generate a valid random folding for the protein."""
        while True:
            index = random.randint(1, self.sequence_length - 2)
            x, y = self.struct[index]
            tail = [coord.copy() for coord in self.struct[index:]]

            # Determine if diagonal move is allowed
            dist = utils.get_dist(self.struct[index - 1], self.struct[index + 1])
            diag_move = math.isclose(dist, math.sqrt(2))

            # Shift tail and previous monomer to origin
            tail = [[cx - x, cy - y] for cx, cy in tail]
            previous = [self.struct[index - 1][0] - x, self.struct[index - 1][1] - y]

            method = random.randint(1, 8) if diag_move else random.randint(1, 7)
            tail = utils.tail_fold(tail, method, previous)

            # Shift tail back to correct position
            tail = [[cx + x, cy + y] for cx, cy in tail]
            new_struct = self.struct[:index] + tail

            self.n_foldings.append(1)
            if utils.is_valid_struct(new_struct):
                return new_struct
