"""
Protein folding simulation class
Refactored for readability and modularity
@author: Tommaso Giacometti
"""
import copy
import math
import random
from typing import List

from tqdm.auto import tqdm

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

    def _initialize_sequence(self, config: utils.Configuration) -> None:
        """Initialize protein sequence, converting to HP if needed."""
        if utils.is_valid_sequence(config.sequence):
            self.sequence = config.sequence
        else:
            self.sequence = utils.hp_sequence_transform(config.sequence)
        self.sequence_length = len(self.sequence)

    def _initialize_structure(self, config: utils.Configuration) -> None:
        """Initialize protein structure."""
        if not config.use_struct:
            self.fold = utils.linear_fold(self.sequence)
        else:
            if len(config.fold) != self.sequence_length:
                raise AssertionError("Sequence and structure lengths do not match")
            self.fold = config.fold

        if not utils.is_valid_fold(self.fold):
            raise AssertionError(
                "Invalid structure: not a self-avoiding walk (SAW) "
                "or distances between consecutive points differ from 1"
            )

    def _initialize_parameters(self, config: utils.Configuration) -> None:
        """Load simulation parameters from configuration."""
        self.do_annealing = config.do_annealing
        self.starting_temperature = config.temperature
        self.n_steps = config.n_steps
        self.do_gif = config.do_gif

    def _initialize_tracking(self) -> None:
        """Initialize tracking lists for energy, temperature, compactness, and GIF frames."""
        self.gif_struct: List[List[List[int]]] = []
        self.gif_steps: List[int] = []
        self.min_energy_fold = copy.deepcopy(self.fold)
        self.max_compactness_fold = copy.deepcopy(self.fold)
        self.energy_evolution: List[float] = [self.get_energy()]
        self.temperature_evolution: List[float] = [self.starting_temperature]
        self.compactness_evolution: List[int] = [self.get_compactness()]
        # self.n_foldings: int = 0

    def evolution(self) -> None:
        """Let the system evolve for `self.steps` steps using the Metropolis algorithm.
        For GIF: show all frames for the first 10 steps, then sample to avoid too many frames.
        """
        temperature = self.starting_temperature

        for step in tqdm(range(self.n_steps), desc="Evolution"):

            # Annealing temperature decreases
            if self.do_annealing and temperature > 0.002:
                temperature = self.starting_temperature * (1 - step / self.n_steps)

            # Metropolis rule (folding happens here)
            self._step_metropolis(temperature, step)

            # Save temperature
            self.temperature_evolution.append(temperature)

    def _step_metropolis(self, temperature: float, step: int) -> None:
        """Perform a single step of the Metropolis evolution and record accepted folds."""
        current_energy = self.get_energy()
        current_fold = [list(coord) for coord in self.fold]

        self.fold = self.random_fold()
        new_energy = self.get_energy()

        if new_energy < current_energy:
            accepted = True
            self.energy_evolution.append(new_energy)
        else:
            if self._is_accepted(current_energy, new_energy, temperature):
                accepted = True
                self.energy_evolution.append(new_energy)
            else:
                accepted = False
                self.fold = current_fold
                self.energy_evolution.append(current_energy)

        # Update min energy fold
        if new_energy < min(self.energy_evolution[:-1]):
            self.min_energy_fold = [list(coord) for coord in self.fold]

        # Update max compactness fold
        compactness = self.get_compactness()
        self.compactness_evolution.append(compactness)
        if compactness > max(self.compactness_evolution[:-1]):
            self.max_compactness_fold = [list(coord) for coord in self.fold]

        # Save to GIF only if fold was accepted
        if self.do_gif and accepted:
            if step < 10 or step % max(1, self.n_steps // 100) == 0:
                self.gif_struct.append([list(coord) for coord in self.fold])
                self.gif_steps.append(step)

    def _is_accepted(self, current: float, new: float, temperature: float) -> bool:
        """Return True if the new structure is accepted by the Metropolis criterion."""
        delta_e = new - current
        prob = math.exp(-delta_e / temperature)
        return prob >= random.uniform(0, 1)

    def get_energy(self, e: float = 1.0) -> float:
        """Compute the total energy of the protein structure."""
        count_hh = sum(self.get_neighbors(i).count('H') for i in range(self.sequence_length))
        return -e * count_hh / 2

    def get_compactness(self) -> int:
        """Compute the compactness of the protein structure."""
        return sum(len(self.get_neighbors(i)) for i in range(self.sequence_length))

    def get_neighbors(self, i: int) -> str:
        """Return string of H/P neighbors for the i-th monomer (excluding backbone)."""
        x, y = self.fold[i]
        candidates = [[x - 1, y], [x + 1, y], [x, y - 1], [x, y + 1]]
        if i > 0 and self.fold[i - 1] in candidates:
            candidates.remove(self.fold[i - 1])
        if i < self.sequence_length - 1 and self.fold[i + 1] in candidates:
            candidates.remove(self.fold[i + 1])
        neighbor_indices = [j for j, struct_pos in enumerate(self.fold) if struct_pos in candidates]
        return ''.join(self.sequence[j] for j in neighbor_indices)

    def random_fold(self) -> List[List[int]]:
        """Generate a valid random folding for the protein."""
        max_tries = 10 ** 9
        for _ in range(max_tries):
            if self.sequence_length < 3:
                raise ValueError("Sequence too short for folding")
            index = random.randint(1, self.sequence_length - 2)
            x, y = self.fold[index]
            tail = [list(coord) for coord in self.fold[index:]]

            # Determine if diagonal move is allowed
            dist = utils.get_distance(self.fold[index - 1], self.fold[index + 1])
            diag_move = math.isclose(dist, math.sqrt(2))

            # Shift tail and previous monomer to origin
            tail = [[cx - x, cy - y] for cx, cy in tail]
            previous = [self.fold[index - 1][0] - x, self.fold[index - 1][1] - y]

            method = random.randint(1, 8) if diag_move else random.randint(1, 7)
            tail = utils.tail_fold(tail, method, previous)

            # Shift tail back to correct position
            tail = [[cx + x, cy + y] for cx, cy in tail]
            new_struct = self.fold[:index] + tail

            if utils.is_valid_fold(new_struct):
                # self.n_foldings += 1
                return new_struct

        raise RuntimeError(f"Could not find a valid fold after {max_tries} attempts")
