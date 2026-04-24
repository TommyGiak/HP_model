"""
SimulationTracker: records and maintains the history of a protein folding simulation.
@author: Tommaso Giacometti, Alessandro Quirile
"""
import copy
from typing import List


class SimulationTracker:
    """
    Observes and records the state of a protein folding simulation over time.

    Single Responsibility: observation and record-keeping only —
    no simulation logic, no fold generation, no physical calculations.

    Stores:
    - Energy, compactness, and temperature histories
    - Best (minimum energy) and most compact fold snapshots
    - Optional GIF frame snapshots

    Best-value tracking is O(1) per step (cached scalars, not list scans).

    Parameters
    ----------
    initial_fold : List[List[int]]
        Starting lattice fold.
    initial_energy : float
        Energy of the starting fold.
    initial_compactness : int
        Compactness of the starting fold.
    initial_temperature : float
        Starting simulation temperature.
    do_gif : bool
        Whether to collect GIF frames during the simulation.
    """

    def __init__(
            self,
            initial_fold: List[List[int]],
            initial_energy: float,
            initial_compactness: int,
            initial_temperature: float,
            do_gif: bool,
    ) -> None:
        # Evolution histories (index 0 = initial state, index i+1 = after step i)
        self.energy_evolution: List[float] = [initial_energy]
        self.compactness_evolution: List[int] = [initial_compactness]
        self.temperature_evolution: List[float] = [initial_temperature]

        # Best snapshots
        self.min_energy_fold: List[List[int]] = copy.deepcopy(initial_fold)
        self.max_compactness_fold: List[List[int]] = copy.deepcopy(initial_fold)

        # Cached best scalars — avoids O(n) min/max over the full history each step
        self._min_energy: float = initial_energy
        self._max_compactness: int = initial_compactness

        # GIF frame collection
        self.do_gif: bool = do_gif
        self.gif_struct: List[List[List[int]]] = []
        self.gif_steps: List[int] = []

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def record(
            self,
            step: int,
            fold: List[List[int]],
            energy: float,
            compactness: int,
            temperature: float,
            accepted: bool,
            n_steps: int,
    ) -> None:
        """
        Record all metrics for one completed simulation step.

        Parameters
        ----------
        step : int
            Current step index (0-based).
        fold : List[List[int]]
            Current lattice fold after the step.
        energy : float
            Energy of the current fold.
        compactness : int
            Compactness of the current fold.
        temperature : float
            Temperature at this step.
        accepted : bool
            Whether the Metropolis move was accepted.
        n_steps : int
            Total number of simulation steps (used for GIF sampling).
        """
        self.energy_evolution.append(energy)
        self.compactness_evolution.append(compactness)
        self.temperature_evolution.append(temperature)

        self._update_best_energy(fold, energy)
        self._update_best_compactness(fold, compactness)

        if self.do_gif:
            self._collect_gif_frame(step, fold, n_steps)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _update_best_energy(self, fold: List[List[int]], energy: float) -> None:
        """Update the minimum-energy snapshot if the current energy is lower."""
        if energy < self._min_energy:
            self._min_energy = energy
            self.min_energy_fold = [list(c) for c in fold]

    def _update_best_compactness(self, fold: List[List[int]], compactness: int) -> None:
        """Update the maximum-compactness snapshot if the current value is higher."""
        if compactness > self._max_compactness:
            self._max_compactness = compactness
            self.max_compactness_fold = [list(c) for c in fold]

    def _collect_gif_frame(
            self, step: int, fold: List[List[int]], n_steps: int
    ) -> None:
        """
        Collect a GIF frame: every step for the first 10, then at ~100 evenly spaced steps.
        """
        if step < 10 or step % max(1, n_steps // 100) == 0:
            self.gif_struct.append([list(c) for c in fold])
            self.gif_steps.append(step)
