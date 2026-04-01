"""
Simulation: orchestrates the protein folding Metropolis evolution.
@author: Tommaso Giacometti, Alessandro Quirile
"""
from typing import Tuple

from tqdm.auto import tqdm

import utils
from fold_sampler import FoldSampler
from metropolis import metropolis_rule
from protein import Protein
from tracker import SimulationTracker


class Simulation:
    """
    Orchestrates the Metropolis-based protein folding simulation with optional annealing.

    Single Responsibility: run the simulation loop.
    Delegates fold generation to FoldSampler (Dependency Inversion — depends on the
    abstraction of a sampler, not on fold geometry) and record-keeping to SimulationTracker.

    Parameters
    ----------
    protein : Protein
        The protein to fold.
    config : utils.Configuration
        Simulation configuration (steps, annealing, temperature, GIF).
    sampler : FoldSampler, optional
        Fold generation strategy. Defaults to FoldSampler().
        Can be replaced with any object implementing `.sample(protein)` for extensibility.
    """

    def __init__(self, protein: Protein, config: utils.Configuration, sampler: FoldSampler = None) -> None:
        self.protein = protein
        self.n_steps = config.n_steps
        self.do_annealing = config.do_annealing
        self.starting_temperature = config.temperature

        self._sampler = sampler if sampler is not None else FoldSampler()

        self.tracker = SimulationTracker(
            initial_fold=protein.fold,
            initial_energy=protein.get_energy(),
            initial_compactness=protein.get_compactness(),
            initial_temperature=config.temperature,
            do_gif=config.do_gif,
        )

    def run(self) -> None:
        """Run the full simulation for n_steps Metropolis steps."""
        temperature = self.starting_temperature

        for step in tqdm(range(self.n_steps), desc="Evolution"):
            temperature = self._update_temperature(temperature, step)
            accepted, energy, compactness = self._metropolis_step(temperature)

            self.tracker.record(
                step=step,
                fold=self.protein.fold,
                energy=energy,
                compactness=compactness,
                temperature=temperature,
                accepted=accepted,
                n_steps=self.n_steps,
            )

    def _update_temperature(self, temperature: float, step: int) -> float:
        """Apply linear annealing schedule if enabled."""
        if self.do_annealing and temperature > 0.002:
            return self.starting_temperature * (1 - step / self.n_steps)
        return temperature

    def _metropolis_step(self, temperature: float) -> Tuple[bool, float, int]:
        """
        Attempt a fold move and apply the Metropolis acceptance criterion.

        Returns
        -------
        Tuple[bool, float, int]
            (accepted, energy_after_step, compactness_after_step)
        """
        current_energy = self.protein.get_energy()
        current_fold = [list(c) for c in self.protein.fold]

        self.protein.fold = self._sampler.sample(self.protein)
        new_energy = self.protein.get_energy()

        accepted = metropolis_rule(current_energy, new_energy, temperature)

        if not accepted:
            self.protein.fold = current_fold
            return False, current_energy, self.protein.get_compactness()

        return True, new_energy, self.protein.get_compactness()
