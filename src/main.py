"""
@author: Tommaso Giacometti
"""

import random

import plots
import utils
from protein_class import Protein


def print_simulation_setup(config):
    """Display the main simulation configuration."""
    print("[HP Model] Simulation setup")
    print(f"Sequence:           {config.sequence}")
    print(f"Sequence length:    {len(config.sequence)}")
    print(f"Structure:          {'Linear' if not config.use_struct else 'Non-linear'}")
    print(f"Folding steps:      {config.n_steps}")
    print(f"Annealing:          {config.do_annealing}")
    print(f"Temperature:        {config.temperature if config.do_annealing else 'Undefined'}")


def generate_plots(protein, config):
    """Generate all plots and optionally create a GIF."""
    plots.view(protein, tit='Final configuration', filename="final.png")
    plots.view_min_energy(protein, filename="min_energy.png")
    plots.view_max_compactness(protein, filename="max_compactness.png")
    plots.plot_energy(protein, filename="energy_evolution.png", annealing=config.do_annealing)
    plots.plot_compactness(protein, filename="compactness_evolution.png", annealing=config.do_annealing)

    if config.do_gif:
        plots.create_gif(protein, filename="evolution.gif")


if __name__ == "__main__":
    config = utils.Configuration("config.yaml")
    random.seed(config.seed)

    protein = Protein(config)
    plots.view(protein, tit="Initial configuration", filename="initial.png")

    print_simulation_setup(config)

    protein.evolution()

    generate_plots(protein, config)
