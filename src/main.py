"""
Entry point for the HP model protein folding simulation.
@author: Tommaso Giacometti, Alessandro Quirile
"""
import random

import plots
from config import Configuration
from protein import Protein
from simulation import Simulation

if __name__ == "__main__":
    config = Configuration("config.yaml")
    random.seed(config.seed)
    config.print_simulation_setup()

    protein = Protein(config)
    plots.plot_fold(protein, tit="Initial configuration", filename="initial.png")

    simulation = Simulation(protein, config)
    simulation.run()

    plots.generate_plots(protein, simulation.tracker, config)
