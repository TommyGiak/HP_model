"""
Entry point for the HP model protein folding simulation.
@author: Tommaso Giacometti, Alessandro Quirile
"""
import random

import plots
import utils
from protein import Protein
from simulation import Simulation

if __name__ == "__main__":
    config = utils.Configuration("config.yaml")
    random.seed(config.seed)

    protein = Protein(config)
    plots.view(protein, tit="Initial configuration", filename="initial.png")

    utils.print_simulation_setup(config)

    simulation = Simulation(protein, config)
    simulation.run()

    plots.generate_plots(protein, simulation.tracker, config)
