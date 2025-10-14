"""
@author: Tommaso Giacometti
"""
import random
import time

import plots
import utils
from protein_class import Protein

if __name__ == '__main__':
    config = utils.Configuration("config.yaml")  # class to get the save the configuration from the file

    # Random seed setting
    random.seed(config.seed)

    protein = Protein(config)

    # plots.view(protein, tit='Initial configuration')

    start = time.time()

    print('--------------------')
    print('Evolution started...')
    protein.evolution()
    print('Evolution ended')
    print(f'It took {time.time() - start:.3f} seconds')
    print('---------------')

    plots.view(protein, save=True, tit='Final configuration')
    plots.view_min_en(protein)
    plots.view_max_comp(protein)
    plots.plot_energy(protein, save=True)
    plots.plot_compactness(protein, save=True)

    # plt.show()

    if config.gif:
        plots.create_gif(protein)
