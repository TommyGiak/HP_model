"""
@author: Tommaso Giacometti
"""
import argparse
import configparser
import random
import time

import matplotlib.pyplot as plt

import plots
import utils
from protein_class import Protein

if __name__ == '__main__':
    # Parser to get from terminal the configuration file
    parser = argparse.ArgumentParser()

    # the filename is optional, the default is the config.txt
    parser.add_argument('configuration_file',
                        help='file from which takes the configuration',
                        default='config.txt',
                        nargs='?')

    args = parser.parse_args()  # get the args written in the terminal
    filename = args.configuration_file  # assign filename

    # setting the parameters from configuration file
    configuration = configparser.ConfigParser()
    configuration.read(filename)

    config = utils.Configuration(configuration)  # class to get the save the configuration from the file

    # Random seed setting
    random.seed(config.seed)
    # print(f'The random seed used is {config.seed}')

    prot = Protein(config)  #  Protein class

    plots.view(protein=prot, tit='Initial configuration')  # plot the structure of the protein

    start = time.time()

    print('--------------------')
    print('Evolution started...')
    prot.evolution()
    print('Evolution ended')
    print(f'It took {time.time() - start:.3f} seconds')
    print('---------------')

    plots.view(protein=prot, save=False, tit='Final configuration')
    plots.view_min_en(protein=prot)
    plots.view_max_comp(protein=prot)
    plots.plot_energy(protein=prot, avg=10)
    plots.plot_compactness(protein=prot, avg=10)

    plt.show()

    if config.gif:
        plots.create_gif(protein=prot)
