# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""
from protein_class import Protein
import random
import time
import configparser
import argparse
import utils
import matplotlib.pyplot as plt
import plots


random.seed(123) # random seed for replicability

# Parser to get from terminal the configuration file
parser = argparse.ArgumentParser()
# the filename is optional, the default is the config.txt
parser.add_argument('configuration_file', help='file from which takes the configuration', default = 'config.txt', nargs='?')

args = parser.parse_args() # get the args written in the terminal
filename = args.configuration_file # assign filename

# setting the parameters from configuration file
configuration = configparser.ConfigParser()
configuration.read(filename)

config = utils.Configuration(configuration) # class to get the save the configuration from the file

prot = Protein(config) # Protein class 

plots.view(protein=prot, tit='Initial configuration') # plot the structure of the protein

start = time.time() 

prot.evolution() # evolve the protein with folds foldings
# various plots:
plots.view(protein=prot, save=False, tit='Final configuration')
plots.view_min_en(protein=prot)
plots.view_max_comp(protein=prot)
plots.plot_energy(protein=prot, avg=10)
plots.plot_compactness(protein=prot, avg=10)

print(f'It took {time.time()-start:.3f} seconds')

plt.show() # to let the let plots on screen at the end

if config.gif:
    plots.create_gif(protein=prot)