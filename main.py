from protein_class import Protein
import random
import time
import configparser
import argparse
import utils
import matplotlib.pyplot as plt


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

prot.view(tit='Initial configuration') # plot the structure of the protein

start = time.time() 

prot.evolution() # evolve the protein with folds foldings
# various plots:
prot.view(save=False, tit='Final configuration')
prot.view_min_en()
prot.view_max_comp()
prot.plot_energy(avg=10)
prot.plot_compactness(avg=10)

print(f'It took {time.time()-start:.3f} seconds')

plt.show() # to let the let plots on screen at the end

if config.gif:
    prot.create_gif()