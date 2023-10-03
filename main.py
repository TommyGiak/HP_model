from protein_class import Protein
import random
import time
import configparser
import argparse
import json
import matplotlib.pyplot as plt


random.seed(0)

# Parser to get from terminal the configuration file
parser = argparse.ArgumentParser()
# the filename is optional, the default is the config.txt
parser.add_argument('configuration_file', help='file from which takes the configuration', default = 'config.txt', nargs='?')

args = parser.parse_args() # get the args written in the terminal
filename = args.configuration_file # assign filename

# setting the parameters from configuration file
config = configparser.ConfigParser()
config.read(filename)

seq = config['SEQUENCE']['sequence'] # selected sequence
folds = config['PROCESS'].getint('folding_steps') # number of folds
use_struct = config['optional'].getboolean('use_structure') # if use the structure present in config file or use linear structure

if use_struct:
    struct = config['optional']['structure']
    struct = json.loads(struct)
    prot = Protein(seq,struct=struct)
else:
    prot = Protein(seq)

prot.view(tit='Initial configuration') # plot the structure of the protein

start = time.time() 

prot.evolution(steps=folds) # evolve the protein with folds foldings
prot.view(save=False, tit='Final configuration') # plot the final structure of the protein
prot.view_min_en()
prot.view_max_comp()
prot.plot_energy(avg=10)
prot.plot_compactness(avg=10)

plt.show()

print(f'It took {time.time()-start:.3f} seconds')
