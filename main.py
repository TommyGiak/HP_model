from protein_class import Protein
import random
import time
import configparser


random.seed(0)

config = configparser.ConfigParser()
config.read('config.txt')

seq = config['SEQUENCE']['sequence']
folds = config['PROCESS'].getint('folding_steps')
use_struct = config['optional'].getboolean('use_structure')

if use_struct:
    struct = config['optional']['structure']
    struct = list(struct)
    prot = Protein(seq,struct=struct)
else:
    prot = Protein(seq)

prot.view() # plot the structure of the protein

start = time.time() 

prot.evolution(steps=folds) # evolve the protein with folds foldings
prot.view() # plot the final structure of the protein
prot.view_min_en()
prot.view_max_comp()
prot.plot_energy(avg=10)
prot.plot_compactness(avg=10)

print(f'It took {time.time()-start:.3f} seconds')
