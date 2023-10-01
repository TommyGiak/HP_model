from protein_class import Protein
import random
import time


random.seed(0)
# structure for tests 
a = [[0, 0],[0, 1],[1, 1],[1, 2],[1, 3],[2, 3],[2, 2],
     [2, 1],[2, 0],[2, -1],[1, -1],[0, -1],[-1, -1],
     [-1, -2],[-1, -3],[-1, -4],[-1, -5]]
seq = 'MAGIIKKQILKHFPKSCDNFNLLHPIFQRHAHEQDTKMHEIYKGNITPQLNKNTLKTSAATDVWAVYFSQFWIDYEGMKSGKGRPISFVDSFPLSIWIC'
prot = Protein(seq) 
prot.view() # plot the structure of the protein

start = time.time()

folds = 2000
prot.evolution(steps=folds) # evolve the protein with folds foldings
prot.view() # plot the final structure of the protein
prot.view_min_en()
prot.view_max_comp()
prot.plot_energy(avg=10)
prot.plot_compactness(avg=10)

print(f'It took {time.time()-start:.3f} seconds')
