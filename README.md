# HP model for protein folding

In my project I will implement in Python the HP model for protein folding.
The HP model is a simply model that can help to understand basic folding behaviours using Monte Carlo simulation over the free energy of the bounds using only two category of aminoacids, H (hydrophobic) and P (polar) ([click here for more info](https://pubs.acs.org/doi/10.1021/ma00200a030)).\
I want to create a command line application that can take the sequence of a protein (even with the 20 distinct aminoacids) and can run a HP model simulation of the folding process of the protein at given temperatures using, optionally, annealing algorithms and giving as outputs the energy evolution of the system, the structure of the protein the native structure energy and the compactness.\
This project can be used to have a first impact of the behaviour of a protein and can be used to study the transitions to native states of proteins as function of temperature and bounding energy. Can be also done different test and comparison, for example seeing the different behaviours of similar proteins when two random adjacent aminoacids are switched.

## Run the code

From the _terminal_ create an empty folder and clone this repository using the following command:

```shell
git clone https://github.com/TommyGiak/HP_model.git
```

Once the github repository is cloned run:

```shell
python main.py
```

## Parameters setting

### Insert the protein sequence

The protein sequence can be written in the _config.txt_ file at the _sequence_ variable. The sequence is caps sensitive so the letters must be upper case only. The sequence doesn't need the quotation marks " or ' .

### Change the number of folding steps

The number of folding steps can be setted in the _config.txt_ file at the _folding\_step_ variable.

## Repository strurcture

Up to now the repository contain three pyhton files: _protein_class.py_, _utils.py_, _test.py_.

- _main.py_: run the scripts
- _protein_class.py_: contain a class named `Protein` which contains and save the protein information and implement all the main function for the evolution of the system, including the plot functions.
- _utils.py_: contain different functions to validate the structures of the proteins and to support the evolution of the protein.
- _test.py_: contain the test functions to test all the code.

## Algorithm for the protein folding

The algorithm for the protein folding is implemented in the `Protein` class in the _protein_class.py_ file, plus the function `diagonal_move` and `tail_fold` in the _utils.py_ file.\
Each folding step involve the following steps:

1. choose a random monomer in the protein: a random integer from $0$ to $l-1$ where $l$ is the lenght of the protein sequence. The sampled monomer will be the starting point for the movement of the protein.
1. sample a random integer in $[1,8]$ to randomly select a type of move, the movement are implemented in the function [`tail_fold`](https://github.com/TommyGiak/HP_model/blob/main/utils.py) in the _utils.py_ file. The movement are: 1 = 90° clockwise rotation, 2 = 90° anticlockwise rotation, 3 = 180° rotaion, 4 = x-axis refletion, 5 = y-axis reflection, 6 = 1 and 3 quadrant bisector symmetry, 7 = 2 and 4 quadrant bisector symmetry and 8 = movement on a digaonal of a random monomer. The movement can also be choosen a priori, if it is not specified a random one is selected.
1. the new protein is validated: if the protein sequence overlap (or the distance between neighbours is different from one) the process restart from the step 1.
1. if the protein structure is valid the new structure (the folded protein) is passed.

## Acceptance of the structure

Once a new structure is generated by a folding step, the energy of this new structure is computed to accept it following the Metropolis algorithm.\
If the energy of the new protein structure is less than the previous one, the new structure is always accepted, if instead the energy is grater the new structure is accepted with probability:

```math
p = e^{-\frac{\Delta E}{k_bT}}
```

where $\Delta E > 0$.
