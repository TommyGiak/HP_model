# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""
from protein_class import Protein
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import numpy as np
import utils

def view(protein : Protein, save = True, tit = None):
    '''
    Function to plot the protein structure with matplotlib.
    As first argument the protein class instance of the desired protein is needed.
    Title can be optionally inserted.
    If save == True the plot will be also saved as pdf.
    '''
    x = [] # x coordinates of the monomers (ordered)
    y = [] # y coordinates of the monomers (ordered)
    
    fig, ax = plt.subplots()
    for i in range(protein.n):
        x.append(protein.struct[i][0])
        y.append(protein.struct[i][1])    
    ax.plot(x,y, alpha = 0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$'+protein.seq[i]+'$', s=20, color = 'red')
    ax.set_xlim(min(x)-6,max(x)+6)
    ax.set_ylim(min(y)-6,max(y)+6)
    ax.grid(alpha=0.2)
    if tit is not None:
        ax.set_title(tit)
    en = protein.energy()
    comp = protein.compactness()
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp/(max(protein.comp_evo)+10e-15):.2f}' # the +10e-15 is used for numerical stability (avoid division by 0)
    ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    plt.show(block=False)
    if save:
        plt.savefig("data/prot_view.png", format="png", bbox_inches="tight", dpi = 200)


def view_min_en(protein, save = True):
    '''
    Function to plot the protein structure founded whit less energy with matplotlib.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf
    '''
    x = [] # x coordinates of the monomers (ordered)
    y = [] # y coordinates of the monomers (ordered)
    
    fig, ax = plt.subplots()
    for i in range(protein.n):
        x.append(protein.min_en_struct[i][0])
        y.append(protein.min_en_struct[i][1])    
    ax.plot(x,y, alpha = 0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$'+protein.seq[i]+'$', s=20, color = 'red')
    ax.set_xlim(min(x)-6,max(x)+6)
    ax.set_ylim(min(y)-6,max(y)+6)
    ax.grid(alpha=0.2)
    en = min(protein.en_evo)
    comp = protein.comp_evo[protein.en_evo.index(en)]
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp/(max(protein.comp_evo)+10e-15):.2f}' # + 10e-15 for numerical stability (avoid division by 0)
    ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    ax.set_title('Min energy structure')
    plt.show(block=False)
    if save:
        plt.savefig("data/min_energy_view.png", format="png", bbox_inches="tight", dpi = 200)

            
def view_max_comp(protein, save = True):
    '''
    Function to plot the protein structure founded whit the max compactness.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf

    '''
    x = [] # x coordinates of the monomers (ordered)
    y = [] # y coordinates of the monomers (ordered)
    
    fig, ax = plt.subplots()
    for i in range(protein.n):
        x.append(protein.max_comp_struct[i][0])
        y.append(protein.max_comp_struct[i][1])    
    ax.plot(x,y, alpha = 0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$'+protein.seq[i]+'$', s=20, color = 'red')
    ax.set_xlim(min(x)-6,max(x)+6)
    ax.set_ylim(min(y)-6,max(y)+6)
    ax.grid(alpha=0.2)
    ax.set_title('Max compactness structure')
    comp = max(protein.comp_evo)
    en = protein.en_evo[protein.comp_evo.index(comp)]
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp/(max(protein.comp_evo)+10e-15):.2f}' # the +10e-15 is used for numerical stability (avoid division by 0)
    ax.text(0.01,0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01,0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    plt.show(block=False)
    if save:
        plt.savefig("data/max_compactness_view.png", format="png", bbox_inches="tight", dpi = 200)

    

def plot_energy(protein, avg : int = 1, save = True) -> None:
    '''
    plot the energy evolution of the system.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf

    Parameters
    ----------
    avg : int, optional
        The energy will be averaged every avg steps. The default is 1.

    Returns
    -------
    Plot
    '''
    en_evo = np.array(protein.en_evo[1:].copy())
    T = np.array(protein.T[1:])
    x = np.arange(0, len(en_evo), avg)
    try:
        en_evo = en_evo.reshape(-1,avg).mean(axis=1)
        T = T.reshape(-1,avg).mean(axis=1)
    except:
        print(f'Mean procedure skipped since the number of time steps is not a multiple of the avarage required: {avg}')
    fig, ax = plt.subplots()
    ax.set_title(f'Energy evolution of the system averaged by {avg} time steps')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Energy', color = 'b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, en_evo, color ='b')
    ax_tw = ax.twinx()
    ax_tw.set_ylabel('T', color = 'r')
    ax_tw.tick_params(axis='y', labelcolor='r')
    ax_tw.plot(x, T, color = 'r')
    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig("data/energy_evolution.png", format="png", bbox_inches="tight", dpi = 200)
    
    
def plot_compactness(protein, avg : int = 1, save = True) -> None:
    '''
    plot the compactness evolution of the system.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf

    Parameters
    ----------
    avg : int, optional
        The compactness will be averaged every avg steps. The default is 1.

    Returns
    -------
    Plot
    '''
    comp = np.array(protein.comp_evo[1:].copy())
    comp = comp/max(comp)
    T = np.array(protein.T[1:])
    x = np.arange(0, len(comp), avg)
    try:
        comp = comp.reshape(-1,avg).mean(axis=1)
        T = T.reshape(-1,avg).mean(axis=1)
    except:
        print(f'Mean procedure skipped since the number of time steps is not a multiple of the avarage required: {avg}')
    fig, ax = plt.subplots()
    ax.set_title(f'Compactness evolution of the system averaged by {avg} time steps')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Compactness', color = 'b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, comp, color ='b')
    ax_tw = ax.twinx()
    ax_tw.set_ylabel('T', color = 'r')
    ax_tw.tick_params(axis='y', labelcolor='r')
    ax_tw.plot(x, T, color = 'r')
    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig("data/compactness_evolution.png", format="png", bbox_inches="tight", dpi = 200)


def create_gif(protein):
    '''
    Function to create a gif of the evolution process.\n
    As first argument the protein class instance of the desired protein is needed.\n
    The gif will have more o less 100 frames with fps = 5.
    The gif will not be visible with the other plots but it will be saved in the /data folder.\n
    To control the creation or not of the gif you can set TRUE or FALSE for the variable 'create_gif' in the configuration file.
    '''
    print('---------------')
    print('Creating gif...')

    fig, ax = plt.subplots()

    writer =  PillowWriter(fps=5)

    with writer.saving(fig, 'data/evo.gif', 200):

        for i,structure in enumerate(protein.gif_struct):
            utils.progress_bar(i+1, len(protein.gif_struct))

            x = []
            y = []

            for i in range(protein.n):
                x.append(structure[i][0])
                y.append(structure[i][1])
            ax.plot(x,y, alpha = 0.5)
            for i, coord in enumerate(structure):
                ax.scatter(x[i], y[i], marker='$'+protein.seq[i]+'$', s=20, color = 'red')
            ax.set_xlim(min(x)-6,max(x)+6)
            ax.set_ylim(min(y)-6,max(y)+6)
            ax.grid(alpha=0.2)

            writer.grab_frame()
            ax.clear()
    
    print('Gif created')
    print('-----------')