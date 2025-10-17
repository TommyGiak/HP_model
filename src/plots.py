"""
@author: Tommaso Giacometti
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from tqdm.auto import tqdm

from protein_class import Protein


def view(protein: Protein, save=True, tit=None, filename="protein.png"):
    '''
    Function to plot the protein structure with matplotlib.
    As first argument the protein class instance of the desired protein is needed.
    Title can be optionally inserted.
    If save == True the plot will be also saved as pdf.
    '''
    x = []  # x coordinates of the monomers (ordered)
    y = []  # y coordinates of the monomers (ordered)

    fig, ax = plt.subplots()
    for i in range(protein.sequence_length):
        x.append(protein.struct[i][0])
        y.append(protein.struct[i][1])
    ax.plot(x, y, alpha=0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')
    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)
    if tit is not None:
        ax.set_title(tit)
    en = protein.energy()
    comp = protein.compactness()
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp / (max(protein.compactness_evolution) + 10e-15):.2f}'  # the +10e-15 is used for numerical stability (avoid division by 0)
    ax.text(0.01, 0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def view_min_energy(protein, save=True, filename="min_energy.png"):
    '''
    Function to plot the protein structure founded whit less energy with matplotlib.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf
    '''
    x = []  # x coordinates of the monomers (ordered)
    y = []  # y coordinates of the monomers (ordered)

    fig, ax = plt.subplots()
    for i in range(protein.sequence_length):
        x.append(protein.min_energy_structure[i][0])
        y.append(protein.min_energy_structure[i][1])
    ax.plot(x, y, alpha=0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')
    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)
    en = min(protein.energy_evolution)
    comp = protein.compactness_evolution[protein.energy_evolution.index(en)]
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp / (max(protein.compactness_evolution) + 10e-15):.2f}'  #  + 10e-15 for numerical stability (avoid division by 0)
    ax.text(0.01, 0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    ax.set_title('Min energy structure')
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def view_max_compactness(protein, save=True, filename="max_compactness.png"):
    '''
    Function to plot the protein structure founded whit the max compactness.
    As first argument the protein class instance of the desired protein is needed.
    The plot can be saved with save = True as pdf

    '''
    x = []  # x coordinates of the monomers (ordered)
    y = []  # y coordinates of the monomers (ordered)

    fig, ax = plt.subplots()
    for i in range(protein.sequence_length):
        x.append(protein.max_compactness_structure[i][0])
        y.append(protein.max_compactness_structure[i][1])
    ax.plot(x, y, alpha=0.5)
    for i, coord in enumerate(protein.struct):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')
    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)
    ax.set_title('Max compactness structure')
    comp = max(protein.compactness_evolution)
    en = protein.energy_evolution[protein.compactness_evolution.index(comp)]
    string = f'Energy: {en}'
    string_comp = f'Compactness: {comp / (max(protein.compactness_evolution) + 10e-15):.2f}'  # the +10e-15 is used for numerical stability (avoid division by 0)
    ax.text(0.01, 0.99, string, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, string_comp, ha='left', va='top', transform=ax.transAxes)
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def plot_energy(protein, avg: int = 10, save=True, filename="energy_evolution.png", annealing=True) -> None:
    '''
    Plot the energy evolution of the system.
    If annealing is False, temperature curve is not shown.
    '''
    en_evo = np.array(protein.energy_evolution[1:].copy())
    x = np.arange(0, len(en_evo), avg)

    fig, ax = plt.subplots()
    try:
        en_evo = en_evo.reshape(-1, avg).mean(axis=1)
    except:
        print(f'Mean procedure skipped for energy: {avg} steps')

    ax.set_title(f'Energy evolution of the system averaged by {avg} steps')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Energy', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, en_evo, color='b')

    if annealing:
        T = np.array(protein.temperature_evolution[1:])
        try:
            T = T.reshape(-1, avg).mean(axis=1)
        except:
            pass
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color='r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T, color='r')

    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def plot_compactness(protein, avg: int = 10, save=True, filename="compactness_evolution.png", annealing=True) -> None:
    '''
    Plot the compactness evolution of the system.
    If annealing is False, temperature curve is not shown.
    '''
    comp = np.array(protein.compactness_evolution[1:].copy())
    comp = comp / max(comp)
    x = np.arange(0, len(comp), avg)

    fig, ax = plt.subplots()
    try:
        comp = comp.reshape(-1, avg).mean(axis=1)
    except:
        print(f'Mean procedure skipped for compactness: {avg} steps')

    ax.set_title(f'Compactness evolution of the system averaged by {avg} steps')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Compactness', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, comp, color='b')

    if annealing:
        T = np.array(protein.temperature_evolution[1:])
        try:
            T = T.reshape(-1, avg).mean(axis=1)
        except:
            pass
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color='r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T, color='r')

    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def create_gif(protein, filename="evolution.gif") -> None:
    """
    Create a GIF of the evolution process using saved frames.
    Displays the real simulation step in the title.
    Progress is shown with tqdm.
    """
    fig, ax = plt.subplots()
    writer = PillowWriter(fps=2)

    with writer.saving(fig, f'output/{filename}', 200):
        for frame, step_real in tqdm(zip(protein.gif_struct, protein.gif_steps), total=len(protein.gif_struct),
                                     desc="Creating GIF"):
            x = [coord[0] for coord in frame]
            y = [coord[1] for coord in frame]

            ax.plot(x, y, alpha=0.5)
            for j, coord in enumerate(frame):
                ax.scatter(x[j], y[j], marker='$' + protein.sequence[j] + '$', s=20, color='red')

            ax.set_xlim(min(x) - 6, max(x) + 6)
            ax.set_ylim(min(y) - 6, max(y) + 6)
            ax.grid(alpha=0.2)
            ax.set_title(f"Evolution step {step_real}/{protein.steps}", fontsize=10)

            fig.canvas.draw()
            writer.grab_frame()
            ax.clear()
