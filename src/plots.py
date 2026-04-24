"""
Plotting utilities for protein folding simulation results.
@author: Tommaso Giacometti
"""
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from tqdm.auto import tqdm

from protein import Protein
from tracker import SimulationTracker


def generate_plots(
        protein: Protein,
        tracker: SimulationTracker,
        config,
) -> None:
    """Generate all output plots and optionally create the evolution GIF."""
    plot_fold(protein, tracker, tit="Final configuration", filename="final.png")
    plot_min_energy_fold(protein, tracker, filename="min_energy.png")
    plot_max_compactness_fold(protein, tracker, filename="max_compactness.png")
    plot_energy(tracker, filename="energy_evolution.png", annealing=config.do_annealing)
    plot_compactness(tracker, filename="compactness_evolution.png", annealing=config.do_annealing)

    if config.do_gif:
        create_gif(protein, tracker, filename="evolution.gif")


def plot_fold(
        protein: Protein,
        tracker: Optional[SimulationTracker] = None,
        save: bool = True,
        tit: Optional[str] = None,
        filename: str = "protein.png",
) -> None:
    """
    Plot the current protein fold.

    Parameters
    ----------
    protein : Protein
        Protein to visualize.
    tracker : SimulationTracker, optional
        If provided, used to compute the normalized compactness label.
        If None (e.g., before simulation starts), compactness ratio defaults to 1.0.
    save : bool
        Save the figure to output/.
    tit : str, optional
        Plot title.
    filename : str
        Output filename.
    """
    x = [protein.fold[i][0] for i in range(protein.sequence_length)]
    y = [protein.fold[i][1] for i in range(protein.sequence_length)]

    fig, ax = plt.subplots()
    ax.plot(x, y, alpha=0.5)
    for i in range(protein.sequence_length):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')

    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)

    if tit is not None:
        ax.set_title(tit)

    energy = protein.get_energy()
    compactness = protein.get_compactness()

    if tracker is not None:
        max_comp = max(tracker.compactness_evolution) + 1e-14
    else:
        max_comp = compactness + 1e-14  # ratio → 1.0 before simulation

    ax.text(0.01, 0.99, f'Energy: {energy}', ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, f'Compactness: {compactness / max_comp:.2f}', ha='left', va='top', transform=ax.transAxes)

    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def plot_min_energy_fold(
        protein: Protein,
        tracker: SimulationTracker,
        save: bool = True,
        filename: str = "min_energy.png",
) -> None:
    """Plot the minimum-energy fold recorded during the simulation."""
    fold = tracker.min_energy_fold
    x = [fold[i][0] for i in range(protein.sequence_length)]
    y = [fold[i][1] for i in range(protein.sequence_length)]

    fig, ax = plt.subplots()
    ax.plot(x, y, alpha=0.5)
    for i in range(protein.sequence_length):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')

    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)
    ax.set_title('Min energy structure')

    energy = min(tracker.energy_evolution)
    idx = tracker.energy_evolution.index(energy)
    compactness = tracker.compactness_evolution[idx]
    max_comp = max(tracker.compactness_evolution) + 1e-14

    ax.text(0.01, 0.99, f'Energy: {energy}', ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, f'Compactness: {compactness / max_comp:.2f}', ha='left', va='top', transform=ax.transAxes)

    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def plot_max_compactness_fold(
        protein: Protein,
        tracker: SimulationTracker,
        save: bool = True,
        filename: str = "max_compactness.png",
) -> None:
    """Plot the maximum-compactness fold recorded during the simulation."""
    fold = tracker.max_compactness_fold
    x = [fold[i][0] for i in range(protein.sequence_length)]
    y = [fold[i][1] for i in range(protein.sequence_length)]

    fig, ax = plt.subplots()
    ax.plot(x, y, alpha=0.5)
    for i in range(protein.sequence_length):
        ax.scatter(x[i], y[i], marker='$' + protein.sequence[i] + '$', s=20, color='red')

    ax.set_xlim(min(x) - 6, max(x) + 6)
    ax.set_ylim(min(y) - 6, max(y) + 6)
    ax.grid(alpha=0.2)
    ax.set_title('Max compactness structure')

    compactness = max(tracker.compactness_evolution)
    idx = tracker.compactness_evolution.index(compactness)
    energy = tracker.energy_evolution[idx]
    max_comp = compactness + 1e-14

    ax.text(0.01, 0.99, f'Energy: {energy}', ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.95, f'Compactness: {compactness / max_comp:.2f}', ha='left', va='top', transform=ax.transAxes)

    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


# ------------------------------------------------------------------
# Evolution plots
# ------------------------------------------------------------------

def plot_energy(
        tracker: SimulationTracker,
        avg: int = 10,
        save: bool = True,
        filename: str = "energy_evolution.png",
        annealing: bool = True,
) -> None:
    """Plot the energy evolution over simulation steps, averaged in windows of `avg`."""
    en_evo = np.array(tracker.energy_evolution[1:].copy())
    if len(en_evo) == 0:
        print("Energy evolution is empty, skipping plot.")
        return

    x = np.arange(0, len(en_evo), avg)

    fig, ax = plt.subplots()
    try:
        en_evo_avg = en_evo.reshape(-1, avg).mean(axis=1)
    except ValueError:
        print(f'Mean procedure skipped for energy: array not divisible by {avg}')
        en_evo_avg = en_evo
        x = np.arange(len(en_evo_avg))

    ax.set_title(f'Energy evolution averaged over {avg}-step windows')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Energy', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, en_evo_avg, color='b', label='Energy')

    min_idx = np.argmin(en_evo_avg)
    ax.plot(x[min_idx], en_evo_avg[min_idx], 'yo', markersize=4,
            label=f'Min energy = {en_evo_avg[min_idx]:.2f}')
    ax.plot(x[-1], en_evo_avg[-1], 'go', markersize=4,
            label=f'Final energy = {en_evo_avg[-1]:.2f}')
    ax.legend()

    if annealing:
        T = np.array(tracker.temperature_evolution[1:])
        try:
            T_avg = T.reshape(-1, avg).mean(axis=1)
        except ValueError:
            T_avg = T
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color='r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T_avg, color='r', label='Temperature')

    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


def plot_compactness(
        tracker: SimulationTracker,
        avg: int = 10,
        save: bool = True,
        filename: str = "compactness_evolution.png",
        annealing: bool = True,
) -> None:
    """Plot the normalized compactness evolution over simulation steps."""
    comp = np.array(tracker.compactness_evolution[1:].copy(), dtype=float)
    if len(comp) == 0:
        print("Compactness evolution is empty, skipping plot.")
        return

    c_max = comp.max()
    if c_max > 0:
        comp /= c_max
    
    x = np.arange(0, len(comp), avg)

    fig, ax = plt.subplots()
    try:
        comp_avg = comp.reshape(-1, avg).mean(axis=1)
    except ValueError:
        print(f'Mean procedure skipped for compactness: array not divisible by {avg}')
        comp_avg = comp
        x = np.arange(len(comp_avg))

    ax.set_title(f'Compactness evolution averaged over {avg}-step windows')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Compactness', color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.plot(x, comp_avg, color='b', label='Compactness')

    max_idx = np.argmax(comp_avg)
    ax.plot(x[max_idx], comp_avg[max_idx], 'yo', markersize=4,
            label=f'Max compactness = {comp_avg[max_idx]:.2f}')
    ax.plot(x[-1], comp_avg[-1], 'go', markersize=4,
            label=f'Final compactness = {comp_avg[-1]:.2f}')
    ax.legend()

    if annealing:
        T = np.array(tracker.temperature_evolution[1:])
        try:
            T_avg = T.reshape(-1, avg).mean(axis=1)
        except ValueError:
            T_avg = T
        ax_tw = ax.twinx()
        ax_tw.set_ylabel('T', color='r')
        ax_tw.tick_params(axis='y', labelcolor='r')
        ax_tw.plot(x, T_avg, color='r', label='Temperature')

    fig.tight_layout()
    plt.show(block=False)
    if save:
        plt.savefig(f"output/{filename}", format="png", bbox_inches="tight", dpi=200)


# ------------------------------------------------------------------
# GIF
# ------------------------------------------------------------------

def create_gif(
        protein: Protein,
        tracker: SimulationTracker,
        filename: str = "evolution.gif",
) -> None:
    """
    Create an animated GIF of accepted fold moves collected during the simulation.

    Parameters
    ----------
    protein : Protein
        Used to read the sequence for monomer labels.
    tracker : SimulationTracker
        Source of GIF frames (gif_struct, gif_steps) and total step count.
    filename : str
        Output GIF filename (saved inside output/).
    """
    n_steps = tracker.gif_steps[-1] if tracker.gif_steps else 0

    fig, ax = plt.subplots()
    writer = PillowWriter(fps=2)

    with writer.saving(fig, f'output/{filename}', 200):
        for frame, step_real in tqdm(
                zip(tracker.gif_struct, tracker.gif_steps),
                total=len(tracker.gif_struct),
                desc="Creating GIF",
        ):
            x = [coord[0] for coord in frame]
            y = [coord[1] for coord in frame]

            ax.plot(x, y, alpha=0.5)
            for j in range(protein.sequence_length):
                ax.scatter(x[j], y[j], marker='$' + protein.sequence[j] + '$', s=20, color='red')

            ax.set_xlim(min(x) - 6, max(x) + 6)
            ax.set_ylim(min(y) - 6, max(y) + 6)
            ax.grid(alpha=0.2)
            ax.set_title(f"Evolution step {step_real}/{n_steps}", fontsize=10)

            fig.canvas.draw()
            writer.grab_frame()
            ax.clear()
