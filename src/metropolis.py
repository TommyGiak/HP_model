"""Metropolis algorithm for stochastic optimization."""
import math
import random


def metropolis_rule(current_energy: float, new_energy: float, temperature: float) -> bool:
    """Decide if a new state should be accepted."""
    if new_energy < current_energy:
        return True
    if temperature <= 0:
        return False
    delta_e = new_energy - current_energy
    prob = math.exp(-delta_e / temperature)
    return prob >= random.uniform(0, 1)
