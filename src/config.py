import yaml
import random


class Configuration:
    """
    Class to parse and store parameters from a YAML configuration file.
    """

    def __init__(self, filename: str) -> None:
        with open(filename, 'r') as f:
            config = yaml.safe_load(f)

        # Sequence
        self.sequence = config['sequence'].replace(' ', '').replace('\n', '')

        # Structure options
        structure_cfg = config.get('structure', {})
        self.use_struct = structure_cfg.get('use_structure')
        self.fold = structure_cfg.get('coordinates') if self.use_struct else None

        # Simulation options
        sim_cfg = config.get('simulation', {})
        self.n_steps = sim_cfg.get('folding_steps')
        self.do_annealing = sim_cfg.get('annealing')
        self.temperature = sim_cfg.get('temperature')

        # Plot options
        plot_cfg = config.get('plot', {})
        self.do_gif = plot_cfg.get('create_gif')

        # Seed
        seed_val = config.get('seed', None)
        if seed_val is None or seed_val == 'None':
            seed_val = random.randint(0, 10000)
        self.seed = int(seed_val)

    def print_simulation_setup(self):
        """Display the main simulation configuration."""
        print("[HP Model] Simulation setup")
        print(f"Sequence:           {self.sequence}")
        print(f"Sequence length:    {len(self.sequence)}")
        print(f"Structure:          {'Non-linear' if self.use_struct else 'Linear'}")
        print(f"Folding steps:      {self.n_steps}")
        print(f"Annealing:          {self.do_annealing}")
        print(f"Temperature:        {self.temperature}")
        print(f"Seed:               {self.seed}")
