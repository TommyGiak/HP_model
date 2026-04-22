"""
@author: Tommaso Giacometti, Alessandro Quirile
"""


def convert_to_hp(seq: str) -> str:
    """
    Transform a protein sequence of standard amino acids into the HP model sequence.

    Parameters
    ----------
    seq : str
        Sequence containing the 20 standard amino acids (RNDQEHKSTACGILMFPWYV),
        all must be uppercase letters.

    Returns
    -------
    str
        The sequence converted into H (hydrophobic) and P (polar).

    Raises
    ------
    ValueError
        If any amino acid in the sequence is not recognized.
    """
    polar = set('RNDQEHKST')  # polar amino acids
    hydrophobic = set('ACGILMFPWYV')  # hydrophobic amino acids

    hp_seq = []

    for amino_acid in seq:
        if amino_acid in polar:
            hp_seq.append('P')
        elif amino_acid in hydrophobic:
            hp_seq.append('H')
        else:
            raise ValueError(f"Amino acid '{amino_acid}' not recognized")

    return ''.join(hp_seq)
