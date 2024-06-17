from enum import Enum

import periodictable as pt
from pysmiles import read_smiles

"""
chemistry-related utility functions
"""


class StateOfMatter(str, Enum):
    """
    enum for allowed states of matter
    """
    LIQUID = "LIQUID"
    SOLID = "SOLID"
    GAS = "GAS"
    PLASMA = "PLASMA"
    MULTI = "MULTI"


def calculate_mw(smiles: str) -> float:
    """
    calculate the molecular mass for a given SMILES

    :param smiles: the given molecular SMILES
    :return: molecular mass in g/mol
    """
    mol = read_smiles(smiles, explicit_hydrogen=True)
    mw = 0
    for _, element in mol.nodes(data='element'):
        e = pt.elements.symbol(element)
        mw += e.mass
    return mw


def estimate_property_dummy(smiles: str) -> dict[str, float]:
    """
    a dummy estimator for molecular properties based on molecular mass

    :param smiles:
    :return: a dictionary of molecular properties
    """
    mw = calculate_mw(smiles)
    if mw < 200 and all(m not in smiles for m in ['Li', 'Na', 'K', 'Mg', 'Ca', 'Pd']):
        form = StateOfMatter.LIQUID
        density = 1.0
    else:
        form = StateOfMatter.SOLID
        density = 1.4
    return {"mw": mw, "state_of_matter": form, "density": density}
