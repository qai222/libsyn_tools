import functools
import json
import os
import random
import re
from typing import Union

from .chem import StateOfMatter

FilePath = Union[str, os.PathLike]


def json_parse_constant(arg):
    c = {"-Infinity": -float("inf"), "Infinity": float("inf"), "NaN": float("nan")}
    return c[arg]


def json_dump(o, fn: FilePath):
    with open(fn, "w") as f:
        json.dump(o, f, indent=2)


def json_load(fn: FilePath):
    with open(fn, "r") as f:
        o = json.load(f, parse_constant=json_parse_constant)
    return o


@functools.cache
def parse_chemscraper_output(scraper_output: FilePath) -> dict[str, tuple[str | None, float | None]]:
    """
    parse the output from [ChemScraper](https://github.com/qai222/ChemScraper)

    :param scraper_output: the json file from ChemScraper
    :return: d[<molecular smiles>] -> (state of matter, density), values are `None` if not found from the scraper
    """
    scraper_output = json_load(scraper_output)
    output_data = dict()
    for smi, data in scraper_output.items():
        try:
            form_string = data['form']
            if "liquid" in form_string.lower():
                form = StateOfMatter.LIQUID
            else:
                form = StateOfMatter.SOLID
        except KeyError:
            form = None
        try:
            density_string = data['density']
            try:
                density_string = re.findall("\d+\.\d+\s*g\/mL", density_string)[0]
                density_string = density_string.replace("g/mL", "").strip()
            except IndexError:
                density_string = re.findall("\d+\.\d+\s*g\/cm", density_string)[0]
                density_string = density_string.replace("g/cm", "").strip()
            density = float(density_string)
        except KeyError:
            density = None
        output_data[smi] = (form, density)
    return output_data


def parse_sparrow_routes(routes_file: FilePath, sample_seed: int = None, sample_n_target: int = None):
    """
    this method converts the "routes.json" file from sparrow/askcos to a list of reaction smiles

    :param sample_n_target: if randomly sample the target list
    :param sample_seed:
    :param routes_file: the file path
    :return: a list of reagent-free reaction smiles
    """
    routes = json_load(routes_file)
    routes = {k: routes[k] for k in routes if len(routes[k]['Reactions']) > 1}  # exclude orphans

    if sample_seed is None:
        sample_seed = 42

    if sample_n_target:
        random.seed(sample_seed)
        routes = {k: routes[k] for k in random.sample(sorted(routes.keys()), k=sample_n_target)}
    else:
        routes = {k: routes[k] for k in sorted(routes.keys())}

    reaction_smis = []
    for target, data in routes.items():
        for r in data['Reactions']:
            r_smi = r['smiles']
            if r_smi.startswith(">>"):
                continue
            reaction_smis.append(r_smi)
    reaction_smis = sorted(set(reaction_smis))
    return reaction_smis
