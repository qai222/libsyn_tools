import os

this_dir = os.path.dirname(os.path.abspath(__file__))

ROUTES_DIR = os.path.join(this_dir, 'example_routes_data')

ROUTES_PATH_FDA = os.path.join(ROUTES_DIR, 'routes_FDA.json')
"""
This case study select compounds from FDA approved drug lists of 2019-2022.
"""
# TODO add more info
assert os.path.isfile(ROUTES_PATH_FDA)

ROUTES_PATH_VS = os.path.join(ROUTES_DIR, 'routes_VS.json')
"""
35 of the 60 compounds from section _Selection and synthesis of candidate hits_ in 
[Sadybekov, Arman A., et al. "Synthon-based ligand discovery in virtual libraries of over 11 billion compounds." 
Nature 601.7893 (2022): 452-459.](https://www.nature.com/articles/s41586-021-04220-9).
"""
assert os.path.isfile(ROUTES_PATH_VS)
