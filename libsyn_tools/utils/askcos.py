from __future__ import annotations

import json

import requests

from .constants import ASKCOS_URL


def query_askcos_condition_rec(reaction_smiles: str, return_query: bool = False) -> dict | tuple[dict, dict]:
    q = dict(
        headers={"Content-type": "application/json"},
        url=f"{ASKCOS_URL}/api/context-recommender/v2/predict/FP/call-sync",
        data=f'{{"smiles": "{reaction_smiles}", "n_conditions": 3}}',
    )
    response = requests.post(**q)
    response = response.content.decode('utf8')
    response = json.loads(response)
    assert response['status_code'] == 200
    if return_query:
        return response, q
    else:
        return response


def drawing_url(smiles: str, size=80) -> str:
    smiles = smiles.replace("+", "%2b")
    q = f"/api/draw/?smiles={smiles}&transparent=true&annotate=false&size={size}"
    q = ASKCOS_URL + q
    return q
