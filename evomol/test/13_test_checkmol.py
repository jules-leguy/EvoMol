import os

from evomol.evaluation_entropy import extract_checkmol, extract_shingles
from evomol.molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles


print(extract_checkmol(MolGraph(MolFromSmiles("CC(=O)CC(=O)C"))))

print(extract_shingles(MolGraph(MolFromSmiles("CC(=O)CC(=O)C")), 1))