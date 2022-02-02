from rdkit.Chem.rdmolfiles import MolFromSmiles

from evomol.evaluation_entropy import extract_checkmol, extract_shingles
from evomol.molgraphops.molgraph import MolGraph

print(extract_checkmol(MolGraph(MolFromSmiles("CC(=O)CC(=O)C"))))

print(extract_shingles("CC(=O)CC(=O)C", 1))
