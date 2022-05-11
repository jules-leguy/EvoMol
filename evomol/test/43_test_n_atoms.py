from rdkit.Chem import MolFromSmiles

from evomol.molgraphops.molgraph import MolGraph

print(MolGraph(MolFromSmiles("CCl")).get_n_atoms())
