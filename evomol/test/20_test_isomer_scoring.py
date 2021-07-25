from evomol.evaluation import IsomerGuacaMolEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles

s = IsomerGuacaMolEvaluationStrategy("C9H8O4")

print(s.evaluate_individual(MolGraph(MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("CC(=O)NC1=CC=C(C=C1)O"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("Ck"))))