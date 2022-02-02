from rdkit.Chem.rdmolfiles import MolFromSmiles

from evomol import OPTEvaluationStrategy
from evomol.evaluation import GaussianWrapperEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph

smi = "c1nncnn1"
smi = "C1=C2C=C3CC1N1C2C31"
# smi = "CO"

s = OPTEvaluationStrategy("gap")

score, scores = s.evaluate_individual(MolGraph(MolFromSmiles(smi)))
print(s.keys())
print(score)
print(scores)

s2 = GaussianWrapperEvaluationStrategy(mu=1.75, sigma=0.125, evaluation_strategies=[s])

score, scores = s2.evaluate_individual(MolGraph(MolFromSmiles(smi)))
print(s2.keys())
print(score)
print(scores)
