from evomol import OPTEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles

s_obabel = OPTEvaluationStrategy("gap", n_jobs=2, working_dir_path="/home/jleguy/dft_comput",
                                 MM_program="obabel")

s_rdkit = OPTEvaluationStrategy("gap", n_jobs=2, working_dir_path="/home/jleguy/dft_comput",
                                 MM_program="rdkit")

smiles = ["C", "S=CC=S"]

result_obabel = []
result_rdkit = []

for smi in smiles:

    result_obabel.append(s_obabel.evaluate_individual(MolGraph(MolFromSmiles(smi))))
    result_rdkit.append(s_rdkit.evaluate_individual(MolGraph(MolFromSmiles(smi))))


print(result_obabel)
print(result_rdkit)