from evomol import EntropyContribEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles

s = EntropyContribEvaluationStrategy(500, 3, descriptor_key="gen_scaffolds")

s.compute_record_scores_init_pop(
    population=[MolGraph(MolFromSmiles("CCC")),
                MolGraph(MolFromSmiles("CC")),
                MolGraph(MolFromSmiles("NN"))]
)

s.set_to_be_replaced_idx([1])

print(s.scores)
print(s.evaluate_individual_exact(MolGraph(MolFromSmiles("CC"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("CC"))))
print(s.evaluate_individual_exact(MolGraph(MolFromSmiles("CCCC"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("CCCC"))))

print()

s = EntropyContribEvaluationStrategy(500, 3, descriptor_key="gen_scaffolds")

s.compute_record_scores_init_pop(
    population=[MolGraph(MolFromSmiles("C")),
                MolGraph(MolFromSmiles("CC")),
                MolGraph(MolFromSmiles("CCC"))]
)

s.set_to_be_replaced_idx([1])

print(s.scores)
print(s.evaluate_individual_exact(MolGraph(MolFromSmiles("CC"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("CC"))))
print(s.evaluate_individual_exact(MolGraph(MolFromSmiles("CCCC"))))
print(s.evaluate_individual(MolGraph(MolFromSmiles("CCCC"))))
