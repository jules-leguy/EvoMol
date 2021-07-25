from rdkit.Chem.rdmolfiles import MolFromSmiles

from evomol import OPTEvaluationStrategy, LinearCombinationEvaluationStrategy, run_model, QEDEvaluationStrategy, \
    CLScoreEvaluationStrategy, SAScoreEvaluationStrategy, GaussianWrapperEvaluationStrategy
from evomol.evaluation_dft import SharedLastComputation
from evomol.molgraphops.molgraph import MolGraph

shared_last_computation = SharedLastComputation()

s_homo = OPTEvaluationStrategy("homo", shared_last_computation=shared_last_computation, MM_program="rdkit")
s_lumo = OPTEvaluationStrategy("lumo", shared_last_computation=shared_last_computation, MM_program="rdkit")
s_homo_m1 = OPTEvaluationStrategy("homo-1", shared_last_computation=shared_last_computation, MM_program="rdkit")

s_total = LinearCombinationEvaluationStrategy(
    evaluation_strategies=[
        GaussianWrapperEvaluationStrategy([s_homo], -4.5244369922635, 1),
        GaussianWrapperEvaluationStrategy([s_lumo], -1.004100108345, 1),
        GaussianWrapperEvaluationStrategy([s_homo_m1], -6.4466492321955, 1)
    ],
    coefs=[1, 1, 1]
)
#
# score, scores = s_total.evaluate_individual(MolGraph(MolFromSmiles("S1C=CSC1=C2SC=CS2")))
#
# print(s_total.keys())
# print(score)
# print(scores)

run_model(
    {
        "obj_function": s_total,
        "io_parameters": {
            "model_path": "22_test_shared_values",
            "save_n_steps": 1
        }
    }
)