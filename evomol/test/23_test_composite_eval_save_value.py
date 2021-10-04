from rdkit.Chem import MolFromSmiles

from evomol import OppositeWrapperEvaluationStrategy, MeanEvaluationStrategyComposite, CLScoreEvaluationStrategy, \
    QEDEvaluationStrategy, GaussianWrapperEvaluationStrategy, LinearCombinationEvaluationStrategy, \
    SAScoreEvaluationStrategy, PenalizedLogPEvaluationStrategy, run_model
from evomol.molgraphops.molgraph import MolGraph

eval_strat = MeanEvaluationStrategyComposite([
    OppositeWrapperEvaluationStrategy([
        LinearCombinationEvaluationStrategy([
            PenalizedLogPEvaluationStrategy(),
            GaussianWrapperEvaluationStrategy([SAScoreEvaluationStrategy()], 3, 1, True)
        ], coefs=[0.5, 0.5]),
    ]),
    QEDEvaluationStrategy()
])

pop = []
eval_strat.compute_record_scores_init_pop(pop)

ind = MolGraph(MolFromSmiles("CN"))
ind2 = MolGraph(MolFromSmiles(("CONF")))
ind3 = MolGraph(MolFromSmiles(("NC1CNNNNNNNNNC(O)NNNNNCC[SH]=[SH]NNNNONNNNNNNNNN1")))

print(eval_strat.keys())
eval1 = eval_strat.evaluate_individual(ind)
print(eval1)
eval2 = eval_strat.evaluate_individual(ind2)
print(eval2)
eval3 = eval_strat.evaluate_individual(ind3)
print(eval3)

eval_strat.record_ind_score(0, eval1[0], eval1[1], ind)
print(eval_strat.get_population_scores())
eval_strat.record_ind_score(1, eval2[0], eval2[1], ind2)
print(eval_strat.get_population_scores())
eval_strat.record_ind_score(2, eval3[0], eval3[1], ind2)
print(eval_strat.get_population_scores())

run_model({
    "obj_function": eval_strat
})