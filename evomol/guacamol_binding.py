from os.path import join
from typing import Optional, List
import numpy as np
from .evaluation import EvaluationStrategy
from guacamol.goal_directed_benchmark import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction
from .stopcriterion import KthScoreMaxValue


class GuacamolEvaluationStrategy(EvaluationStrategy):

    def __init__(self, scoring_function, benchmark_name):
        self.scoring_function = scoring_function
        self.benchmark_name = benchmark_name

    def keys(self):
        return [self.benchmark_name]

    def evaluate_individual(self, individual):
        score = self.scoring_function.score(individual.to_aromatic_smiles())
        return score, [score]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score


class ChemPopAlgGoalDirectedGenerator(GoalDirectedGenerator):
    """
    Binding EvoMol population algorithm with the GuacaMol benchmark
    """

    def __init__(self, pop_alg, guacamol_init_top_100, init_pop_path, output_save_path):
        """
        :param pop_alg: PopAlg instance
        :param guacamol_init_top_100: whether the starting dataset must be the 100 best scoring molecules of the given
        dataset for each property
        :param init_pop_path: initial population path
        """
        self.pop_alg = pop_alg
        self.guacamol_init_top_100 = guacamol_init_top_100
        self.init_pop_path = init_pop_path
        self.output_save_path = output_save_path
        self.curr_benchmark_id = -1

    def _get_benchmark_name(self, curr_benchmark_id):
        """
        Returning the name of the current benchmark based on the number of calls to the generate_optimized_molecules
        method
        :param curr_benchmark_id:
        :return:
        """

        benchmark_names_list = ["Celecoxib", "Troglitazone", "Thiothixene", "Aripiprazole", "Albuterol", "Mestranol",
                                "C11H24", "C9H10N2O2PF2Cl", "Median molecules 1", "Median molecules 2",
                                "Osimertinib MPO", "Fexofenadine MPO", "Ranolazine MPO", "Perindopril MPO",
                                "Amlodipine MPO", "Sitagliptin MPO", "Zaleplon MPO", "Valsartan SMARTS", "Scaffold Hop",
                                "Deco Hop"]

        return benchmark_names_list[curr_benchmark_id]

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int,
                                     starting_population: Optional[List[str]] = None, name=None) -> List[str]:

        # Updating benchmark id
        self.curr_benchmark_id += 1

        # Extracting benchmark name
        curr_benchmark_name = self._get_benchmark_name(self.curr_benchmark_id)

        # Setting folder to save the results
        self.pop_alg.output_folder_path = join(self.output_save_path, curr_benchmark_name)

        # Setting GuacaMol evaluation function
        evaluation_strategy = GuacamolEvaluationStrategy(scoring_function, curr_benchmark_name)
        self.pop_alg.evaluation_strategy = evaluation_strategy
        self.pop_alg.mutation_strategy.evaluation_strategy = evaluation_strategy

        # Setting additional stop criterion, stopping the execution when best possible score is obtained
        additional_stop_criterion = KthScoreMaxValue(1, round=3)
        self.pop_alg.stop_criterion_strategy.set_additional_strategy(additional_stop_criterion)
        self.pop_alg.stop_criterion_strategy.set_pop_alg_instance(self.pop_alg)

        # PopAlg instance initialization
        self.pop_alg.initialize()

        # Population initialization
        if self.guacamol_init_top_100:

            # Extracting the top 100 SMILES for the property from ChEMBL and setting it at initial population
            # From https://github.com/BenevolentAI/guacamol_baselines/blob/master/graph_ga/goal_directed_generation.py
            with open(self.init_pop_path, "r") as f:

                smiles_list = f.readlines()
                scores = [scoring_function.score(s) for s in smiles_list]
                top_100_smiles = np.array(smiles_list)[np.argsort(scores)[::-1][:100]]
                self.pop_alg.load_pop_from_smiles_list(smiles_list=top_100_smiles)
        else:
            self.pop_alg.load_pop_from_smiles_list(smiles_list=["C"])

        # Running EvoMol
        self.pop_alg.run()

        # Extracting best individuals
        ind_to_return_indices = np.argsort(self.pop_alg.curr_total_scores)[::-1].flatten()[:number_molecules]
        output_population = []
        for ind_idx in ind_to_return_indices:
            output_population.append(self.pop_alg.pop[ind_idx].to_aromatic_smiles())

        # Returning optimized population
        return output_population
