from os.path import join
from typing import Optional, List
import numpy as np

from .evaluation import EvaluationStrategy, EvaluationStrategyComposite
from guacamol.goal_directed_benchmark import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction
from .stopcriterion import KthScoreMaxValue


class UndefinedGuacaMolEvaluationStrategy(EvaluationStrategy):
    """
    Class representing GuacaMol evaluation strategy not defined yet
    """

    def __init__(self, name):
        super().__init__()
        self.name = name

    def keys(self):
        return []

    def evaluate_individual(self, individual, to_replace_idx=None):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")

    def compute_record_scores_init_pop(self, population):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")

    def record_ind_score(self, idx, new_total_score, new_scores, new_individual):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")

    def get_population_scores(self):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")

    def end_step_population(self, pop):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")

    def get_additional_population_scores(self):
        raise RuntimeError("Undefined GuacaMol evaluation strategy being used")


class GuacamolEvaluationStrategy(EvaluationStrategy):

    def __init__(self, scoring_function, benchmark_name):
        super().__init__()
        self.scoring_function = scoring_function
        self.benchmark_name = benchmark_name

    def keys(self):
        return [self.benchmark_name]

    def evaluate_individual(self, individual, to_replace_idx=None):
        super().evaluate_individual(individual, to_replace_idx)
        score = self.scoring_function.score(individual.to_aromatic_smiles())
        return score, [score]


def define_GuacaMol_evaluation_strategies(evaluation_strategy, defined_evaluation_strategy):
    """
    Setting all UndefinedGuacaMolEvaluationStrategie instances contained in a EvaluationStrategyComposite to the given
    defined value
    """

    if isinstance(evaluation_strategy, EvaluationStrategyComposite):
        for i in range(len(evaluation_strategy.evaluation_strategies)):
            if isinstance(evaluation_strategy.evaluation_strategies[i], UndefinedGuacaMolEvaluationStrategy):
                evaluation_strategy.evaluation_strategies[i] = defined_evaluation_strategy
            elif isinstance(evaluation_strategy.evaluation_strategies[i], EvaluationStrategyComposite):
                define_GuacaMol_evaluation_strategies(evaluation_strategy.evaluation_strategies[i],
                                                      defined_evaluation_strategy)


def is_or_contains_undefined_GuacaMol_evaluation_strategy(evaluation_strategy):
    """
    Returns whether the given evaluation strategy is or contain an undefined GuacaMol evaluation strategy
    """

    if isinstance(evaluation_strategy, UndefinedGuacaMolEvaluationStrategy):
        return True

    elif isinstance(evaluation_strategy, EvaluationStrategyComposite):
        for i in range(len(evaluation_strategy.evaluation_strategies)):
            if is_or_contains_undefined_GuacaMol_evaluation_strategy(evaluation_strategy.evaluation_strategies[i]):
                return True

    return False


def get_GuacaMol_benchmark_parameter(evaluation_strategy):
    """
    Returning the GuacaMol benchmark parameter describing the set of benchmarks found in the first found
    UndefinedGuacaMolEvaluationStrategy
    """

    if isinstance(evaluation_strategy, UndefinedGuacaMolEvaluationStrategy):
        return evaluation_strategy.name
    elif isinstance(evaluation_strategy, EvaluationStrategyComposite):
        for i in range(len(evaluation_strategy.evaluation_strategies)):
            retrieved_value = get_GuacaMol_benchmark_parameter(evaluation_strategy.evaluation_strategies[i])
            if retrieved_value is not None:
                return retrieved_value

    return None


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

    def _generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int, name: str,
                                      starting_population: Optional[List[str]] = None) -> List[str]:
        """
        Identical to self.generate_optimized_molecules but when using the modified GuacaMol version that allows for
        parallel optimization of objectives.
        """
        instance = self.pop_alg.copy_instance_with_parameters()

        # Updating benchmark id
        self.curr_benchmark_id += 1

        # Extracting benchmark name
        curr_benchmark_name = self._get_benchmark_name(self.curr_benchmark_id)

        # Setting folder to save the results
        instance.output_folder_path = join(self.output_save_path, name)

        # Extracting GuacaMol evaluation function
        guacamol_evaluation_strategy = GuacamolEvaluationStrategy(scoring_function, name)

        # Merging the evaluation strategy of the PopAlg instance to the GuacaMol objective
        if isinstance(instance.evaluation_strategy, UndefinedGuacaMolEvaluationStrategy):
            instance.evaluation_strategy = guacamol_evaluation_strategy
        else:
            define_GuacaMol_evaluation_strategies(instance.evaluation_strategy, guacamol_evaluation_strategy)

        # Updating mutation strategy evaluator
        instance.mutation_strategy.evaluation_strategy = instance.evaluation_strategy

        # Setting additional stop criterion, stopping the execution when best possible score is obtained
        instance.kth_score_to_record_key = curr_benchmark_name
        # instance.kth_score_to_record_key = name
        additional_stop_criterion = KthScoreMaxValue(1, round=3)
        instance.stop_criterion_strategy.set_additional_strategy(additional_stop_criterion)
        instance.stop_criterion_strategy.set_pop_alg_instance(instance)

        # Setting kth score to record
        instance.kth_score_to_record = number_molecules

        # PopAlg instance initialization
        instance.initialize()

        # Population initialization
        if self.guacamol_init_top_100:

            # Extracting the top 100 SMILES for the property from ChEMBL and setting it at initial population
            # From https://github.com/BenevolentAI/guacamol_baselines/blob/master/graph_ga/goal_directed_generation.py
            with open(self.init_pop_path, "r") as f:

                smiles_list = f.readlines()
                scores = [scoring_function.score(s) for s in smiles_list]
                top_100_smiles = np.array(smiles_list)[np.argsort(scores)[::-1][:100]]
                instance.load_pop_from_smiles_list(smiles_list=top_100_smiles)
        else:
            instance.load_pop_from_smiles_list(smiles_list=["C"])

        # Running EvoMol
        instance.run()

        # Extracting the vector containing the guacamol  objective property value for all individuals
        if instance.kth_score_to_record_key == "total":
            obj_prop_vector = instance.curr_total_scores
        else:
            obj_prop_vector = instance.curr_scores[instance.kth_score_to_record_idx]

        # Extracting best individuals
        ind_to_return_indices = np.argsort(obj_prop_vector)[::-1].flatten()[:number_molecules]
        output_population = []
        for ind_idx in ind_to_return_indices:
            output_population.append(instance.pop[ind_idx].to_aromatic_smiles())

        # Returning optimized population
        return output_population

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int,
                                     starting_population: Optional[List[str]] = None) -> List[str]:

        instance = self.pop_alg.copy_instance_with_parameters()

        # Updating benchmark id
        self.curr_benchmark_id += 1

        # Extracting benchmark name
        curr_benchmark_name = self._get_benchmark_name(self.curr_benchmark_id)

        # Setting folder to save the results
        instance.output_folder_path = join(self.output_save_path, curr_benchmark_name)

        # Extracting GuacaMol evaluation function
        guacamol_evaluation_strategy = GuacamolEvaluationStrategy(scoring_function, curr_benchmark_name)

        # Merging the evaluation strategy of the PopAlg instance to the GuacaMol objective
        if isinstance(instance.evaluation_strategy, UndefinedGuacaMolEvaluationStrategy):
            instance.evaluation_strategy = guacamol_evaluation_strategy
        else:
            define_GuacaMol_evaluation_strategies(instance.evaluation_strategy, guacamol_evaluation_strategy)

        # Updating mutation strategy evaluator
        instance.mutation_strategy.evaluation_strategy = instance.evaluation_strategy

        # Setting additional stop criterion, stopping the execution when best possible score is obtained
        instance.kth_score_to_record_key = curr_benchmark_name
        # instance.kth_score_to_record_key = name
        additional_stop_criterion = KthScoreMaxValue(1, round=3)
        instance.stop_criterion_strategy.set_additional_strategy(additional_stop_criterion)
        instance.stop_criterion_strategy.set_pop_alg_instance(instance)

        # Setting kth score to record
        instance.kth_score_to_record = number_molecules

        # PopAlg instance initialization
        instance.initialize()

        # Population initialization
        if self.guacamol_init_top_100:

            # Extracting the top 100 SMILES for the property from ChEMBL and setting it at initial population
            # From https://github.com/BenevolentAI/guacamol_baselines/blob/master/graph_ga/goal_directed_generation.py
            with open(self.init_pop_path, "r") as f:

                smiles_list = f.readlines()
                scores = [scoring_function.score(s) for s in smiles_list]
                top_100_smiles = np.array(smiles_list)[np.argsort(scores)[::-1][:100]]
                instance.load_pop_from_smiles_list(smiles_list=top_100_smiles)
        else:
            instance.load_pop_from_smiles_list(smiles_list=["C"])

        # Running EvoMol
        instance.run()

        # Extracting the vector containing the guacamol  objective property value for all individuals
        if instance.kth_score_to_record_key == "total":
            obj_prop_vector = instance.curr_total_scores
        else:
            obj_prop_vector = instance.curr_scores[instance.kth_score_to_record_idx]

        # Extracting best individuals
        ind_to_return_indices = np.argsort(obj_prop_vector)[::-1].flatten()[:number_molecules]
        output_population = []
        for ind_idx in ind_to_return_indices:
            output_population.append(instance.pop[ind_idx].to_aromatic_smiles())

        # Returning optimized population
        return output_population
