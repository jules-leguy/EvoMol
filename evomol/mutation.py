import time
from abc import ABC, abstractmethod

import numpy as np

from .evaluation import EvaluationError, RDFiltersEvaluationStrategy, SillyWalksEvaluationStrategy, \
    SAScoreEvaluationStrategy
from .molgraphops.exploration import random_neighbour
from .molgraphops.molgraph import MolGraphBuilder


class MutationError(RuntimeError):
    """
    Exception raised when the mutation fails.
    """

    def __init__(self, desc):
        self.desc = desc

    def __str__(self):
        return self.desc + " (Mutation error)"


class NoImproverError(RuntimeError):
    pass


class MutationStrategy(ABC):
    """
    Interface of the class used to defined the mutations to apply to the individuals of the population.
    """

    @abstractmethod
    def mutate(self, individual, ind_to_replace, curr_score, pop_tabu_list, external_tabu_list,
               generated_ind_recorder):
        """
        Finding an improver of ind_to_replace by mutating individual. The improver cannot be already in the population.
        :param individual: individual to be mutated to find an improver of ind_to_replace
        :param ind_to_replace: individual to be replaced
        :param curr_score: score of the individual to be replaced
        :param pop_tabu_list: list of individuals currently in the population
        :param external_tabu_list: list of tabu SMILES
        :param generated_ind_recorder: instance of PopAlg.GeneratedIndividualsRecorder in which all generated
        individuals are stored
        :return: (improver, mutation string description, total score of the new individual, list of intermediate scores)
        """


class KRandomGraphOpsImprovingMutationStrategy(MutationStrategy):
    """
    Performing a graph operations mutation composed of at most k actions. The mutation is performed at most n_max_try
    to find an improver. If no improver is found, raising a MutationError.
    """

    def __init__(self, k, max_n_try, evaluation_strategy, action_spaces, action_spaces_parameters, problem_type="max",
                 quality_filter=False, silly_molecules_fp_threshold=1, silly_molecules_db=None,
                 sascore_threshold=float("inf"), custom_filter_function=None):
        """

        :param k: max number of successive graph operations
        :param max_n_try: max number of tries to find an improver
        :param evaluation_strategy: EvaluationStrategy instance with an evaluate_individual method
        :param action_spaces: list of ActionSpace instances
        :param action_spaces_parameters: instance of ActionSpace.ActionSpaceParameters
        :param problem_type: Whether it is a maximization ("max") or a minimization ("min") problem
        :param quality_filter: Whether to prevent molecules that do not pass the quality filter to be considered as
        valid improvers (using https://github.com/PatWalters/rd_filters filters)
        :param silly_molecules_fp_threshold: Using Patrick Walters's silly walk program to count the proportion of
        bits in the ECFP4 fingerprint that do not exist in the ChemBL. The molecules with a proportion that is higher
        than the given threshold are discarded (https://github.com/PatWalters/silly_walks).
        :param sascore_threshold : discarding solutions that have a SAScore value above the given threshold
        :param custom_filter_function: custom boolean function that evaluates a smiles and assess whether it is part
        of an acceptable search space or not.
        """
        self.k = k
        self.max_n_try = max_n_try
        self.evaluation_strategy = evaluation_strategy
        self.action_spaces = action_spaces
        self.actionspace_parameters = action_spaces_parameters
        self.problem_type = problem_type
        self.quality_filter = quality_filter
        self.silly_molecules_fp_threshold = silly_molecules_fp_threshold
        self.sascore_threshold = sascore_threshold
        self.custom_filter_function = custom_filter_function

        if self.quality_filter:
            self.rd_filter_eval_strat = RDFiltersEvaluationStrategy()

        if self.silly_molecules_fp_threshold < 1:
            self.silly_walks_eval_strat = SillyWalksEvaluationStrategy(silly_molecules_db)

        if self.sascore_threshold < float("inf"):
            self.sascore_eval_strat = SAScoreEvaluationStrategy()

    def is_improver(self, curr_total_score, mutated_total_score):

        return curr_total_score is None \
               or (mutated_total_score > curr_total_score and self.problem_type == "max") \
               or (mutated_total_score < curr_total_score and self.problem_type == "min") \
               or (mutated_total_score == curr_total_score)

    def mutate(self, individual, ind_to_replace_idx, curr_total_score, pop_tabu_list, external_tabu_list,
               generated_ind_recorder):

        # Drawing the number of actions
        n_actions = int(np.random.choice(np.arange(1, self.k + 1)))

        # Trying max_n_try times to find an improver
        for i in range(self.max_n_try):

            try:

                # Creating QuMolGraphBuilder
                qumol_builder = MolGraphBuilder(self.actionspace_parameters, self.action_spaces, individual)

                # Performing mutation
                mutated_ind, desc = random_neighbour(qumol_builder, n_actions, return_mol_graph=True,
                                                     uniform_action_type=True)

            except Exception as e:
                print(e)
                raise MutationError(individual.to_aromatic_smiles()) from e

            # Computing boolean filter values
            failed_tabu_pop = mutated_ind.to_aromatic_smiles() in pop_tabu_list
            failed_tabu_external = external_tabu_list is not None and \
                                   mutated_ind.to_aromatic_smiles() in external_tabu_list
            failed_quality_filter = self.quality_filter and \
                                    self.rd_filter_eval_strat.evaluate_individual(mutated_ind)[0] == 0
            failed_sillywalks_filter = self.silly_molecules_fp_threshold < 1 and \
                                       self.silly_walks_eval_strat.evaluate_individual(mutated_ind)[
                                           0] > self.silly_molecules_fp_threshold
            failed_custom_filter = self.custom_filter_function is not None and \
                                   not self.custom_filter_function(mutated_ind.to_aromatic_smiles())

            try:
                failed_sascore_filter = self.sascore_threshold < float("inf") and \
                                        self.sascore_eval_strat.evaluate_individual(mutated_ind)[
                                            0] > self.sascore_threshold
            except ZeroDivisionError:
                failed_sascore_filter = True

            # Only evaluating the neighbour if it has not been encountered yet in the population and if it is valid
            # if the filters are enabled and if it is not in the external tabu list
            if not failed_tabu_pop and not failed_tabu_external and not failed_quality_filter and \
                    not failed_sillywalks_filter and not failed_sascore_filter and not failed_custom_filter:

                try:
                    tstart = time.time()

                    # Computing score
                    mutated_total_score, mutated_scores = \
                        self.evaluation_strategy.evaluate_individual(mutated_ind, to_replace_idx=ind_to_replace_idx)

                    evaluation_time = time.time() - tstart

                except Exception as e:

                    evaluation_time = time.time() - tstart

                    generated_ind_recorder.record_individual(individual=mutated_ind,
                                                             total_score=None,
                                                             scores=np.full((len(self.evaluation_strategy.keys(), )),
                                                                            None),
                                                             objective_calls=self.evaluation_strategy.n_calls,
                                                             success_obj_computation=False,
                                                             improver=False,
                                                             obj_computation_time=evaluation_time)

                    raise EvaluationError(str(e) + individual.to_aromatic_smiles() + " " + desc) from e

                # Checking if the mutated individual is an improver
                is_improver = self.is_improver(curr_total_score, mutated_total_score)

                # Recording the mutated individual
                generated_ind_recorder.record_individual(individual=mutated_ind,
                                                         total_score=mutated_total_score,
                                                         scores=mutated_scores,
                                                         objective_calls=self.evaluation_strategy.n_calls,
                                                         success_obj_computation=True,
                                                         improver=is_improver,
                                                         obj_computation_time=evaluation_time)

                # Returning the mutated solution if it is an improver
                if is_improver:
                    return mutated_ind, desc, mutated_total_score, mutated_scores, evaluation_time

            # The mutated solution did not pass the filters
            else:
                generated_ind_recorder.record_individual(individual=mutated_ind,
                                                         total_score=None,
                                                         scores=np.full((len(self.evaluation_strategy.keys(), )), None),
                                                         objective_calls=self.evaluation_strategy.n_calls,
                                                         success_obj_computation=False,
                                                         improver=False,
                                                         failed_tabu_pop=failed_tabu_pop,
                                                         failed_tabu_external=failed_tabu_external,
                                                         failed_rdfilters=failed_quality_filter,
                                                         failed_sillywalks=failed_sillywalks_filter,
                                                         failed_sascore=failed_sascore_filter,
                                                         failed_custom_filter=failed_custom_filter,
                                                         obj_computation_time=None)

        # Raising error if no improver was found
        raise NoImproverError("No improver found")
