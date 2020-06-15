from abc import ABC, abstractmethod

import numpy as np
from .evaluation import EvaluationError
from .molgraphops.molgraph import MolGraphBuilder
from .molgraphops.exploration import random_neighbour


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
    def mutate(self, individual, ind_to_replace, curr_score, pop_tabu_list):
        """
        Finding an improver of ind_to_replace by mutating individual. The improver cannot be already in the population.
        :param individual: individual to be mutated to find an improver of ind_to_replace
        :param ind_to_replace: individual to be replaced
        :param curr_score: score of the individual to be replaced
        :param pop_tabu_list: list of individuals currently in the population
        :return: (improver, mutation string description, total score of the new individual, list of intermediate scores)
        """


class KRandomGraphOpsImprovingMutationStrategy(MutationStrategy):
    """
    Performing a graph operations mutation composed of at most k actions. The mutation is performed at most n_max_try
    to find an improver. If no improver is found, raising a MutationError.
    """

    def __init__(self, k, max_n_try, evaluation_strategy, action_spaces, action_spaces_parameters, problem_type="max"):
        """

        :param k: max number of successive graph operations
        :param max_n_try: max number of tries to find an improver
        :param evaluation_strategy: EvaluationStrategy instance with an evaluate_individual method
        :param action_spaces: list of ActionSpace instances
        :param action_spaces_parameters: instance of ActionSpace.ActionSpaceParameters
        :param plateau: whether individuals of same score are considered improvers
        :param uniform_action_type: If true, the action type is drawn with a uniform law before the action is drawn. If
        false, the action is drawn directly with uniform law among all possible actions
        :param problem_type: Whether it is a maximization ("max") or a minimization ("min") problem
        """
        self.k = k
        self.max_n_try = max_n_try
        self.evaluation_strategy = evaluation_strategy
        self.action_spaces = action_spaces
        self.actionspace_parameters = action_spaces_parameters
        self.problem_type = problem_type

    def is_improver(self, curr_total_score, mutated_total_score):

        return curr_total_score is None \
               or (mutated_total_score > curr_total_score and self.problem_type == "max") \
               or (mutated_total_score < curr_total_score and self.problem_type == "min") \
               or (mutated_total_score == curr_total_score)

    def mutate(self, individual, ind_to_replace, curr_total_score, pop_tabu_list):

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

            # Only evaluating the neighbour if it has not been encountered yet in the population
            if not mutated_ind.to_aromatic_smiles() in pop_tabu_list:

                try:

                    # Computing score
                    mutated_total_score, mutated_scores = self.evaluation_strategy.evaluate_individual(mutated_ind)

                except Exception as e:
                    raise EvaluationError(str(e) + individual.to_aromatic_smiles() + " " + desc) from e

                # Returning mutated individual if it is an improver
                if self.is_improver(curr_total_score, mutated_total_score):
                    return mutated_ind, desc, mutated_total_score, mutated_scores

        # Raising error if no improver was found
        raise NoImproverError("No improver found")
