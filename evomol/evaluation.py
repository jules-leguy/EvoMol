import pickle
from abc import ABC, abstractmethod
from math import exp
from os.path import join

import networkx as nx
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.QED import qed
from rdkit.Chem.rdmolfiles import MolToSmiles, MolFromSmiles

from rdkit.Chem import RDConfig
import os
import sys

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
import numpy as np


class EvaluationError(RuntimeError):

    def __init__(self, desc):
        self.desc = desc

    def __str__(self):
        return self.desc + " (Evaluation error)"


class EvaluationStrategy(ABC):
    """
    Strategy to evaluate the individuals of a population.
    It can be single objective or multi-objective.
    """

    @abstractmethod
    def keys(self):
        """
        Returning a unique list of key(s) describing the evaluator(s)
        :return list of string key(s)
        """
        pass

    @abstractmethod
    def get_population_scores(self):
        """
        Returning the scores of the complete population.
        :return total scores of the population (list), list of list of intermediate scores for each contained evaluator
        """
        pass

    @abstractmethod
    def evaluate_individual(self, individual):
        """
        Evaluation of a given individual.
        :param individual: individual to be evaluated.
        :return: total score of given individual, list of intermediate scores for each contained evaluator
        """
        pass

    @abstractmethod
    def compute_record_scores(self, population):
        """
        Computing and recording the internal scores for the complete population
        :return: None
        """
        pass

    @abstractmethod
    def record_score(self, idx, new_total_score, new_scores):
        """
        Updating the scores of the individual at the given index
        :param idx: index
        :param new_total_score: total score of the individual
        :param new_scores: intermediate scores
        :return:
        """


class GenericFunctionEvaluationStrategy(EvaluationStrategy):
    """
    Evaluating individuals with a given function evaluating a SMILES representation of a molecule
    """

    def __init__(self, evaluation_function, function_name="custom_function"):
        """
        :param evaluation_function: must evaluate a SMILES representation with a value
        :param function_name: name of the function
        """
        self.evaluation_function = evaluation_function
        self.function_name = function_name

    def keys(self):
        return [self.function_name]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def evaluate_individual(self, individual):
        score = self.evaluation_function(individual.to_aromatic_smiles())
        return score, [score]


class ZincNormalizedPLogPEvaluationStrategy(EvaluationStrategy):

    def __init__(self):
        self.scores = None

    def keys(self):
        return ["penalized_logP"]

    def evaluate_individual(self, individual):
        """
        from https://github.com/bowenliu16/rl_graph_generation/blob/master/gym-molecule/gym_molecule/envs/molecule.py
        """
        # normalization constants, statistics from 250k_rndm_zinc_drugs_clean.smi
        logP_mean = 2.4570953396190123
        logP_std = 1.434324401111988
        SA_mean = -3.0525811293166134
        SA_std = 0.8335207024513095
        cycle_mean = -0.0485696876403053
        cycle_std = 0.2860212110245455

        mol_graph = MolFromSmiles(individual.to_aromatic_smiles())

        log_p = Descriptors.MolLogP(mol_graph)
        SA = -sascorer.calculateScore(mol_graph)

        # cycle score
        cycle_list = nx.cycle_basis(nx.Graph(
            Chem.rdmolops.GetAdjacencyMatrix(mol_graph)))
        if len(cycle_list) == 0:
            cycle_length = 0
        else:
            cycle_length = max([len(j) for j in cycle_list])
        if cycle_length <= 6:
            cycle_length = 0
        else:
            cycle_length = cycle_length - 6
        cycle_score = -cycle_length

        normalized_log_p = (log_p - logP_mean) / logP_std
        normalized_SA = (SA - SA_mean) / SA_std
        normalized_cycle = (cycle_score - cycle_mean) / cycle_std

        score = normalized_log_p + normalized_SA + normalized_cycle

        return score, [score]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score

    def get_population_scores(self):
        return self.scores, np.array([self.scores])


class PenalizedLogPEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of penalized logP
        from : https://github.com/google-research/google-research/blob/master/mol_dqn/chemgraph/dqn/py/molecules.py
    """

    def __init__(self):
        self.scores = None

    def keys(self):
        return ["penalized_logP"]

    def get_largest_ring_size(self, molecule):
        """
        Calculates the largest ring size in the molecule.
        Refactored from
        https://github.com/wengong-jin/icml18-jtnn/blob/master/bo/run_bo.py
        Args:
          molecule: Chem.Mol. A molecule.
        Returns:
          Integer. The largest ring size.
        """
        cycle_list = molecule.GetRingInfo().AtomRings()
        if cycle_list:
            cycle_length = max([len(j) for j in cycle_list])
        else:
            cycle_length = 0
        return cycle_length

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def evaluate_individual(self, individual):

        mol_graph = MolFromSmiles(individual.to_aromatic_smiles())

        log_p = Descriptors.MolLogP(mol_graph)
        sas_score = sascorer.calculateScore(mol_graph)
        largest_ring_size = self.get_largest_ring_size(mol_graph)
        cycle_score = max(largest_ring_size - 6, 0)
        score = log_p - sas_score - cycle_score
        return score, [score]

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score

    def get_population_scores(self):
        return self.scores, np.array([self.scores])


class CLScoreEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of CLscore (Bühlmann, Sven, et Jean-Louis Reymond. « ChEMBL-Likeness Score and Database GDBChEMBL ».
    Frontiers in Chemistry 8 (4 février 2020). https://doi.org/10.3389/fchem.2020.00046.)
    Based on https://github.com/reymond-group/GDBChEMBL
    """

    def __init__(self):
        self.scores = None
        self.radius = 3
        self.rooted = True
        self.weighted = True
        self.cut_off = 0.0

        # Loading ChEMBL shingles database
        if self.rooted:
            with open(join(os.environ["SHINGLE_LIBS"],
                           "chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl"), "rb") as pyc:
                self.db_shingles = pickle.load(pyc)
        else:
            with open(join(os.environ["SHINGLE_LIBS"],
                           "chembl_24_1_shingle_scores_log10_nrooted_nchir.pkl"), "rb") as pyc:
                self.db_shingles = pickle.load(pyc)

    def keys(self):
        return ["CLScore"]

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def extract_shingles(self, individual):

        qry_shingles = set()

        radius_constr = self.radius + 1

        # Reloading molecule to make it aromatic
        mol = MolFromSmiles(individual.to_aromatic_smiles())

        for atm_idx in range(individual.mol_graph.GetNumAtoms()):
            for N in range(1, radius_constr):
                bonds = AllChem.FindAtomEnvironmentOfRadiusN(mol, N, atm_idx)

                if not bonds:
                    break

                # the reportedly faster method
                atoms = set()
                for bond_id in bonds:
                    bond = mol.GetBondWithIdx(bond_id)
                    atoms.add(bond.GetBeginAtomIdx())
                    atoms.add(bond.GetEndAtomIdx())

                if self.rooted:
                    new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol, list(atoms), bonds, 0, 0,
                                                                      False, False, atm_idx, True, False, False)
                else:
                    new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol, list(atoms), bonds, 0, 0,
                                                                      False, False, -1, True, False, False)

                qry_shingles.add(new_shingle)

        return qry_shingles

    def evaluate_individual(self, individual):
        """
        Based on https://github.com/reymond-group/GDBChEMBL
        :param individual:
        :return:
        """

        # Extracting shingles
        qry_shingles = self.extract_shingles(individual)

        # calculate shingle count averaged score
        avg_score = 0
        if qry_shingles:
            sum_scores = 0
            # using log10 of shingle frequency
            if self.weighted:
                for shingle in qry_shingles:
                    # if key not present, add 0 per default
                    sum_scores += self.db_shingles.get(shingle, 0)
            # working binary (i.e. if present -> count++ )
            else:
                for shingle in qry_shingles:
                    if shingle in self.db_shingles:
                        sum_scores += 1
            avg_score = sum_scores / len(qry_shingles)

        if self.cut_off == 0.0 or self.cut_off <= avg_score:
            return avg_score, [avg_score]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score


class SAScoreEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of SAScore.
    Ertl, Peter, et Ansgar Schuffenhauer. « Estimation of synthetic accessibility score of drug-like molecules based
    on molecular complexity and fragment contributions ». Journal of Cheminformatics 1, no 1 (10 juin 2009): 8.
    https://doi.org/10.1186/1758-2946-1-8.
    Returning the opposite of the value so that the metric can be maximized
    """

    def keys(self):
        return ["SAScore"]

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def evaluate_individual(self, individual):
        if individual is None:
            return None
        else:

            mol_graph = MolFromSmiles(individual.to_aromatic_smiles())
            score = sascorer.calculateScore(mol_graph)

            return score, [score]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score


class NormalizedSAScoreEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of SAScore.
    Ertl, Peter, et Ansgar Schuffenhauer. « Estimation of synthetic accessibility score of drug-like molecules based
    on molecular complexity and fragment contributions ». Journal of Cheminformatics 1, no 1 (10 juin 2009): 8.
    https://doi.org/10.1186/1758-2946-1-8.
    Returning the normalized [0, 1] SAScore
    """

    def __init__(self):
        self.sascore_evaluation = SAScoreEvaluationStrategy()

    def keys(self):
        return ["SAScore"]

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def evaluate_individual(self, individual):
        if individual is None:
            return None, [None]
        else:
            unnormalized_sascore, _ = self.sascore_evaluation.evaluate_individual(individual)
            score = 1 - (unnormalized_sascore - 1) / 9

            return score, [score]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score


class QEDEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of population with QED score using RDKit implementation.
    (Bickerton, G. Richard, Gaia V. Paolini, Jérémy Besnard, Sorel Muresan, et Andrew L. Hopkins. « Quantifying the
    Chemical Beauty of Drugs ». Nature Chemistry 4, nᵒ 2 (février 2012): 90‑98. https://doi.org/10.1038/nchem.1243.
    """

    def __init__(self):
        self.scores = None

    def keys(self):
        return ["qed"]

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def evaluate_individual(self, individual):
        if individual is None:
            return None
        else:
            mol_graph = MolFromSmiles(individual.to_aromatic_smiles())
            score = qed(mol_graph)
            return score, [score]

    def record_score(self, idx, new_total_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = new_total_score


class LinearCombinationEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation of the population with a linear combination of given evaluation strategies.
    The coefficients are given in a list of same size as the number of strategies.
    """

    def __init__(self, evaluation_strategies, coefs):

        self.evaluation_strategies = evaluation_strategies
        self.coefs = np.array(coefs).reshape(-1, 1)
        self.strat_keys = []
        for strat in self.evaluation_strategies:
            self.strat_keys.append(strat.keys()[0])

    def get_population_scores(self):

        scores = None

        # Creating lists of scores all the scores of population for all evaluation strategies
        for i, strategy in enumerate(self.evaluation_strategies):
            strategy_evaluation, _ = strategy.get_population_scores()

            if i == 0:
                scores = np.full((len(self.evaluation_strategies), len(strategy_evaluation)), np.nan)

            scores[i] = strategy_evaluation

        total_scores = np.sum(scores * self.coefs, axis=0)

        return total_scores, scores

    def evaluate_individual(self, individual):
        scores = []
        for strategy in self.evaluation_strategies:

            strat_score, _ = strategy.evaluate_individual(individual)
            scores.append(strat_score)

        scores = np.array(scores).reshape(-1, 1)

        return np.sum(scores * self.coefs), scores

    def compute_record_scores(self, population):
        for strategy in self.evaluation_strategies:
            strategy.compute_record_scores(population)

    def record_score(self, idx, new_total_score, new_scores):
        for i, strategy in enumerate(self.evaluation_strategies):
            strategy.record_score(idx, new_scores[i], None)

    def keys(self):
        return self.strat_keys


class ProductSigmLinEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation strategy returning the product of multiple scores after passing them through a linear function and
    a sigmoid function.
    Each score for an individual x is computed as sigm(lin(x)), with specified coefficient.
    """

    def __init__(self, evaluation_strategies, a, b, l):
        """
        Initialization of evaluation strategies and coefficients. All parameters must be lists of the same size.
        :param evaluation_strategies: list of evaluation strategies.
        :param a: list of a coefficients for the linear functions applied to each score in the form ax+b
        :param b: list of b coefficients for the linear functions applied to each score in the form ax+b
        :param l: list of lambda coefficient for the sigmoid functions applied to each score
        """

        self.evaluation_strategies = evaluation_strategies

        # Recording parameters
        self.a = a
        self.b = b
        self.l = l

        # Population initialization
        self.scores = None

        # Recording keys
        self.strat_keys = []
        for strat in self.evaluation_strategies:
            self.strat_keys.append(strat.keys()[0])

    def compute_record_scores(self, population):
        for strategy in self.evaluation_strategies:
            strategy.compute_record_scores(population)

    def get_population_scores(self):

        scores = None

        # Creating lists of scores all the scores of population for all evaluation strategies
        for i, strategy in enumerate(self.evaluation_strategies):
            strategy_evaluation, _ = strategy.get_population_scores()

            if i == 0:
                scores = np.full((len(self.evaluation_strategies), len(strategy_evaluation)), np.nan)

            scores[i] = strategy_evaluation

        total_scores = []
        for curr_ind_scores in scores.T:
            total_scores.append(self._compute_total_score(curr_ind_scores))

        return total_scores, scores

    def _compute_total_score(self, strat_scores):

        tmp_scores = []

        for i, curr_strat_score in enumerate(strat_scores):
            tmp_scores.append(1 / (1 + exp(self.l[i] * (self.a[i] * curr_strat_score + self.b[i]))))

        print("TMP scores : " + str(tmp_scores))

        return np.prod(tmp_scores)

    def evaluate_individual(self, individual):

        strat_scores = []
        for i, eval_strat in enumerate(self.evaluation_strategies):
            # Scoring the individual with the current strategy
            strat_score, _ = eval_strat.evaluate_individual(individual)
            strat_scores.append(strat_score)

        # Computing total score
        total_score = self._compute_total_score(strat_scores)

        # Returning the product of scores
        return total_score, strat_scores

    def record_score(self, idx, new_total_score, new_scores):
        for i, strat in enumerate(self.evaluation_strategies):
            self.evaluation_strategies[i].record_score(idx, new_scores[i], None)

    def keys(self):
        return self.strat_keys


def scores_to_scores_dict(total_scores, scores, keys):
    # Creation of dictionary containing the scores for each evaluator
    step_scores_dict = {}
    for i, k in enumerate(keys):
        step_scores_dict[k] = scores[i]
    step_scores_dict["total"] = total_scores

    return step_scores_dict
