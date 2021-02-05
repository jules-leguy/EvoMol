import time
from abc import ABC, abstractmethod
from os.path import join, exists

import numpy as np
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.rdmolfiles import MolFromSmiles


class StopCriterionStrategy(ABC):
    """
    Abstract base class of the class that define the stop criterion for the population-based algorithm
    """

    def __init__(self):
        self.pop_alg = None

    def set_pop_alg_instance(self, pop_alg):
        self.pop_alg = pop_alg

    @abstractmethod
    def time_to_stop(self, output_folder_path):
        pass

    def write_stop(self, filepath, message):
        with open(filepath, "w") as f:
            f.write(message)


class MultipleStopCriterionsStrategy(StopCriterionStrategy):

    def __init__(self, stop_criterions):
        super().__init__()
        self.stop_criterions = []
        for stop_criterion in stop_criterions:
            self.stop_criterions.append(stop_criterion)

        self.additional_strategy = []

    def set_additional_strategy(self, additional_stop_strategy):
        self.additional_strategy = [additional_stop_strategy]

    def set_pop_alg_instance(self, pop_alg):
        super(MultipleStopCriterionsStrategy, self).set_pop_alg_instance(pop_alg)
        for stop_criterion in self.stop_criterions + self.additional_strategy:
            stop_criterion.set_pop_alg_instance(pop_alg)

    def time_to_stop(self, output_folder_path):

        for stop_criterion in self.stop_criterions + self.additional_strategy:
            if stop_criterion.time_to_stop(output_folder_path):
                return True

        return False


class KthScoreMaxValue(StopCriterionStrategy):
    """
    Stopping the algorithm if the kth (best) score has reached maximum value of the score
    """

    def __init__(self, max_value, round=None):
        super().__init__()
        self.max_value = max_value
        self.round = round

    def time_to_stop(self, output_folder_path):
        # Extracting the current kth score
        current_kth_score = self.pop_alg.kth_score_history[0]

        if not np.isnan(current_kth_score):

            n_max_scores = 0
            for score in self.pop_alg.curr_total_scores:

                if self.round is None and score >= self.max_value:
                    n_max_scores += 1
                elif self.round is not None and round(score, self.round) >= self.max_value:
                    n_max_scores += 1

            print("n max scores : " + str(n_max_scores))

            if self.round is None:
                test = current_kth_score >= self.max_value
            else:
                test = round(current_kth_score, self.round) >= self.max_value

            if test and output_folder_path:
                self.write_stop(join(output_folder_path, "stop.txt"), "Kth score has max value")

            return test


class KStepsStopCriterionStrategy(StopCriterionStrategy):
    """
    Stopping the algorithm if a given K number of steps is reached
    """

    def __init__(self, n_steps):
        super().__init__()
        self.n_steps = n_steps

    def time_to_stop(self, output_folder_path):
        test = self.n_steps <= self.pop_alg.curr_step_id

        if test and output_folder_path:
            self.write_stop(join(output_folder_path, "stop.txt"), "Max number of steps reached")

        return test


class TimeStopCriterionStrategy(StopCriterionStrategy):
    """
    Stopping the algorithm after a given time
    """

    def __init__(self, max_duration):
        super().__init__()
        self.max_duration = max_duration
        self.start_time = time.time()

    def time_to_stop(self, output_folder_path):

        test = time.time() - self.start_time > self.max_duration

        if test and output_folder_path:
            self.write_stop(join(output_folder_path, "stop.txt"), "Max. time reached (" + str(self.max_duration) + " s)")

        return test


class RadicalFoundStopCriterionStrategy(StopCriterionStrategy):

    def __init__(self):
        super().__init__()

    def _is_radical(self, mol):
        return NumRadicalElectrons(MolFromSmiles(mol.to_aromatic_smiles())) != 0

    def time_to_stop(self, output_folder_path):

        test = False

        for ind in self.pop_alg.pop:
            if ind is not None and self._is_radical(ind):
                test = True
                pass

        if test and output_folder_path:
            self.write_stop(join(output_folder_path, "stop.txt"), "Radical found")

        return test

class FileStopCriterion(StopCriterionStrategy):
    """
    Stopping the algorithm if the given file exists
    """

    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath

    def time_to_stop(self, output_folder_path):

        test = exists(self.filepath)

        if test and output_folder_path:
            self.write_stop(join(output_folder_path, "stop.txt"),
                            "User stop")

        return test
