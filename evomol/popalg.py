import csv
import time
from collections import deque
from os import makedirs
from os.path import dirname, join

import numpy as np
from evomol.evaluation import EvaluationError, scores_to_scores_dict
from evomol.molgraphops.molgraph import MolGraph
from evomol.mutation import NoImproverError, MutationError
from rdkit.Chem.rdmolfiles import MolFromSmiles


class NoMoreIndToMutate(Exception):
    pass


class PopAlg:
    """
    Class running the population algorithm defined by the given strategies.
    """

    def __init__(self, evaluation_strategy, mutation_strategy, stop_criterion_strategy,
                 output_folder_path="EvoMol_model/", pop_max_size=1000, k_to_replace=10, save_n_steps=100,
                 print_n_steps=1, kth_score_to_record=1, record_history=False, problem_type="max"):
        """
        :param evaluation_strategy: EvaluationStrategy instance to evaluate individuals
        :param mutation_strategy: MutationStrategy instance to mutate solutions and find improvers
        :param stop_criterion_strategy: StopCriterionStrategy instance to stop the search when a condition is reached
        :param output_folder_path: Path of folder where the data is recorded (default : "EvoMol_model/")
        :param pop_max_size: Maximum population size (default : 1000)
        :param k_to_replace: Number of individuals to be replaced each step (default : 10)
        :param save_n_steps: Frequency of saving the model (default : 100)
        :param print_n_steps: Frequency of printing the results (default : 1)
        :param kth_score_to_record: Kth score to be recorded for premature stop
        :param record_history: Whether to record history of actions (necessary to draw exploration trees)
        :param problem_type: Whether it is a maximization ("max") or minimization ("min") problem. (default : "max")
        """

        # Loading problem type
        self.problem_type = problem_type

        # Loading strategy modules
        self.evaluation_strategy = evaluation_strategy
        self.mutation_strategy = mutation_strategy
        self.stop_criterion_strategy = stop_criterion_strategy

        # Updating PopAlg instance in StopCriterionStrategy instance
        self.stop_criterion_strategy.set_pop_alg_instance(self)

        # Saving population's max size and the number of individuals to replace each step
        self.pop_max_size = pop_max_size
        self.k_to_replace = k_to_replace

        # Frequency of saving and printing results
        self.save_n_steps = save_n_steps
        self.print_n_steps = print_n_steps

        # Recording output files paths
        self.output_folder_path = output_folder_path

        # Kth score to record
        self.kth_score_to_record = kth_score_to_record

        # History recording parameter
        self.record_history = record_history

        # Attributes initialization
        self.pop = None
        self.pop_tabu_list = None
        self.step_traces = None
        self.curr_step_id = None
        self.errors = None
        self.curr_total_scores = None
        self.curr_scores = None
        self.kth_score_history = None
        self.n_success_mut = None
        self.n_fail_mut = None
        self.actions_history = None
        self.removed_actions_score_smi_tuple = None
        self.timestamp_start = None

    def initialize(self):
        """
        Initialization of EvoMol with starting values.
        This method MUST BE CALLED BEFORE running the algorithm.
        :return:
        """

        # Initialization of population
        self.pop = list(np.full((self.pop_max_size,), None))

        # Initialization of the dictionaries containing the smiles of former and current individuals as keys
        self.pop_tabu_list = list(np.full((self.pop_max_size,), None))

        # Initialization of the dictionary containing the traces of steps of the algorithm
        self.step_traces = {
            'scores': {},
            'n_replaced': [],
            'additional_values': {},
            'timestamps': []
        }

        # Initialization of keys in the self.step_traces dict declared by the evaluation strategy instance
        for k in self.evaluation_strategy.keys() + ["total"]:
            for stat in ["mean", "med", "min", "max", "std"]:
                self.step_traces["scores"][k + "_" + stat] = []

        # Initialization of the step counter.
        self.curr_step_id = 0

        # Initialization of errors list
        self.errors = []
        self.curr_total_scores = None
        self.curr_scores = None
        self.timestamp_start = None

        self.kth_score_history = deque(maxlen=500)

        self.n_success_mut = np.zeros(self.pop_max_size, dtype=np.int)
        self.n_fail_mut = np.zeros(self.pop_max_size, dtype=np.int)

        self.actions_history = list(np.full(self.pop_max_size, None))
        self.removed_actions_score_smi_tuple = {}

        # Computing start timestamp
        self.timestamp_start = time.time()

    def load_pop_from_smiles_list(self, smiles_list, atom_mutability=True):
        """
        Loading the population from the given smiles list.
        Setting the internal variables to their values
        :param smiles_list: list of SMILES
        :param atom_mutability: whether the core of the molecules of the starting population can be modified
        :return:
        """

        # Iterating over all the given smiles
        for i, smi in enumerate(smiles_list):

            # Loading QuMolGraph object
            self.pop[i] = MolGraph(MolFromSmiles(smi), sanitize_mol=True, mutability=atom_mutability)

            # Saving smiles in the tabu dictionary and in action history initialization
            self.pop_tabu_list[i] = self.pop[i].to_aromatic_smiles()
            self.actions_history[i] = self.pop[i].to_aromatic_smiles()

        # Evaluation of the population
        print("Computing descriptors at initialization...")
        self.evaluation_strategy.compute_record_scores(self.pop)

    def save(self):
        """
        Saving the data to the files
        :return:
        """
        if self.output_folder_path:

            # Creating directories if they don't exist
            makedirs(dirname(join(self.output_folder_path, "file")), exist_ok=True)

            # ### Steps data ###
            csv_array = []
            for k, v in self.step_traces["scores"].items():
                csv_array.append([k] + v)
            csv_array.append(["n_replaced"] + self.step_traces["n_replaced"])
            csv_array.append(["timestamps"] + self.step_traces["timestamps"])

            with open(join(self.output_folder_path, 'steps.csv'), "w") as f:
                writer = csv.writer(f)
                for row in np.array(csv_array).T:
                    writer.writerow(row)

            # ### Last step population data ###
            csv_array = []

            # Mutation success history
            n_success_mut_str = []
            n_fail_mut_str = []
            for i, ind in enumerate(self.pop):
                n_success_mut_str.append(str(self.n_success_mut[i]))
                n_fail_mut_str.append(str(self.n_fail_mut[i]))

            csv_array.append(["smiles"] + self.pop_tabu_list)

            # Mutation success and failures
            csv_array.append(["n_success_mut"] + n_success_mut_str)
            csv_array.append(["n_failures_mut"] + n_fail_mut_str)

            # Scores data
            self.curr_total_scores, self.curr_scores = self.evaluation_strategy.get_population_scores()
            step_scores_dict = scores_to_scores_dict(self.curr_total_scores,
                                                     self.curr_scores,
                                                     self.evaluation_strategy.keys())

            for k, scores_list in step_scores_dict.items():
                scores_list_np = np.full((self.pop_max_size,), None)
                scores_list_np[:len(scores_list)] = scores_list
                csv_array.append([k] + list(scores_list_np))

            # Action history data
            csv_array.append(["history_data"] + self.actions_history)

            with open(join(self.output_folder_path, 'pop.csv'), "w") as f:
                writer = csv.writer(f)
                for row in np.array(csv_array).T:
                    writer.writerow(row)

            # ### Removed individuals actions recording ###
            if self.record_history:
                with open(join(self.output_folder_path, 'removed_ind_act_history.csv'), "w") as f:

                    writer = csv.writer(f)
                    writer.writerow(["history_data", "total"] + self.evaluation_strategy.keys() + ["smiles"])

                    for removed_act_history in self.removed_actions_score_smi_tuple.keys():
                        if removed_act_history != "":
                            total_score = self.removed_actions_score_smi_tuple[removed_act_history][0]
                            scores = self.removed_actions_score_smi_tuple[removed_act_history][1]
                            smi = self.removed_actions_score_smi_tuple[removed_act_history][2]

                            writer.writerow([removed_act_history, total_score] + list(scores) + [smi])

            # ### Errors data ###
            with open(join(self.output_folder_path, 'errors.csv'), "w") as f:
                writer = csv.writer(f)
                writer.writerow(["step", "error"])
                for error in self.errors:
                    writer.writerow(error)

    def record_step_data(self):

        # Extracting scores dictionary containing the scores for each objective
        step_scores_dict = scores_to_scores_dict(self.curr_total_scores, self.curr_scores,
                                                 self.evaluation_strategy.keys())

        # Saving statistics on current step scores
        for k, scores_list in step_scores_dict.items():
            step_mean = np.mean(scores_list)
            step_min = np.min(scores_list)
            step_max = np.max(scores_list)
            step_med = np.median(scores_list)
            step_std = np.std(scores_list)

            self.step_traces["scores"][k + "_mean"].append(step_mean)
            self.step_traces["scores"][k + "_med"].append(step_med)
            self.step_traces["scores"][k + "_min"].append(step_min)
            self.step_traces["scores"][k + "_max"].append(step_max)
            self.step_traces["scores"][k + "_std"].append(step_std)

            if self.curr_step_id % self.print_n_steps == 0:
                print(k + "_mean : " + str("%.5f" % step_mean))
                print(k + "_med : " + str("%.5f" % step_med))
                print(k + "_std : " + str("%.5f" % step_std))
                print(k + "_min : " + str("%.5f" % step_min))
                print(k + "_max : " + str("%.5f" % step_max))

    def evaluate_pop_record_step_data(self, n_replaced, record_step_data=True):

        # Population evaluation
        self.curr_total_scores, self.curr_scores = self.evaluation_strategy.get_population_scores()

        # Updating the history of the kth score
        if len(self.curr_total_scores) >= self.kth_score_to_record:
            kth_score = np.partition(self.curr_total_scores, -self.kth_score_to_record)[-self.kth_score_to_record]
            self.kth_score_history.appendleft(kth_score)
        else:
            self.kth_score_history.appendleft(np.nan)

        if record_step_data:
            # Recording step data
            self.record_step_data()

            # Saving step timestamp
            self.step_traces["timestamps"].append(time.time() - self.timestamp_start)

            # Saving the number of replaced individuals
            self.step_traces["n_replaced"].append(n_replaced)

    def select_to_be_replaced(self):

        # Computing number of defined individuals
        n_defined_ind = len(self.curr_total_scores)

        # Computing priority order of undefined individuals
        undefined_priority_order = np.arange(n_defined_ind, self.pop_max_size)

        defined_priority_order = None

        if self.problem_type == "max":
            defined_priority_order = np.argsort(self.curr_total_scores)
        elif self.problem_type == "min":
            defined_priority_order = np.argsort(self.curr_total_scores)[::-1]

        # Computing complete order
        to_be_replaced_indices = list(undefined_priority_order) + list(defined_priority_order)

        return to_be_replaced_indices[:self.k_to_replace]

    def sort_to_be_mutated(self):

        # Computing max of valid individual (that have not reached the max number of mutation failure)
        ind_valid_idx = np.argwhere(self.n_fail_mut[:len(self.curr_total_scores)] < float("inf")).flatten()

        # Extracting the scores of the valid individual
        scores_valid = np.array(self.curr_total_scores)[ind_valid_idx]

        to_be_mutated_in_order_mask = None

        # Sorting in descending order the scores of the valid individuals
        if self.problem_type == "max":
            to_be_mutated_in_order_mask = ind_valid_idx[np.argsort(scores_valid)[::-1].flatten()]
        elif self.problem_type == "min":
            to_be_mutated_in_order_mask = ind_valid_idx[np.argsort(scores_valid).flatten()]

        return to_be_mutated_in_order_mask

    def run(self):
        """
        Running the algorithm
        :return:
        """

        try:

            print("Start pop algorithm")
            if self.curr_step_id == 0:
                # Evaluation of population at initialization and recording data
                self.evaluate_pop_record_step_data(n_replaced=0, record_step_data=True)
            else:
                # Evaluation of population at hot start step
                self.evaluate_pop_record_step_data(n_replaced=None, record_step_data=False)

            # Running the algorithm while the stop criterion is not reached
            while not self.stop_criterion_strategy.time_to_stop(self.output_folder_path):

                print("new step")

                if self.curr_step_id % self.print_n_steps == 0:
                    print("step : " + str(self.curr_step_id))

                # Selecting the individuals to be replaced
                curr_to_be_replaced_indices = self.select_to_be_replaced()

                # Selecting the individuals to be mutated in priority order
                to_be_mutated_indices = self.sort_to_be_mutated()

                n_mutated_tries = 0
                replaced_during_step = []

                if self.problem_type == "max":
                    print("best : " + str(self.pop[np.argmax(self.curr_total_scores)]))
                elif self.problem_type == "min":
                    print("best : " + str(self.pop[np.argmin(self.curr_total_scores)]))

                try:

                    # Iterating over all individuals to be replaced
                    for curr_to_be_replaced_idx in curr_to_be_replaced_indices:

                        replacement_successful = False

                        while not replacement_successful:

                            # Individual to mutate is available
                            if n_mutated_tries < len(to_be_mutated_indices):

                                curr_to_be_mutated_idx = to_be_mutated_indices[n_mutated_tries]

                                try:

                                    # Extracting the score of the individual to be replaced if defined
                                    if curr_to_be_replaced_idx < len(self.curr_total_scores):
                                        curr_to_be_replaced_total_score = self.curr_total_scores[
                                            curr_to_be_replaced_idx]
                                        curr_to_be_replaced_scores = self.curr_scores.T[curr_to_be_replaced_idx]

                                    else:
                                        curr_to_be_replaced_total_score = None
                                        curr_to_be_replaced_scores = None

                                    curr_to_be_replaced_smiles = self.pop_tabu_list[curr_to_be_replaced_idx]

                                    # Trying to perform mutation
                                    mutated_ind, mutation_desc, mutated_total_score, mutated_scores = \
                                        self.mutation_strategy.mutate(
                                            individual=self.pop[curr_to_be_mutated_idx],
                                            ind_to_replace=self.pop[curr_to_be_replaced_idx],
                                            curr_total_score=curr_to_be_replaced_total_score,
                                            pop_tabu_list=self.pop_tabu_list)

                                    # Saving the new individual in the pop smiles list
                                    self.pop_tabu_list[curr_to_be_replaced_idx] = mutated_ind.to_aromatic_smiles()

                                    if self.record_history:

                                        if self.actions_history[curr_to_be_replaced_idx] is not None:
                                            # Recording the history of actions, score and additional values for the
                                            # replaced individual
                                            self.removed_actions_score_smi_tuple[
                                                self.actions_history[curr_to_be_replaced_idx]] = \
                                                (curr_to_be_replaced_total_score,
                                                 curr_to_be_replaced_scores,
                                                 curr_to_be_replaced_smiles)

                                        # Recording the history of actions for new individual
                                        self.actions_history[curr_to_be_replaced_idx] = \
                                            self.actions_history[curr_to_be_mutated_idx] + "|" + mutation_desc

                                    # Replacing individual
                                    self.pop[curr_to_be_replaced_idx] = mutated_ind
                                    self.n_success_mut[curr_to_be_replaced_idx] = 0
                                    self.n_fail_mut[curr_to_be_replaced_idx] = 0

                                    # Updating score
                                    self.evaluation_strategy.record_score(curr_to_be_replaced_idx, mutated_total_score,
                                                                          mutated_scores)

                                    # Recording success
                                    replaced_during_step.append((curr_to_be_replaced_idx, mutated_ind))
                                    replacement_successful = True
                                    self.n_success_mut[curr_to_be_mutated_idx] += 1

                                except NoImproverError as err:
                                    self.n_fail_mut[curr_to_be_mutated_idx] += 1
                                    self.errors.append([self.curr_step_id, str(err)])
                                except MutationError as err:
                                    self.n_fail_mut[curr_to_be_mutated_idx] += 1
                                    self.errors.append([self.curr_step_id, str(err)])
                                except EvaluationError as err:
                                    self.n_fail_mut[curr_to_be_mutated_idx] += 1
                                    self.errors.append([self.curr_step_id, str(err)])

                                finally:
                                    n_mutated_tries += 1

                            # No more individual is available to be mutated
                            else:
                                raise NoMoreIndToMutate()

                except NoMoreIndToMutate:

                    self.errors.append([self.curr_step_id, "No more individual to be mutated"])

                    # Saving information if no individual was replaced during step
                    if len(replaced_during_step) == 0:

                        print("No replacement occurred")
                        self.errors.append([self.curr_step_id, "No replacement occured"])

                # Recording the number of replaced individuals for the step
                n_replaced = len(replaced_during_step)

                # Evaluation of new population and recording step data
                self.evaluate_pop_record_step_data(n_replaced=n_replaced)

                if self.curr_step_id % self.save_n_steps == 0:
                    self.save()

                # Updating curr step id
                self.curr_step_id += 1

            print("Stopping : stop condition reached")

        except KeyboardInterrupt:
            print("Stopping : interrupted by user")
        except MemoryError:
            print("Stopping : no memory available")
        finally:

            # Saving algorithm result
            self.save()
