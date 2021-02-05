from os.path import join

from .evaluation import EvaluationStrategy, GenericFunctionEvaluationStrategy, QEDEvaluationStrategy, \
    NormalizedSAScoreEvaluationStrategy, CLScoreEvaluationStrategy, SAScoreEvaluationStrategy, \
    PenalizedLogPEvaluationStrategy, ZincNormalizedPLogPEvaluationStrategy, LinearCombinationEvaluationStrategy, \
    ProductSigmLinEvaluationStrategy, ProductEvaluationStrategy, SigmLinWrapperEvaluationStrategy
from .evaluation_dft import OPTEvaluationStrategy
from .evaluation_entropy import EntropyContribEvaluationStrategy
from .molgraphops.default_actionspaces import generic_action_space
from .mutation import KRandomGraphOpsImprovingMutationStrategy
from .popalg import PopAlg
from .stopcriterion import MultipleStopCriterionsStrategy, FileStopCriterion, KStepsStopCriterionStrategy
from guacamol.assess_goal_directed_generation import assess_goal_directed_generation
from .guacamol_binding import ChemPopAlgGoalDirectedGenerator, is_or_contains_undefined_GuacaMol_evaluation_strategy, \
    GuacamolEvaluationStrategy, UndefinedGuacaMolEvaluationStrategy, get_GuacaMol_benchmark_parameter


def _is_describing_multi_objective_function(param_eval):
    """
    Determining whether the given parameter describes a multi-objective function
    :param param_eval:
    :return:
    """
    return isinstance(param_eval, dict)


def _is_describing_custom_function(param_eval):
    """
    Determining whether the given parameter describes a custom evaluation function
    :param param_eval:
    :return:
    """
    return callable(param_eval) or (isinstance(param_eval, tuple) and callable(param_eval[0]))


def _is_describing_implemented_function(param_eval):
    """
    Determines whether the given parameter describes an already implemented evaluation function
    :param param_eval:
    :return:
    """

    return param_eval in ["qed", "sascore", "norm_sascore", "plogp", "norm_plogp", "clscore", "homo", "lumo",
                          "entropy_gen_scaffolds", "entropy_ifg", "entropy_shg_1", "entropy_checkmol"] \
           or param_eval.startswith("guacamol")


def _build_evaluation_strategy_from_custom_function(obj_fun_param):
    """
    Building a proper EvaluationStrategy from a custom evaluation function
    :param obj_fun_param:
    :return:
    """

    if isinstance(obj_fun_param, tuple):
        strat = GenericFunctionEvaluationStrategy(evaluation_function=obj_fun_param[0], function_name=obj_fun_param[1])
    elif callable(obj_fun_param):
        strat = GenericFunctionEvaluationStrategy(evaluation_function=obj_fun_param)

    return strat


def _build_evaluation_strategy_from_implemented_function(param_eval, explicit_IO_parameters_dict,
                                                         explicit_search_parameters_dict):
    """
    Building a proper EvaluationStrategy from a string description of an implemented function
    :param param_eval:
    :param kwargs:
    :return:
    """

    if param_eval == "qed":
        strat = QEDEvaluationStrategy()
    elif param_eval == "sascore":
        strat = SAScoreEvaluationStrategy()
    elif param_eval == "norm_sascore":
        strat = NormalizedSAScoreEvaluationStrategy()
    elif param_eval == "plogp":
        strat = PenalizedLogPEvaluationStrategy()
    elif param_eval == "norm_plogp":
        strat = ZincNormalizedPLogPEvaluationStrategy()
    elif param_eval == "clscore":
        strat = CLScoreEvaluationStrategy()
    elif param_eval == "homo":
        strat = OPTEvaluationStrategy("homo",
                                      working_dir_path=explicit_IO_parameters_dict["dft_working_dir"],
                                      cache_files=explicit_IO_parameters_dict["dft_cache_files"])
    elif param_eval == "lumo":
        strat = OPTEvaluationStrategy("lumo",
                                      working_dir_path=explicit_IO_parameters_dict["dft_working_dir"],
                                      cache_files=explicit_IO_parameters_dict["dft_cache_files"])
    elif param_eval == "entropy_ifg":
        strat = EntropyContribEvaluationStrategy(explicit_search_parameters_dict["n_max_desc"],
                                                 pop_size_max=explicit_search_parameters_dict["pop_max_size"],
                                                 descriptor_key="ifg")
    elif param_eval == "entropy_gen_scaffolds":
        strat = EntropyContribEvaluationStrategy(explicit_search_parameters_dict["n_max_desc"],
                                                 pop_size_max=explicit_search_parameters_dict["pop_max_size"],
                                                 descriptor_key="gen_scaffolds")
    elif param_eval == "entropy_shg_1":
        strat = EntropyContribEvaluationStrategy(explicit_search_parameters_dict["n_max_desc"],
                                                 pop_size_max=explicit_search_parameters_dict["pop_max_size"],
                                                 descriptor_key="shg_1")
    elif param_eval == "entropy_checkmol":
        strat = EntropyContribEvaluationStrategy(explicit_search_parameters_dict["n_max_desc"],
                                                 pop_size_max=explicit_search_parameters_dict["pop_max_size"],
                                                 descriptor_key="checkmol")

    # Parameter is a GuacaMol evaluation so the evaluation strategy is undefined for now
    elif param_eval.startswith("guacamol"):
        strat = UndefinedGuacaMolEvaluationStrategy(name=param_eval)


    return strat


def _build_evaluation_strategy_from_single_objective(param_eval, explicit_IO_parameters_dict,
                                                     explicit_search_parameters_dict):
    """
    Building a proper EvaluationStrategy from an evaluation parameter describing a single objective function
    :param param_eval:
    :return:
    """

    # Parameter is an already defined EvaluationStrategy
    if isinstance(param_eval, EvaluationStrategy):
        return param_eval

    # Parameter describes a custom function
    if _is_describing_custom_function(param_eval):
        return _build_evaluation_strategy_from_custom_function(param_eval)

    # Parameter describes an implemented function
    elif _is_describing_implemented_function(param_eval):
        return _build_evaluation_strategy_from_implemented_function(param_eval, explicit_IO_parameters_dict,
                                                                    explicit_search_parameters_dict)


def _build_evaluation_strategy_from_multi_objective(param_eval, explicit_IO_parameters_dict,
                                                    explicit_search_parameters_dict):
    """
    Building a proper EvaluationStrategy from an evaluation parameter describing a multi-objective function
    :param param_eval:
    :param obj_fun_kwargs:
    :return:
    """

    if param_eval["type"] in ["linear_combination", "product", "sigm_lin", "product_sigm_lin"]:

        # Building evaluation strategies
        functions_desc = param_eval["functions"]
        evaluation_strategies = []
        for function_desc in functions_desc:

            if _is_describing_multi_objective_function(function_desc):
                evaluation_strategies.append(_build_evaluation_strategy_from_multi_objective(function_desc,
                                                                                             explicit_IO_parameters_dict,
                                                                                             explicit_search_parameters_dict))
            else:
                evaluation_strategies.append(_build_evaluation_strategy_from_single_objective(function_desc,
                                                                                              explicit_IO_parameters_dict,
                                                                                              explicit_search_parameters_dict))

        if param_eval["type"] == "linear_combination":
            return LinearCombinationEvaluationStrategy(evaluation_strategies, coefs=param_eval["coef"])
        elif param_eval["type"] == "product":
            return ProductEvaluationStrategy(evaluation_strategies)
        elif param_eval["type"] == "sigm_lin":
            return SigmLinWrapperEvaluationStrategy(evaluation_strategies, a=param_eval["a"], b=param_eval["b"],
                                                    l=param_eval["lambda"])
        elif param_eval["type"] == "product_sigm_lin":
            return ProductSigmLinEvaluationStrategy(evaluation_strategies, a=param_eval["a"], b=param_eval["b"],
                                                    l=param_eval["lambda"])


def _parse_objective_function_strategy(parameters_dict, explicit_IO_parameters_dict, explicit_search_parameters_dict):
    """
    Parsing objective function data in the parameters of the model
    Returning an EvaluationStrategy instance.
    :param parameters_dict:
    :return:
    """

    # Extracting objective function parameters
    param_eval = parameters_dict["obj_function"]

    # Parameter is an already defined EvaluationStrategy
    if isinstance(param_eval, EvaluationStrategy):
        eval_strat = param_eval

    # Parameter is a multi-objective function
    elif _is_describing_multi_objective_function(param_eval):
        eval_strat = _build_evaluation_strategy_from_multi_objective(param_eval, explicit_IO_parameters_dict,
                                                                     explicit_search_parameters_dict)

    # Parameter is a single objective function
    else:
        eval_strat = _build_evaluation_strategy_from_single_objective(param_eval, explicit_IO_parameters_dict,
                                                                      explicit_search_parameters_dict)

    return eval_strat


def _parse_action_space(parameters_dict):
    """
    Parsing action space parameters.
    Returning a tuple (ActionSpace instance, ActionSpace.ActionSpaceParameters instance)
    :param parameters_dict:
    :return:
    """

    input_param_action_space = parameters_dict[
        "action_space_parameters"] if "action_space_parameters" in parameters_dict else {}

    explicit_action_space_parameters = {
        "atoms": input_param_action_space["atoms"] if "atoms" in input_param_action_space else "C,N,O,F,P,S,Cl,Br",
        "max_heavy_atoms": input_param_action_space[
            "max_heavy_atoms"] if "max_heavy_atoms" in input_param_action_space else 38,
        "substitution": input_param_action_space[
            "substitution"] if "substitution" in input_param_action_space else True,
        "cut_insert": input_param_action_space["cut_insert"] if "cut_insert" in input_param_action_space else True,
        "move_group": input_param_action_space["move_group"] if "move_group" in input_param_action_space else True,
        "use_rd_filters": input_param_action_space["use_rd_filters"] if "use_rd_filters" in input_param_action_space else False}

    symbols_list = explicit_action_space_parameters["atoms"].split(",")

    action_spaces, action_spaces_parameters = \
        generic_action_space(atom_symbols_list=symbols_list,
                             max_heavy_atoms=explicit_action_space_parameters["max_heavy_atoms"],
                             substitution=explicit_action_space_parameters["substitution"],
                             cut_insert=explicit_action_space_parameters["cut_insert"],
                             move_group=explicit_action_space_parameters["move_group"])

    for parameter in input_param_action_space:
        if parameter not in explicit_action_space_parameters:
            raise RuntimeWarning("Unrecognized parameter in 'action_space_parameters' : " + str(parameter))

    return action_spaces, action_spaces_parameters, explicit_action_space_parameters


def _parse_mutation_parameters(explicit_search_parameters, evaluation_strategy, action_spaces,
                               action_spaces_parameters, search_space_parameters):
    """
    Parsing mutation parameters
    :param parameters_dict:
    :return: Returning a MutationStrategy instance
    """

    mutation_strategy = KRandomGraphOpsImprovingMutationStrategy(k=explicit_search_parameters["mutation_max_depth"],
                                                                 max_n_try=explicit_search_parameters[
                                                                     "mutation_find_improver_tries"],
                                                                 evaluation_strategy=evaluation_strategy,
                                                                 action_spaces=action_spaces,
                                                                 action_spaces_parameters=action_spaces_parameters,
                                                                 problem_type=explicit_search_parameters[
                                                                     "problem_type"],
                                                                 quality_filter=search_space_parameters["use_rd_filters"])

    return mutation_strategy


def _extract_explicit_search_parameters(parameters_dict):
    """
    Building a complete dictionary of search parameters
    :param parameters_dict:
    :return: dictionary of search parameters
    """

    input_search_parameters = parameters_dict[
        "optimization_parameters"] if "optimization_parameters" in parameters_dict else {}

    explicit_search_parameters = {
        "problem_type": input_search_parameters["problem_type"] if "problem_type" in input_search_parameters else "max",
        "pop_max_size": input_search_parameters["pop_max_size"] if "pop_max_size" in input_search_parameters else 1000,
        "k_to_replace": input_search_parameters["k_to_replace"] if "k_to_replace" in input_search_parameters else 10,
        "selection": input_search_parameters["selection"] if "selection" in input_search_parameters else "best",
        "max_steps": input_search_parameters["max_steps"] if "max_steps" in input_search_parameters else 1500,
        "mutable_init_pop": input_search_parameters[
            "mutable_init_pop"] if "mutable_init_pop" in input_search_parameters else True,
        "guacamol_init_top_100": input_search_parameters[
            "guacamol_init_top_100"] if "guacamol_init_top_100" in input_search_parameters else True,
        "mutation_max_depth": input_search_parameters[
            "mutation_max_depth"] if "mutation_max_depth" in input_search_parameters else 2,
        "mutation_find_improver_tries": input_search_parameters[
            "mutation_find_improver_tries"] if "mutation_find_improver_tries" in input_search_parameters else 50,
        "n_max_desc": input_search_parameters["n_max_desc"] if "n_max_desc" in input_search_parameters else 3000000,
        "shuffle_init_pop": input_search_parameters["shuffle_init_pop"] if "shuffle_init_pop" in input_search_parameters else False}

    for parameter in input_search_parameters:
        if parameter not in explicit_search_parameters:
            raise RuntimeWarning("Unrecognized parameter in 'search_parameters' : " + str(parameter))

    return explicit_search_parameters


def _extract_explicit_IO_parameters(parameters_dict):
    """
    Building a complete dictionary of IO parameters
    :param parameters_dict:
    :return: dictionary of search parameters
    """

    input_IO_parameters = parameters_dict["io_parameters"] if "io_parameters" in parameters_dict else {}
    explicit_IO_parameters = {
        "model_path": input_IO_parameters["model_path"] if "model_path" in input_IO_parameters else "EvoMol_model/",
        "smiles_list_init_path": input_IO_parameters[
            "smiles_list_init_path"] if "smiles_list_init_path" in input_IO_parameters else None,
        "smiles_list_init": input_IO_parameters["smiles_list_init"] if "smiles_list_init" in input_IO_parameters else None,
        "external_tabu_list": input_IO_parameters["external_tabu_list"] if "external_tabu_list" in input_IO_parameters else None,
        "save_n_steps": input_IO_parameters["save_n_steps"] if "save_n_steps" in input_IO_parameters else 100,
        "print_n_steps": input_IO_parameters["print_n_steps"] if "print_n_steps" in input_IO_parameters else 1,
        "dft_working_dir": input_IO_parameters[
            "dft_working_dir"] if "dft_working_dir" in input_IO_parameters else "/tmp/",
        "dft_cache_files": input_IO_parameters[
            "dft_cache_files"] if "dft_cache_files" in input_IO_parameters else [],
        "record_history": input_IO_parameters["record_history"] if "record_history" in input_IO_parameters else False,
        "record_all_generated_individuals": input_IO_parameters["record_all_generated_individuals"] if "record_all_generated_individuals" in input_IO_parameters else False}

    for parameter in input_IO_parameters:
        if parameter not in explicit_IO_parameters:
            raise RuntimeWarning("Unrecognized parameter in 'io_parameters' : " + str(parameter))

    return explicit_IO_parameters


def _parse_stop_criterion_strategy(explicit_search_parameters_dict, explicit_IO_parameters_dict):
    """
    Building the StopCriterionStrategy instance
    :param explicit_search_parameters_dict:
    :param explicit_IO_parameters_dict:
    :return:
    """

    stop_criterion_strategy = MultipleStopCriterionsStrategy(
        [KStepsStopCriterionStrategy(explicit_search_parameters_dict["max_steps"]),
         FileStopCriterion(join(explicit_IO_parameters_dict["model_path"], "stop_execution"))])

    return stop_criterion_strategy


def _read_smiles_list_from_file(smiles_list_path):
    """
    Reading a list of SMILES from a file
    :param smiles_list_path: smiles list filepath
    :return: list of SMILES
    """

    with open(smiles_list_path, "r", newline='') as f:
        smiles_list = f.readlines()
    return smiles_list


def _build_instance(evaluation_strategy, mutation_strategy, stop_criterion_strategy, explicit_search_parameters_dict,
                    explicit_IO_parameters_dict):
    """
    Building PopAlg instance
    :param evaluation_strategy: EvaluationStrategy instance
    :param mutation_strategy: MutationStrategy instance
    :param stop_criterion_strategy: MultipleStopCriterionsStrategy instance
    :param explicit_search_parameters_dict: dictionary of search parameters
    :param explicit_IO_parameters_dict: dictionary of IO parameters
    :return:
    """

    pop_alg = PopAlg(

        evaluation_strategy=evaluation_strategy,
        mutation_strategy=mutation_strategy,
        stop_criterion_strategy=stop_criterion_strategy,
        pop_max_size=explicit_search_parameters_dict["pop_max_size"],
        selection=explicit_search_parameters_dict["selection"],
        output_folder_path=explicit_IO_parameters_dict["model_path"],
        save_n_steps=explicit_IO_parameters_dict["save_n_steps"],
        print_n_steps=explicit_IO_parameters_dict["print_n_steps"],
        k_to_replace=explicit_search_parameters_dict["k_to_replace"],
        problem_type=explicit_search_parameters_dict["problem_type"],
        record_history=explicit_IO_parameters_dict["record_history"],
        record_all_generated_individuals=explicit_IO_parameters_dict["record_all_generated_individuals"],
        shuffle_init_pop=explicit_search_parameters_dict["shuffle_init_pop"],
        external_tabu_list=explicit_IO_parameters_dict["external_tabu_list"]
    )

    # Setting the instance for the stop criterion
    pop_alg.stop_criterion_strategy.pop_alg = pop_alg

    # Loading initial population if it is not the GuacaMol special case
    if not is_or_contains_undefined_GuacaMol_evaluation_strategy(evaluation_strategy):

        # Initialization of the PopAlg instance
        pop_alg.initialize()

        # Reading smiles list from "smiles_list" attribute
        if explicit_IO_parameters_dict["smiles_list_init"] is not None:
            pop_alg.load_pop_from_smiles_list(explicit_IO_parameters_dict["smiles_list_init"])

        # Reading smiles list from the file pointed by "smiles_list_init_path" attribute
        else:

            # Initialization of the population with a single methane molecule
            if explicit_IO_parameters_dict["smiles_list_init_path"] is None:
                pop_alg.load_pop_from_smiles_list(smiles_list=["C"],
                                                  atom_mutability=explicit_search_parameters_dict["mutable_init_pop"])
            else:
                pop_alg.load_pop_from_smiles_list(
                    smiles_list=_read_smiles_list_from_file(explicit_IO_parameters_dict["smiles_list_init_path"]),
                    atom_mutability=explicit_search_parameters_dict["mutable_init_pop"])

    return pop_alg


def run_model(parameters_dict):
    # Extracting search parameters
    explicit_search_parameters_dict = _extract_explicit_search_parameters(parameters_dict)

    # Extracting IO parameters
    explicit_IO_parameters_dict = _extract_explicit_IO_parameters(parameters_dict)

    # Building objective function
    evaluation_strategy = _parse_objective_function_strategy(parameters_dict,
                                                             explicit_IO_parameters_dict=explicit_IO_parameters_dict,
                                                             explicit_search_parameters_dict=explicit_search_parameters_dict)

    # Building action space
    action_spaces, action_spaces_parameters, explicit_action_space_parameters = _parse_action_space(parameters_dict)

    # Building mutation strategy
    mutation_strategy = _parse_mutation_parameters(explicit_search_parameters=explicit_search_parameters_dict,
                                                   evaluation_strategy=evaluation_strategy,
                                                   action_spaces=action_spaces,
                                                   action_spaces_parameters=action_spaces_parameters,
                                                   search_space_parameters=explicit_action_space_parameters)

    # Building stop criterion strategy
    stop_criterion_strategy = _parse_stop_criterion_strategy(
        explicit_search_parameters_dict=explicit_search_parameters_dict,
        explicit_IO_parameters_dict=explicit_IO_parameters_dict)

    # Building instance
    pop_alg = _build_instance(evaluation_strategy=evaluation_strategy,
                              mutation_strategy=mutation_strategy,
                              stop_criterion_strategy=stop_criterion_strategy,
                              explicit_search_parameters_dict=explicit_search_parameters_dict,
                              explicit_IO_parameters_dict=explicit_IO_parameters_dict)

    # GuacaMol special case
    if is_or_contains_undefined_GuacaMol_evaluation_strategy(evaluation_strategy):

        model_generator = ChemPopAlgGoalDirectedGenerator(pop_alg=pop_alg,
                                                          guacamol_init_top_100=explicit_search_parameters_dict[
                                                              "guacamol_init_top_100"],
                                                          init_pop_path=explicit_IO_parameters_dict[
                                                              "smiles_list_init_path"],
                                                          output_save_path=explicit_IO_parameters_dict["model_path"])

        # Extracting proper set of benchmarks
        benchmark_parameter = get_GuacaMol_benchmark_parameter(evaluation_strategy)
        benchmark_key = benchmark_parameter.split("_")[1]

        assess_goal_directed_generation(model_generator,
                                        json_output_file=join(explicit_IO_parameters_dict["model_path"],
                                                              "output_GuacaMol.json"),
                                        benchmark_version=benchmark_key)

    else:
        pop_alg.run()
        return pop_alg
