from evomol import run_model


def count_nitrogen(smi):
    n_nitrogen \
        = 0
    for c in smi:
        if c == "N" or c == "n":
            n_nitrogen += 1
    return n_nitrogen


eval_functions = ["qed", "plogp", "norm_plogp", "sascore", "norm_sascore", "clscore", "homo", "lumo", (count_nitrogen,
                                                                                                       "count_N"),
                  {"type": "linear_combination",
                   "functions": [(count_nitrogen, "count_N"), "qed", "clscore"],
                   "coefs": [0.5, 0.5, 0]},
                  {"type": "product_sigm_lin",
                   "functions": ["homo", "clscore"],
                   "a": [-1, -1],
                   "b": [-7, 1.5],
                   "lambda": [1, 10]},
                  ]

for i, eval_fun in enumerate(eval_functions):
    run_model({
        "obj_function": eval_fun,
        "search_parameters": {
            "pop_max_size": 10,
            "k_to_replace": 2,
            "max_steps": 10 if eval_fun == "homo" or eval_fun == "lumo" or isinstance(eval_fun, dict) else 50,
            "problem_type": "min" if eval_fun == "sascore" or eval_fun == "lumo" else "max"
        },
        "action_space_parameters": {
            "atoms": "C,N,O,F" if eval_fun == "homo" or eval_fun == "lumo" or isinstance(eval_fun,
                                                                                         dict) else "C,N,O,F,P,S,Cl,Br"
        },
        "io_parameters": {
            "model_path": '1_test_evaluation_functions/' + str(i) + "_" + str(eval_fun),
            "dft_working_dir": "/home/jleguy/dft_comput/",
            "dft_cache_files": ["/home/jleguy/Documents/these/prod/data/00_datasets/DFT/cache_OPT.json"]
        }
    })
