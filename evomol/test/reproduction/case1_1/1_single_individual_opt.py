from os.path import join

from evomol import run_model

output_results_prefix_path = "results/"

params = [
    {"obj_function": "qed",
     "atoms": "C,N,O,F,P,S,Cl,Br",
     "name": "QED"},

    {"obj_function": "plogp",
     "atoms": "C,N,O,F,P,S,Cl,Br",
     "name": "plogP"},

    {"obj_function": "norm_plogp",
     "atoms": "C,N,O,F,P,S,Cl,Br",
     "name": "norm_plogP"},

    {"obj_function": "plogp",
     "atoms": "C,N,O,F",
     "name": "plogP_CNOF"},

    {"obj_function": "norm_plogp",
     "atoms": "C,N,O,F",
     "name": "norm_plogP_CNOF"},

    {"obj_function": "norm_sascore",
     "atoms": "C,N,O,F,P,S,Cl,Br",
     "name": "norm_sascore"},

    {"obj_function": "clscore",
     "atoms": "C,N,O,F,P,S,Cl,Br",
     "name": "clscore"}
]


def run_model_params(params):
    run_model({
        "obj_function": params["obj_function"],
        "search_parameters": {
            "pop_max_size": 1,
            "k_to_replace": 1,
        },
        "action_space_parameters": {
            "atoms": params["atoms"]
        },
        "io_parameters": {
            "model_path": join(output_results_prefix_path, params["name"]),
            "save_n_steps": 1,
        }
    })


for exp_param in params:
    run_model_params(exp_param)
