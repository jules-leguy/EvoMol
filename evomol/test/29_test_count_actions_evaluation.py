from evomol import run_model

run_model({
    "obj_function": {
        "type": "linear_combination",
        "functions": [
            "qed",
            "n_perturbations"
        ],
        "coef": [1, -1]
    },
    "optimization_parameters": {
        "max_steps": 200,
        "pop_max_size": 1000,
        "mutation_max_depth": 1,
        "selection": "random_weighted"
    },
    "io_parameters": {
        "model_path": "28_count_actions/",
        "record_history": True,
        "record_all_generated_individuals": True
    },
})
