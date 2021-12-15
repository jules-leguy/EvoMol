from evomol import run_model

run_model({
    "obj_function": {
        "type": "linear_combination",
        "functions": [
            {
                "type": "opposite",
                "function": "qed"
            },
            "n_perturbations"
        ],
        "coef": [1, -1]
    },
    "optimization_parameters": {
        "max_steps": 100,
        "pop_max_size": 1000,
        "mutation_max_depth": 1,
        "selection": "random_weighted"
    },
    "io_parameters": {
        "model_path": "30_test1/",
        "record_history": True,
        "record_all_generated_individuals": True
    },
})


run_model({
    "obj_function": {
        "type": "linear_combination",
        "functions": [
            {
                "type": "opposite",
                "functions": ["qed"]
            },
            "n_perturbations"
        ],
        "coef": [1, -1]
    },
    "optimization_parameters": {
        "max_steps": 100,
        "pop_max_size": 1000,
        "mutation_max_depth": 1,
        "selection": "random_weighted"
    },
    "io_parameters": {
        "model_path": "30_test2/",
        "record_history": True,
        "record_all_generated_individuals": True
    },
})
