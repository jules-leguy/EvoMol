from evomol import run_model

run_model({
    "obj_function": {
            "type": "linear_combination",
            "coef": [1, 0.1],
        "functions": ["guacamol_v2",
                      "entropy_shg_1"]
    },
    "optimization_parameters": {
        "max_steps": 100,
        "pop_max_size": 1000,
        "guacamol_init_top_100": False},
    "io_parameters": {
        "model_path": "14_calls/"
    },
})