from evomol import run_model

run_model(
    {"obj_function": {
        "type": "linear_combination",
        "coef": [0.5, 0.5],
        "functions": ["entropy_ifg",
                      "entropy_gen_scaffolds"]},
        "optimization_parameters": {
            "pop_max_size": 100,
            "k_to_replace": 10,
            "max_steps": 200
        },
        "io_parameters": {
            "model_path": '5_test_entropy/'
        }
    }
)
