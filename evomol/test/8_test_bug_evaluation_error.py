import os

from evomol import run_model

output_model_path = "8_test_bug_evaluation_error"

pop_size = 1000
k_to_replace = 10
max_steps = 700
mutation_max_depth = 3

proportion_qed = 0.001


# Generate coefficients

def gen_coef(prop_qed):
    return [prop_qed, (1 - prop_qed) / 2, (1 - prop_qed) / 2]


coefficients = []
coefs = gen_coef(proportion_qed)

run_model(
    {
        "obj_function": {
            "type": "linear_combination",
            "coef": coefs,
            "functions": ["qed",
                          "entropy_ifg",
                          "entropy_gen_scaffolds"]},
        "optimization_parameters": {
            "pop_max_size": pop_size,
            "k_to_replace": k_to_replace,
            "max_steps": max_steps,
            "mutation_max_depth": mutation_max_depth
        },
        "io_parameters": {
            "model_path": output_model_path,
            "record_history": True
        }
    })
