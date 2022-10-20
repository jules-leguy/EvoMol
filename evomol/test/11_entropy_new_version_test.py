from evomol import run_model
from evomol.plot_exploration import exploration_graph

coefs = [1, 1000]

path = "tests/11_entropy_new_version"

run_model(
    {
        "obj_function": {
        "type": "linear_combination",
        "coef": coefs,
        "functions": ["qed",
                      "entropy_checkmol"]},
        "optimization_parameters": {
            "pop_max_size": 1000,
            "k_to_replace": 10,
            "max_steps": 50,
            "mutation_max_depth": 3
        },
        "io_parameters": {
            "model_path": path,
            "record_history": True
        }
    })

exploration_graph(model_path=path, layout="neato", prop_to_study_key="qed")