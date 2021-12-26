from evomol import run_model

model_path = "32_debut_plot_tree/"

# run_model({
#     "obj_function": "qed",
#     "optimization_parameters": {
#         "max_steps": 100},
#     "io_parameters": {
#         "model_path": model_path,
#         "record_history": True
#     }
# })

from evomol.plot_exploration import exploration_graph
exploration_graph(model_path=model_path, layout="neato")
