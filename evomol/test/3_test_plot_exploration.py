from evomol import run_model
from evomol.plot_exploration import exploration_graph


# Plotting large exploration tree with colors
model_path = "3_plot_exploration_tree/"

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": 150},
    "io_parameters": {
        "model_path": model_path,
        "record_history": True
    }
})

exploration_graph(model_path=model_path, layout="neato")


# Plotting small exploration tree with images and actions

model_path = "3_plot_exploration_tree_images/"

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": 10,
        "pop_max_size": 10,
        "k_to_replace": 2,
        "mutable_init_pop": False
    },
    "io_parameters": {
        "model_path": model_path,
        "record_history": True,
        "smiles_list_init_path": "acetylsalicylic_acid.smi"
    }
})


exploration_graph(model_path=model_path, layout="dot", draw_actions=True, plot_images=True, draw_scores=True,
                  root_node="O=C(C)Oc1ccccc1C(=O)O", legend_scores_keys_strat=["total"], mol_size=0.3,
                  legend_offset=(-0.007, -0.05), figsize=(20, 20/1.5), legends_font_size=13)