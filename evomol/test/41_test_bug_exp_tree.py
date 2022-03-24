from evomol.plot_exploration import exploration_graph

path = "41_test_bug_exp_tree/"

exploration_graph(model_path=path, draw_actions=True, plot_images=True, draw_scores=True,
                  root_node="C", legend_scores_keys_strat=["total"], mol_size_inches=0.06, mol_size_px=(200, 200),
                  legend_offset=(-0.01, -0.03), figsize=(15, 15 / 1.5), legends_font_size=9)
