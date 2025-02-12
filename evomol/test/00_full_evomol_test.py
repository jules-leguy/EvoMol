from rdkit.Chem import Lipinski
from rdkit.Chem.rdmolfiles import MolFromSmiles
from evomol import run_model

silly_molecules_reference_db_path = "/home/jleguy/Documents/postdoc/data/EvoMol_data/ECF4_filters/complete_ChEMBL_ecfp4_dict.json"


### 1 : test objective functions

def count_nitrogen(smi):
    n_nitrogen \
        = 0
    for c in smi:
        if c == "N" or c == "n":
            n_nitrogen += 1
    return n_nitrogen


def hetero_atoms_proportion(smiles):
    return Lipinski.NumHeteroatoms(MolFromSmiles(smiles)) / Lipinski.HeavyAtomCount(MolFromSmiles(smiles))


eval_fun = ["qed", "plogp", "norm_plogp", "sascore", "norm_sascore", "clscore", (count_nitrogen, "count_N"),
            "rediscovery_CC(=O) OC1=CC=CC=C1C(=O)O",
            {"type": "linear_combination",
             "functions": [(count_nitrogen, "count_N"), "qed", "clscore"],
             "coef": [0.5, 0.5, 0]},
            {"type": "linear_combination",
             "functions": ["isomer_C8H9NO2", "qed", "clscore"],
             "coef": [1, 1, 1]},
            {
                "type": "product",
                "coef": [1, 1],
                "functions": ["qed",
                              {
                                  "type": "linear_combination",
                                  "coef": [1, 1],
                                  "functions": ["plogp", "sascore"]
                              }]
            },
            {
                "type": "linear_combination",
                "functions": [
                    "qed",
                    "sillywalks_proportion"
                ],
                "coef": [1, -1]
            },
            {
                "type": "linear_combination",
                "coef": [1, 1000, 1000, 1000, 1000, 1000],
                "functions": ["qed", "entropy_checkmol", "entropy_ifg", "entropy_gen_scaffolds", "entropy_shg_1",
                              "entropy_ecfp4"]
            },
            {
                # The objective function is a linear combination of three sub-functions
                "type": "linear_combination",
                "coef": [0.5, 0.25, 0.25],
                "functions": [
                    {
                        # Centering a Gaussian function on the 0.8 QED value
                        "type": "gaussian",
                        "function": "qed",
                        "mu": 0.8,
                        "sigma": 0.5,
                        "normalize": True
                    },
                    {
                        # The parameters of the sigmoid function are chosen so that the function value is above 0.99 when the CLScore value is above 4
                        "type": "sigm_lin",
                        "function": "clscore",
                        "a": -1,
                        "b": 3.5,
                        "lambda": 10
                    },
                    {
                        # Centering a Gaussian function on the 0.7 proportion of hetero atoms
                        "type": "gaussian",
                        "function": {
                            # Returning 1 - (proportion of heteroatoms) to represent the proportion of carbon atoms
                            "type": "one_minus",
                            "function": (hetero_atoms_proportion, "hetero_atoms_proportion")
                        },
                        "mu": 0.7,
                        "sigma": 0.5,
                        "normalize": True
                    },
                ]
            }
            ]

for i, eval_fun in enumerate(eval_fun):
    run_model({
        "obj_function": eval_fun,
        "optimization_parameters": {
            "pop_max_size": 10,
            "k_to_replace": 2,
            "max_steps": 10 if eval_fun == "homo" or eval_fun == "lumo" or isinstance(eval_fun, dict) else 50,
            "problem_type": "min" if eval_fun == "sascore" or eval_fun == "lumo" else "max"
        },
        "action_space_parameters": {
            "atoms": "C,N,O,F" if eval_fun == "homo" or eval_fun == "lumo" or isinstance(eval_fun,
                                                                                         dict) else "C,N,O,F,P,S,Cl,Br",
            "use_rd_filters": True
        },
        "io_parameters": {
            "model_path": '00_full_test_1_test_evaluation_functions/' + str(i) + "_" + str(eval_fun)[:100],
            "silly_molecules_reference_db_path": silly_molecules_reference_db_path
        }
    })

### 2 : test guacamol

run_model({
    "obj_function": "guacamol_v2",
    "optimization_parameters": {
        "max_steps": 50,
        "pop_max_size": 1000,
        "guacamol_init_top_100": False},
    "io_parameters": {
        "model_path": "00_full_test_2_guacamol/"
    },
})

### 3 : test plot exploration

from evomol.plot_exploration import exploration_graph

#
#
# Plotting large exploration tree with colors
model_path = "00_full_test_3_plot_exploration_tree/"

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": 200},
    "io_parameters": {
        "model_path": model_path,
        "record_history": True
    }
})

exploration_graph(model_path=model_path, layout="neato")

# Plotting small exploration tree with images and actions

model_path = "00_full_test_3_plot_exploration_tree_images/"

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
                  root_node="O=C(C)Oc1ccccc1C(=O)O", legend_scores_keys_strat=["total"], mol_size_inches=0.055,
                  mol_size_px=(200, 200), figsize=(30, 30 / 1.5), legend_offset=(-0.007, -0.05), legends_font_size=13)

### 4 : large entropy calculation test

from random import randint, choice
from evomol import EntropyContribEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph


def gen_rand_ind(smi_size_max):
    smi_size = randint(1, smi_size_max)
    smi = ""
    for i in range(smi_size):
        smi += choice(["C", "O", "N", "P", "S"])

    return smi


def gen_random_pop(pop_size, smi_size_max):
    pop = []

    for i in range(pop_size):
        pop.append(MolGraph(MolFromSmiles(gen_rand_ind(smi_size_max))))

    return pop


def loop_test(pop_size, mol_size_max, n_tests_per_loop, n_loops, max_desc):
    for loop in range(n_loops):

        print("loop : " + str(loop))

        pop = gen_random_pop(pop_size, mol_size_max)
        s = EntropyContribEvaluationStrategy(max_desc, pop_size, descriptor_key="atoms")

        s.compute_record_scores_init_pop(pop)

        for i in range(n_tests_per_loop):

            ind_to_replace_idx = randint(0, pop_size - 1)
            ind_to_add = gen_rand_ind(mol_size_max)

            print("desc count : " + str(s.desc_count))

            print("Replacing " + str(pop[ind_to_replace_idx]) + " by " + str(ind_to_add))

            true_delta = s.true_delta_ent(ind_to_replace_idx, MolGraph(MolFromSmiles(ind_to_add)))
            computed_delta = -s.scores[ind_to_replace_idx] + \
                             s.evaluate_individual(MolGraph(MolFromSmiles(ind_to_add)), ind_to_replace_idx)[0]

            print("loss replace : " + str(s.scores[ind_to_replace_idx]))
            print(
                "gain add : " + str(s.evaluate_individual(MolGraph(MolFromSmiles(ind_to_add)), ind_to_replace_idx)[0]))

            print("true delta : " + str(true_delta))
            print("computed delta : " + str(computed_delta))

            if abs(true_delta - computed_delta) > 0.0001:
                print("Error !")
                exit(1)


loop_test(pop_size=100000, mol_size_max=10, n_tests_per_loop=100, n_loops=10, max_desc=500000)

### 5 : test record all generated

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": 5,
        "pop_max_size": 1,
        "k_to_replace": 1,
        "problem_type": "max",
        "mutation_find_improver_tries": 50,
    },
    "io_parameters": {
        "model_path": "00_full_test_5_test_record_all_generated",
        "smiles_list_init": ["C"],
        "record_history": False,
        "save_n_steps": 1,
        "record_all_generated_individuals": True
    },
    "action_space_parameters": {
        "atoms": "C,N,O,F",
        "max_heavy_atoms": 9,
    }
})

### 6 : test stop at N calls

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": float("inf"),
        "max_obj_calls": 500
    },
    "io_parameters": {
        "model_path": "00_full_test_6_test_n_calls_stop"
    },
})

### 7 : test record failed filters

run_model({
    "obj_function": "qed",
    "action_space_parameters": {
        "sillywalks_threshold": 0.1,
        "sascore_threshold": 4.4
    },
    "optimization_parameters": {
        "max_steps": 200
    },
    "io_parameters": {
        "model_path": "00_full_test_7_test_failed_filters",
        "save_n_steps": 1,
        "silly_molecules_reference_db_path": silly_molecules_reference_db_path,
        "record_all_generated_individuals": True
    }
})
