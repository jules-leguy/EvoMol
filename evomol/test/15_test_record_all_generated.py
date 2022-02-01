import os

from evomol import run_model

cache_files = [os.environ["DATA"] + "/00_datasets/DFT/cache_OD9_step7.json",
               os.environ["DATA"] + "/00_datasets/DFT/cache_OPT_OD9_Marta_filtered.json",
               os.environ["DATA"] + "/00_datasets/DFT/cache_OPT_OD9_0.json",
               os.environ["DATA"] + "/00_datasets/DFT/cache_OPT.json"]

run_model({
    "obj_function": "homo",
    "optimization_parameters": {
        "max_steps": 20,
        "pop_max_size": 1,
        "k_to_replace": 1,
        "problem_type": "max",
        "mutation_find_improver_tries": 50,
    },
    "io_parameters": {
        "model_path": "15_test_record_all_generated",
        "smiles_list_init": ["C"],
        "record_history": False,
        "dft_working_dir": "/home/jleguy/dft_comput",
        "dft_cache_files": [],
        "save_n_steps": 1,
        "record_all_generated_individuals": True
    },
    "action_space_parameters": {
        "atoms": "C,N,O,F",
        "max_heavy_atoms": 9,
    }
})
