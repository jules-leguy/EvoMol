import os

from evomol import run_model


def run_EvoMol(i):
    run_model({
        "obj_function": "guacamol_v2",
        "optimization_parameters": {
            "max_steps": 20,
            "pop_max_size": 1000,
            "guacamol_init_top_100": True
        },
        "io_parameters": {
            "model_path": "42_guacamol_parallel/" + str(i),
            "record_all_generated_individuals": False,
            "save_n_steps": 20,
            "silly_molecules_reference_db_path": os.environ[
                                                     "DATA"] + "/00_datasets/ChEMBL25/complete_ChEMBL_ecfp4_dict.json",
            "smiles_list_init_path": os.environ["DATA"] + "/00_datasets/Guacamol/guacamol_v1_best.smiles"

        },
        "action_space_parameters": {
            "max_heavy_atoms": 50,
            "sillywalks_threshold": 0.0,
        }
    })


run_EvoMol(1)
