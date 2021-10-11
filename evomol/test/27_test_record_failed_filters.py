from evomol import run_model
import os

run_model({
    "obj_function": "qed",
    "action_space_parameters":{
        "sillywalks_threshold": 0.1,
        "sascore_threshold": 4.4
    },
    "io_parameters": {
        "model_path": "27_test_failed_filters",
        "save_n_steps": 1,
        "silly_molecules_reference_db_path": os.environ["DATA"] + "/00_datasets/ChEMBL25/complete_ChEMBL_ecfp4_dict.json",
        "record_all_generated_individuals": True
    }
}
)