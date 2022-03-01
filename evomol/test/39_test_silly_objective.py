import os

from evomol import run_model

run_model({
    "obj_function": {
        "type": "linear_combination",
        "functions": [
            "qed",
            "sillywalks_proportion"
        ],
        "coef": [1, -1]
    },
    "io_parameters": {
        "silly_molecules_reference_db_path": os.environ[
                                                 "DATA"] + "/00_datasets/ChEMBL25/complete_ChEMBL_ecfp4_dict.json",
        "model_path": "39_test_silly_objective"
    }
})
