from evomol import run_model

run_model({
    "obj_function": "qed",
    "io_parameters": {
        "model_path": "37_results0"
    },
    "optimization_parameters": {
        "max_steps": 150
    }
})

run_model({
    "obj_function": "qed",
    "io_parameters": {
        "model_path": "37_results1",
        "smiles_list_init_path": "37_results0/pop.csv"

    },
    "optimization_parameters": {
        "max_steps": 1,
    }
})

run_model({
    "obj_function": "qed",
    "io_parameters": {
        "model_path": "37_results2",
        "smiles_list_init_path": "37_results0/test_smiles.smi"
    },
    "optimization_parameters": {
        "max_steps": 1,
    }
})

run_model({
    "obj_function": "qed",
    "io_parameters": {
        "model_path": "37_results3",
        "smiles_list_init": ["c1ccccc1"]
    },
    "optimization_parameters": {
        "max_steps": 1,
    }
})
