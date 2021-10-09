from evomol import run_model

run_model({
    "obj_function": {
        "type": "mean",
        "functions": ["qed", "sascore"]
    },
    "io_parameters": {
        "smiles_list_init": ["C", "CN", "CNF"],
        "record_all_generated_individuals": True
    }
})