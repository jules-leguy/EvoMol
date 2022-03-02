from evomol import run_model

run_model({
    "obj_function": {
        "type": "linear_combination",
        "functions": ["homo",
                      "entropy_ifg",
                      {
                          "type": "mean",
                          "functions": ["qed", "sascore"]

                      }],
        "coef": [1, 1, 1]
    },
    "io_parameters": {
        "model_path": "38_test_record_time",
        "record_all_generated_individuals": True,
        "smiles_list_init": ["C", "N"]
    },
    "optimization_parameters": {
        "max_steps": 5
    }
})
