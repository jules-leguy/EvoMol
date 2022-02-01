from evomol import run_model

run_model({
    "obj_function": {
        "type": "product",
        "coef": [1, 1],
        "functions": ["qed",
                      {
                          "type": "linear_combination",
                          "coef": [1, 1],
                          "functions": ["plogp", "sascore"]
                      }]
    },
    "optimization_parameters": {
        "max_steps": 500
    },
    "io_parameters": {
        "model_path": "examples/6_test_strategies"
    },
})
