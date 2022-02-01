from evomol import run_model

eval_fun = "qed"

run_model({
    "obj_function": eval_fun,
    "optimization_parameters": {
        "pop_max_size": 10,
        "k_to_replace": 2,
        "max_steps": 10,
        "selection": "random"
    },
    "io_parameters": {
        "model_path": '4_test_selection/'
    }
})
