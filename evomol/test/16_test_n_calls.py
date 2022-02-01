from evomol import run_model

run_model({
    "obj_function": "qed",
    "optimization_parameters": {
        "max_steps": float("inf"),
        "max_obj_calls": 500
    },
    "io_parameters": {
        "model_path": "16_test_n_calls_stop"
    },
})
