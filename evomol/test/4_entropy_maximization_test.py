from evomol import run_model

run_model({
    "obj_function": "entropy_shg_1",
    "optimization_parameters": {
        "max_steps": 100,
        "pop_max_size": 1000},
    "io_parameters": {
        "model_path": "4_entropy/"
    },
})
