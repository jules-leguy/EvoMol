from evomol import run_model

model_path = "34_n_jobs_dft"

run_model({
    "obj_function": "homo",
    "optimization_parameters": {
        "max_steps": 20},
    "io_parameters": {
        "model_path": model_path,
        "record_history": True,
        "dft_n_jobs": 4
    }
})
