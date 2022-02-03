from evomol import run_model

model_path = "34_n_jobs_dft"

run_model({
    "obj_function": "homo",
    "optimization_parameters": {
        "max_steps": 20},
    "io_parameters": {
        "model_path": model_path,
        "record_history": True,
        "dft_n_jobs": 4,
        "dft_method": "test",
        "dft_base": "3-21G*",
        "dft_mem_mb": 32000,
        "dft_working_dir": "/home/jleguy/Documents/these/prod/prog/evomol/evomol/test/tests/dft_test"
    }
})
