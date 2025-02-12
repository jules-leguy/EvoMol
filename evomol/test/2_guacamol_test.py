from evomol import run_model

run_model({
    "obj_function": "guacamol_v2",
    "optimization_parameters": {
        "max_steps": 50,
        "pop_max_size": 1000,
        "guacamol_init_top_100": False},
    "io_parameters": {
        "model_path": "2_guacamol/"
    },
})
