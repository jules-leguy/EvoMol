from evomol import run_model

run_model(
    {
        "obj_function": "isomer_S1F6",
        "optimization_parameters":{
            "max_steps": 150
        },
        "io_parameters": {
            "model_path": "test/21_valence_6",
            "save_n_steps": 1
        }
    }
)

run_model(
    {
        "obj_function": "isomer_S1F6",
        "optimization_parameters":{
            "max_steps": 150
        },
        "action_space_parameters":{
            "sulfur_valence": 2
        },
        "io_parameters": {
            "model_path": "test/21_valence_2",
            "save_n_steps": 1
        }
    }
)