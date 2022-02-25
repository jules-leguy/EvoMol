from evomol import run_model

run_model({
    "obj_function": {
        "type": "abs_difference",
        "functions": ["qed", "sascore"]
    }
})
