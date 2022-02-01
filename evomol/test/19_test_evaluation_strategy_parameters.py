from evomol import run_model, OPTEvaluationStrategy
import os

model_path = "19_test_evaluation_strategy_parameters"

QM9_path = os.environ["DATA"] + "/00_datasets/DFT/QM9/QM9_smi_datasets/train_test_dataset_100000.smi"


init_smi = []

json_cache_main = [os.environ["DATA"] + "/00_datasets/DFT/QM9/cache_QM9.json",
                   os.environ["DATA"] + "/00_datasets/DFT/cache_OD9_step7.json",
                   os.environ["DATA"] + "/00_datasets/DFT/cache_OPT_OD9_Marta_filtered.json",
                   os.environ["DATA"] + "/00_datasets/DFT/cache_OPT_OD9_0.json",
                   os.environ["DATA"] + "/00_datasets/DFT/cache_OPT.json"]

with open(QM9_path, "r") as f:
    smiles_list = f.readlines()[:300]
    for smi in smiles_list:
        init_smi.append(smi.rstrip())

run_model({
    "obj_function": OPTEvaluationStrategy(
                            prop="homo",
                            n_jobs=10,
                            cache_files=json_cache_main,
                            working_dir_path="/tmp/test_19",
                            MM_program="rdkit",
                            cache_behaviour="compute_again_delete_files"
                          ),
    "optimization_parameters": {
        "max_steps": 10,
        "pop_max_size": 300,
        "k_to_replace": 2,
        "mutable_init_pop": False
    },
    "io_parameters": {
        "model_path": model_path,
        "record_history": True,
        "smiles_list_init": init_smi,
        "evaluation_strategy_parameters": {
            "evaluate_init_pop": {"cache_behaviour": "retrieve_OPT_data"},
            "evaluate_new_solution": {"cache_behaviour": "compute_again_delete_files"}
        }
    }
})
