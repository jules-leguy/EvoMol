from rdkit.Chem.rdmolfiles import MolFromSmiles

from evomol import KRandomGraphOpsImprovingMutationStrategy, run_model
from evomol.evaluation import IsomerGuacaMolEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph, MolGraphBuilder


run_model(
    {
        "obj_function": "isomer_S1F6",
        "optimization_parameters":{
            "max_steps": 250
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
            "max_steps": 250
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