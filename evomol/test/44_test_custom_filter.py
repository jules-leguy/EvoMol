from rdkit.Chem import MolFromSmiles

from evomol import run_model
from evomol.molgraphops.molgraph import MolGraph


def filter_fun(smiles):
    return MolGraph(MolFromSmiles(smiles)).get_n_atoms() <= 4


run_model({
    "obj_function": "qed",
    "action_space_parameters": {
        "custom_filter_function": filter_fun
    },
    "io_parameters": {
        "model_path": "44_test_custom_filter",
        "save_n_steps": 1,
        "record_all_generated_individuals": True
    }
})
