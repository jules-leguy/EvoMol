from rdkit.Chem import Lipinski, MolFromSmiles


def hetero_atoms_proportion(smiles):
    return Lipinski.NumHeteroatoms(MolFromSmiles(smiles)) / Lipinski.HeavyAtomCount(MolFromSmiles(smiles))


objective_tree = {
    # The objective function is a linear combination of three sub-functions
    "type": "linear_combination",
    "coef": [0.5, 0.25, 0.25],
    "functions": [
        {
            # Centering a Gaussian function on the 0.8 QED value
            "type": "gaussian",
            "function": "qed",
            "mu": 0.8,
            "sigma": 0.5,
            "normalize": True
        },
        {
            # The parameters of the sigmoid function are chosen so that the function value is above 0.99 when the CLScore value is above 4
            "type": "sigm_lin",
            "function": "clscore",
            "a": -1,
            "b": 3.5,
            "lambda": 10
        },
        {
            # Centering a Gaussian function on the 0.7 proportion of hetero atoms
            "type": "gaussian",
            "function": {
                # Returning 1 - (proportion of heteroatoms) to represent the proportion of carbon atoms
                "type": "one_minus",
                "function": (hetero_atoms_proportion, "hetero_atoms_proportion")
            },
            "mu": 0.7,
            "sigma": 0.5,
            "normalize": True
        },
    ]
}

from evomol import run_model

run_model({
    "obj_function": objective_tree,
    "io_parameters": {
        "model_path": "examples/multiobjective_run"
    }
})
