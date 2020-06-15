# EvoMol

## Install

To install EvoMol, run the following commands in your terminal.

```bash
$ git clone https://github.com/jules-leguy/EvoMol.git     # Clone EvoMol
$ cd EvoMol                                               # Move into EvoMol directory
$ conda env create -f evomol_env.yml                      # Create conda environment
$ conda activate evomolenv                                # Activate environment
$ conda install -c conda-forge rdkit                      # Install RDKit using conda-forge
$ pip install .                                           # Instal EvoMol
```


## Quickstart

Launching a QED optimization for 500 steps. Beware, you need to activate the evomolenv conda environment when you use EvoMol.

```python
from evomol import run_model
run_model({
    "obj_function": "qed",
    "search_parameters": {
        "max_steps": 500
    },
    "io_parameters": {
        "model_path": "examples/1_qed"
    },
})
```

## Model parametrization

To run a model, you need to pass a dictionary describing the run to the run_model function. This dictionary can have up to 4 entries that are described in this section.

### Objective function

The ```"obj_function"``` entry can take the following values.
* Implemented functions (see article) : "qed", "plogp", "norm_plogp", "sascore", "norm_sascore", "clscore", "homo", "lumo".
* A custom function evaluating a SMILES.
* A dictonary describing a multiobjective function containing the following entries.
    * ```"type"``` : "linear_combination" (linear combination of the properties) or "product_sigm_lin" (product of the properties after passing a linear function and a sigmoid function).
    * ```"functions"``` : list of functions (string keys describing implemented functions or custom functions).
    * Specific to the linear combination.
        * ```"coefs"``` : list of coefficients.
    * Specific to the product of sigmoid/linear functions
        * ```"a"``` list of *a* coefficients for the *ax+b* linear function definition.
        * ```"b"``` list of *b* coefficients for the *ax+b* linear function definition.
        * ```"lambda"``` list of λ coefficients for the sigmoid function definition.

