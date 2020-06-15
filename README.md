# EvoMol

## Install

To install EvoMol, run the following commands.

```bash
git clone https://github.com/jules-leguy/EvoMol.git     # Clone EvoMol
cd EvoMol                                               # Move into EvoMol directory
conda env create -f evomol_env.yml                      # Create conda environment
conda activate evomolenv                                # Activate environment
conda install -c conda-forge rdkit                      # Install RDKit using conda-forge
pip install .                                           # Instal EvoMol
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


