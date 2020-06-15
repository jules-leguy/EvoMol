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

To check if the installation is a success, run the following commands.

```
ipython
import evomol
```

No error should appear.