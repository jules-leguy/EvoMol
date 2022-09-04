from evomol.molgraphops.default_actionspaces import generic_action_space

action_spaces, action_spaces_parameters = \
    generic_action_space(atom_symbols_list=["C", "N", "O", "F"],
                         max_heavy_atoms=20,
                         append_atom=False,
                         remove_atom=False,
                         change_bond=False,
                         change_bond_prevent_breaking_creating_bonds=False,
                         substitution=False,
                         cut_insert=False,
                         move_group=True,
                         remove_group=False,
                         remove_group_only_remove_smallest_group=True)

SULFUR_VALENCE = 6

from evomol.molgraphops.molgraph import MolGraphBuilder, MolGraph
from rdkit.Chem import MolFromSmiles
import numpy as np


def enumerate_neighbours(smiles, graph_actions_depth, action_spaces, action_spaces_parameters):
    if graph_actions_depth <= 0:
        return []
    else:

        # Initialization of the list of all neighbour SMILES
        all_smiles = []

        # Initialiation of the list that contains all valid actions for current SMILES
        valid_action_coords_list = []

        molgraph_builder = MolGraphBuilder(action_spaces_parameters, action_spaces,
                                           MolGraph(MolFromSmiles(smiles), sulfur_valence=SULFUR_VALENCE))

        # Computing valid actions
        valid_action_dict = molgraph_builder.get_action_spaces_masks()

        # Iterating over the actions of the different action spaces
        for key, validity in valid_action_dict.items():

            # Recording the id of the valid actions for the current action space
            curr_key_valid_actions = np.nonzero(validity)

            # Iterating over the valid actions for the current action space
            for curr_key_valid_act in curr_key_valid_actions[0]:
                # Adding the current valid action to the list
                valid_action_coords_list.append((key, int(curr_key_valid_act)))

        # Iterating over all valid actions
        for curr_valid_act_coords in valid_action_coords_list:
            # Copying QuMolGraphBuilder
            molgraph_builder_copy = molgraph_builder.copy()

            # Application of current action
            molgraph_builder_copy.execute_action_coords(curr_valid_act_coords)

            # Computing the SMILES of the new neighbour
            new_neighbour_smiles = molgraph_builder_copy.qu_mol_graph.to_aromatic_smiles()

            # Recursive call to obtain the neighbours of new SMILES (in-depth search)
            new_neigbour_neighbours_smiles_list = enumerate_neighbours(new_neighbour_smiles,
                                                                       graph_actions_depth - 1,
                                                                       action_spaces,
                                                                       action_spaces_parameters)

            # Adding new neighbour and its recursive neighbours to the list
            all_smiles.append(new_neighbour_smiles)
            all_smiles.extend(new_neigbour_neighbours_smiles_list)

    return all_smiles


smiles_list = enumerate_neighbours("c1ccc(CO)nc1", 1, action_spaces, action_spaces_parameters)
# smiles_list = enumerate_neighbours("NNC", 1, action_spaces, action_spaces_parameters)
for smi in smiles_list:
    mol = MolGraph(MolFromSmiles(smi))
    mol.draw()
