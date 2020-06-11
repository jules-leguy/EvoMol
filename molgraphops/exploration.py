import numpy as np


def _compute_root_node_id():
    return ""


def _compute_new_edge_name(action_coords):
    """
    Computing the name of the edge created by applying the given action to the given state of the molecular graph
    The edge name is the concatenation of the action type and the action id
    @param action_coords:
    @return:
    """
    return str(action_coords[0]) + "-" + str(action_coords[1])


def _compute_new_node_id(parent_node_id, action_coords):
    """
    Computing the identifier of a node from the action coordinates and the identifier of its parent.
    The node id is the concatenation of the id of its parent and the name of its edge with its parent (action)
    :param parent_node_id:
    :param action_coords
    :return:
    """

    if parent_node_id == _compute_root_node_id():
        separator = ""
    else:
        separator = "_"

    return parent_node_id + separator + _compute_new_edge_name(action_coords)


def random_neighbour(molgraph_builder, depth, return_mol_graph=False, uniform_action_type=False):
    """
    Computing a random neighbour of depth level
    Returning a tuple (SMILES, id) or (mol. graph, id)
    Raising an exception if no neighbour of given depth was found
    @param depth:
    @param molgraph_builder:
    @param return_mol_graph: whether to return a QuMolGraph object or a SMILES
    @param uniform_action_type: If true, the action type is drawn with a uniform law before the action is drawn. If
    false, the action is drawn directly with uniform law among all possible actions
    """

    # Initialization of molecular graph ID
    id = _compute_root_node_id()

    # Copying QuMolGraphBuilder
    molgraph_builder = molgraph_builder.copy()

    for i in range(depth):

        # Valid action list initialization
        valid_action_coords_list = []

        if uniform_action_type:

            # Drawing the action space
            action_space_k = np.random.choice(molgraph_builder.get_action_spaces_keys())

            # Computing the mask of valid actions
            action_space_mask = molgraph_builder.get_valid_mask_from_key(action_space_k)

            # Extracting the set of valid actions
            valid_actions = np.nonzero(action_space_mask)

            # Creating the set of valid action coords
            for valid_act in valid_actions[0]:
                valid_action_coords_list.append((action_space_k, int(valid_act)))

        else:
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

        if valid_action_coords_list:

            # Drawing action to apply
            action_coords = valid_action_coords_list[np.random.choice(np.arange(len(valid_action_coords_list)))]

            # Updating molecule ID
            id = _compute_new_node_id(id, action_coords)

            # Applying action
            molgraph_builder.execute_action_coords(action_coords)

    if return_mol_graph:
        return molgraph_builder.qu_mol_graph, id
    else:
        return molgraph_builder.qu_mol_graph.to_aromatic_smiles(), id
