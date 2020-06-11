from chempopalg.molgraphops.actionspace import ActionSpace, AddAtomActionSpace, SubstituteAtomActionSpace, \
    RemoveAtomActionSpace, ChangeBondActionSpace, CutAtomV2ActionSpace, InsertCarbonAtomV2ActionSpace, \
    MoveFunctionalGroupActionSpace


def generic_action_space(atom_symbols_list, max_heavy_atoms, substitution=True, cut_insert=True, move_group=True):
    """
    Building the action space for given atom list and max molecular size
    :param atom_symbols_list: list of atoms symbols
    :param max_heavy_atoms: max number of heavy atoms in the solutions
    :param substitution: whether the substitution action is active
    :param cut_insert: whether cut atom and insert carbon atom actions are active
    :param move_group: whether move group action is active
    :return: tuple (action_spaces, action_spaces_parameters)
    """

    accepted_atoms = atom_symbols_list
    print("SYMBOLS LIST : " + str(accepted_atoms))

    accepted_substitutions = {}

    for accepted_atom in accepted_atoms:
        if accepted_atom != "H":
            curr_atom_subst = []
            for other_accepted_atom in accepted_atoms:
                if other_accepted_atom != accepted_atom:
                    curr_atom_subst.append(other_accepted_atom)
            accepted_substitutions[accepted_atom] = curr_atom_subst

    parameters = ActionSpace.ActionSpaceParameters(max_heavy_atoms=max_heavy_atoms,
                                                   accepted_atoms=accepted_atoms,
                                                   accepted_substitutions=accepted_substitutions)

    action_spaces = [
        AddAtomActionSpace(keep_connected=True, check_validity=False),
        RemoveAtomActionSpace(keep_connected=True, check_validity=False),
        ChangeBondActionSpace(keep_connected=True, check_validity=True)]

    if substitution:
        action_spaces.append(SubstituteAtomActionSpace(check_validity=False))

    if cut_insert:
        action_spaces.extend([
            CutAtomV2ActionSpace(check_validity=False),
            InsertCarbonAtomV2ActionSpace(check_validity=False)
        ])

    if move_group:
        action_spaces.append(MoveFunctionalGroupActionSpace(check_validity=False))

    return action_spaces, parameters


