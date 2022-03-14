import io
from os import makedirs
from os.path import dirname

import networkx as nx
import numpy as np
from PIL import Image
from rdkit.Chem import Kekulize, FastFindRings, MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdchem import RWMol, Atom, GetPeriodicTable, Bond, BondType
from rdkit.Chem.rdmolfiles import MolToSmiles
from rdkit.Chem.rdmolops import GetAdjacencyMatrix, SanitizeMol


class MolGraph:

    def __init__(self, init_state=None, sanitize_mol=False, mutability=True, sulfur_valence=6):
        """
        Constructor of the QuMolGraph class.
        If a rdchem.Mol object is given, it becomes the initial state of the molecular graph
        :param init_state: Molecule at initialization
        :@param sanitize_mol: Whether to sanitize the molecules at initialization and after changes. This can help to
        improve the uniqueness of canonical smiles (!), but can also have an effect on the location of aromatic double
        bonds.
        :param value of the mutability flag to all atoms of the initial state
        """

        # Creating or reading molecular graph
        if init_state is not None:
            self.mol_graph = RWMol(init_state)
        else:
            self.mol_graph = RWMol()

        # Recording sanitize mol choice
        self.sanitize_mol = sanitize_mol

        # Updating the internal representation
        self.update_mol_representation()

        # Setting mutability of initial state
        for id in range(self.get_n_atoms()):
            self.set_atom_mutability(id, mutability)

        # Initialization of modifications count
        self.n_modifications = 0

        self.sulfur_valence = sulfur_valence

    def get_atom_mutability(self, id):
        """
        Return the boolean value of the mutability of the atom of given id
        :param id:
        :return:
        """
        return self.mol_graph.GetAtomWithIdx(id).GetBoolProp("mutability")

    def set_atom_mutability(self, id, mutability):
        """
        Setting the property of mutability of the atom of given id to the given boolean value
        :param id:
        :param mutability:
        :return:
        """
        self.mol_graph.GetAtomWithIdx(id).SetBoolProp("mutability", mutability)

    def __repr__(self):
        return self.to_aromatic_smiles()

    def update_mol_representation(self):
        """
        Updating internal RDKIT representation of the molecular graph
        This method should be called at each action to the molecular graph
        :return:
        """

        # Sanitizing mol if needed
        if self.sanitize_mol:
            SanitizeMol(self.mol_graph)

            # Clear computed properties
            self.mol_graph.ClearComputedProps()

        # Kekulization of the molecular graph
        Kekulize(self.mol_graph)

        # Setting all atoms to non aromatics
        for i in range(self.mol_graph.GetNumAtoms()):
            self.mol_graph.GetAtomWithIdx(i).SetIsAromatic(False)

        # Setting all bonds to non aromatics
        for i in range(self.mol_graph.GetNumAtoms()):
            for j in range(self.mol_graph.GetNumAtoms()):
                bond = self.mol_graph.GetBondBetweenAtoms(i, j)
                if bond is not None:
                    bond.SetIsAromatic(False)

        # Updating the property cache of atoms
        for i in range(self.mol_graph.GetNumAtoms()):
            self.mol_graph.GetAtomWithIdx(i).UpdatePropertyCache()

        # Updating RDKit representation
        self.mol_graph.UpdatePropertyCache()
        FastFindRings(self.mol_graph)

    def copy(self):
        """
        Returning a deep copy of the current instance
        :return:
        """
        new_qumol_graph = MolGraph(sanitize_mol=self.sanitize_mol, sulfur_valence=self.sulfur_valence)
        new_qumol_graph.n_modifications = self.n_modifications
        new_qumol_graph.mol_graph = RWMol(self.mol_graph, True)
        new_qumol_graph.update_mol_representation()
        return new_qumol_graph

    def add_atom(self, type):
        """
        Adding an atom of given type at the last position in the molecule ordering
        :param type:
        :return: The index of the added atom
        """
        at = Atom(type)
        at.SetBoolProp("mutability", True)
        self.mol_graph.AddAtom(at)
        return at

    def end_atom_addition_procedure(self):
        """
        THIS METHOD MUST BE CALLED AFTER EACH ATOM ADDITION
        @return:
        """

        # Updating internal representation
        self.update_mol_representation()

    def rm_atom(self, id, update_repr=True):
        """
        Removing the atom of given id from the molecule
        :param id:
        :return:
        """

        # Extracting the list of neighbour atoms
        # neigh_at = []
        # for neigh_id in np.argwhere(self.get_adjacency_matrix()[id] == 1).reshape(-1, ):
        #     neigh_at.append(self.mol_graph.GetAtomWithIdx(int(neigh_id)))

        # Removing atom
        self.mol_graph.RemoveAtom(id)

        # # Updating RDKit representation of former neighbour atoms so that implicit H are updated
        # for neigh in neigh_at:
        #     neigh.UpdatePropertyCache()

        # Updating the internal representation
        if update_repr:
            self.update_mol_representation()

    def replace_atom(self, id, new_at_type):
        """
        Replacing the atom of given id with a new atom of given type
        :param id: id of atom that must be replaced
        :param new_at_type: type (atomic symbol) of the new atom
        :return:
        """

        # Changing atomic number
        self.mol_graph.GetAtomWithIdx(id).SetAtomicNum(GetPeriodicTable().GetAtomicNumber(new_at_type))

        # Setting formal charge to 0
        self.mol_graph.GetAtomWithIdx(id).SetFormalCharge(0)

        # Updating the internal representation
        self.update_mol_representation()

    def set_bond(self, from_at, to_at, bond_type_num, update_repr=True):
        """
        Setting the bond of given type between the two specified atoms
        """

        # Extracting current bond between given atoms
        curr_bond = self.mol_graph.GetBondBetweenAtoms(int(from_at), int(to_at))

        if bond_type_num == 0:
            bond_to_set = None
        elif bond_type_num == 1:
            bond_to_set = BondType.SINGLE
        elif bond_type_num == 2:
            bond_to_set = BondType.DOUBLE
        elif bond_type_num == 3:
            bond_to_set = BondType.TRIPLE

        if curr_bond is None:
            self.mol_graph.AddBond(int(from_at), int(to_at), bond_to_set)
        elif bond_to_set is None:
            self.mol_graph.RemoveBond(int(from_at), int(to_at))
        else:
            curr_bond.SetBondType(bond_to_set)

        if update_repr:
            # Updating internal representation
            self.update_mol_representation()

    def add_bond(self, from_at, to_at, update_repr=True):
        """
        Adding a bond between the two atoms of given ids according to the following scheme.
        None -> Simple
        Simple -> Double
        Double -> Triple
        If bonding is impossible because of valence, an exception is returned.
        :param from_at:
        :param to_at:
        :return:
        """

        # Extracting current bond between given atoms
        curr_bond = self.mol_graph.GetBondBetweenAtoms(from_at, to_at)

        if curr_bond is None:
            self.mol_graph.AddBond(from_at, to_at, BondType.SINGLE)
            # print("Adding None -> Single bond between atoms of idx " + str(from_at) + " and " + str(to_at))
        elif Bond.GetBondType(curr_bond) == BondType.SINGLE:
            # print("Adding Single -> Double bond between atoms of idx " + str(from_at) + " and " + str(to_at))
            curr_bond.SetBondType(BondType.DOUBLE)
        elif Bond.GetBondType(curr_bond) == BondType.DOUBLE:
            # print("Adding Double -> Triple bond between atoms of idx " + str(from_at) + " and " + str(to_at))
            curr_bond.SetBondType(BondType.TRIPLE)

        else:
            print("Unknown bond type : " + self.to_smiles())
            raise Exception("Unknown bond type : " + self.to_smiles())

        if update_repr:
            # Updating internal representation
            self.update_mol_representation()

    def rm_bond(self, from_at, to_at, update_repr=True):
        """
        Removing a bond between the two atoms of given ids
        :param from_at:
        :param to_at:
        :return:
        """

        # Extracting current bond
        curr_bond = self.mol_graph.GetBondBetweenAtoms(from_at, to_at)

        if curr_bond is None:
            raise Exception("Trying to remove not existing bond")
        elif Bond.GetBondType(curr_bond) == BondType.TRIPLE:
            curr_bond.SetBondType(BondType.DOUBLE)
            # print("Removing Triple -> Double bond between atoms of idx " + str(from_at) + " and " + str(to_at))
        elif Bond.GetBondType(curr_bond) == BondType.DOUBLE:
            curr_bond.SetBondType(BondType.SINGLE)
            # print("Removing Double -> Simple bond between atoms of idx " + str(from_at) + " and " + str(to_at))
        elif Bond.GetBondType(curr_bond) == BondType.SINGLE:
            # print("Removing Single -> None bond between atoms of idx " + str(from_at) + " and " + str(to_at))
            self.mol_graph.RemoveBond(from_at, to_at)
        else:
            print("Unknown bond type : " + self.to_smiles())
            raise Exception("Unknown bond type : " + self.to_smiles())

        if update_repr:
            # Updating internal representation
            self.update_mol_representation()

    def get_bridge_bonds_matrix(self):
        """
        Returning a boolean matrix of size (n_defined_atoms, n_defined_atoms) representing whether bonds of the
        molecular graph are bridges.
        """

        # Converting the molecular graph to a NetworkX object
        nx_mol_graph = nx.from_numpy_array(self.get_adjacency_matrix())

        # Initialization of the output matrix of bridge bonds
        output_bridges_matrix = np.full((self.get_n_atoms(), self.get_n_atoms()), False)

        # Extracting the list of bridges in the molecular simple graph
        bridges_list = list(nx.bridges(nx_mol_graph))

        for bridge in bridges_list:
            output_bridges_matrix[bridge[0], bridge[1]] = True
            output_bridges_matrix[bridge[1], bridge[0]] = True

        return output_bridges_matrix

    def get_articulation_points_vector(self):
        """
        Returning a boolean vector representing whether the atoms of the molecular graph are articulation points
        (vertices whose removal would create two connected components).
        :return:
        """

        # Articulation points vector initialization
        art_points_vector = np.zeros((self.get_n_atoms(),))

        # Converting the molecular graph to a NetworkX object
        nx_mol_graph = nx.from_numpy_array(self.get_adjacency_matrix())

        # Computing articulation points
        art_points_ids = nx.articulation_points(nx_mol_graph)

        # Setting output vector
        for art_points_id in art_points_ids:
            art_points_vector[art_points_id] = 1

        return art_points_vector

    def get_atom_degree(self, id, as_multigraph):
        """
        Returning the degree of the the atom of given id.
        If as_multigraph is set to False, the degree of the atom is defined as the number of atoms it shares bonds with.
        If as_multigraph is set to True, the degree of the atom considers the types of the bonds and is equivalent to
        the explicit valence.
        :return:
        """

        if as_multigraph:
            return self._get_expl_valence(id)
        else:
            return self.mol_graph.GetAtomWithIdx(id).GetDegree()

    def get_formal_charge_vector(self):
        """
        Returning a vector containing the formal charge for each defined atom
        """

        formal_charge_vector = []

        for i in range(self.get_n_atoms()):
            formal_charge_vector.append(self.get_formal_charge(i))

        return formal_charge_vector

    def get_formal_charge(self, at_idx):
        """
        Returning the formal charge of the atom of given idx
        """

        return self.mol_graph.GetAtomWithIdx(at_idx).GetFormalCharge()

    def get_max_valence(self, atom_type):
        """
        Returning max. valence for atom of given type
        :param atom_type:
        :return:
        """

        if atom_type == "S":
            return self.sulfur_valence
        else:
            return GetPeriodicTable().GetDefaultValence(GetPeriodicTable().GetAtomicNumber(atom_type))

    def _get_n_free_electrons(self, at_idx):
        """
        Returning the number of free electrons of the atom of given id
        :param at_idx:
        :return:
        """
        return self.get_max_valence(
            self.mol_graph.GetAtomWithIdx(at_idx).GetSymbol()) - self._get_expl_valence(at_idx)

    def get_free_electrons_vector(self):
        """
        Returning a vector containing the number of free electrons for each defined atom
        :return:
        """
        free_electrons = []
        for i in range(self.get_n_atoms()):
            free_electrons.append(self._get_n_free_electrons(i))
        return np.array(free_electrons)

    def get_adjacency_matrix(self):
        """
        Returning the adjacency matrix of the molecular graph (defined atoms)
        :return:
        """
        return GetAdjacencyMatrix(self.export_mol())

    def get_atom_type(self, at_idx):
        """
        Returning the type (string) of the atom of given index
        :param at_idx:
        :return:
        """
        return self.mol_graph.GetAtomWithIdx(int(at_idx)).GetSymbol()

    def get_atom_types(self):
        """
        Returning a vector containing the atom types contained in the molecule
        :return:
        """
        atom_types = []
        for i in range(self.get_n_atoms()):
            atom_types.append(self.get_atom_type(i))
        return np.array(atom_types)

    def get_bond_type_num(self, at1_idx, at2_idx):
        """
        Returning the integer bond type formed by the two atoms of given idx. If no bond exists, returning zero.
        @param at1_idx:
        @param at2_idx:
        @return:
        """

        bond = self.mol_graph.GetBondBetweenAtoms(int(at1_idx), int(at2_idx))

        if bond is None:
            return 0
        elif Bond.GetBondType(bond) == BondType.SINGLE:
            return 1
        elif Bond.GetBondType(bond) == BondType.DOUBLE:
            return 2
        elif Bond.GetBondType(bond) == BondType.TRIPLE:
            return 3
        elif Bond.GetBondType(bond) == BondType.QUADRUPLE:
            return 4

    def _get_expl_valence(self, at_idx):
        """
        Returning the number of electrons engaged in non H covalent bond of the atom of given id
        :param at_idx:
        :return:
        """
        at = self.mol_graph.GetAtomWithIdx(at_idx)
        at.UpdatePropertyCache()
        return at.GetExplicitValence()

    def _is_new_bond_possible(self, at1_idx, at2_idx):
        """
        Returning whether a new bond is possible between the two atoms of given ID. Considering the number of free
        electrons of both atoms.
        :param at1_idx:
        :param at2_idx:
        :return:
        """
        return min(self._get_n_free_electrons(at1_idx), self._get_n_free_electrons(at2_idx)) > 0

    def to_smiles(self):
        """
        Returning the SMILES version of the molecule
        :return:
        """
        return MolToSmiles(self.mol_graph)

    def to_aromatic_smiles(self):
        return MolToSmiles(MolFromSmiles(MolToSmiles(MolGraph(MolFromSmiles(self.to_smiles())).mol_graph)))
        # return MolToSmiles(MolFromSmiles(self.to_smiles()))

    def export_mol(self):
        """
        Exporting the molecule as a RDKit Molecule object
        :return:
        """
        return self.mol_graph.GetMol()

    def draw(self, at_idx=False, show=True, size=200, write_to_path=None):
        """
        Drawing the molecule
        :return:
        """
        mol = self.export_mol()
        atoms = mol.GetNumAtoms()

        # Setting the ids as a property if requested
        if at_idx:
            for idx in range(atoms):
                mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))

        # Computing coordinates and making sure the properties are computed
        Compute2DCoords(mol)
        mol.UpdatePropertyCache()

        # Drawing the molecule
        dr = rdMolDraw2D.MolDraw2DCairo(size, size)
        opts = dr.drawOptions()

        # Transparent background if not writing to file
        if write_to_path is None:
            opts.clearBackground = False

        dr.DrawMolecule(mol)
        dr.FinishDrawing()

        # Loading the molecule as a PIL object
        bytes_images = dr.GetDrawingText()
        image = Image.open(io.BytesIO(bytes_images))

        if show:
            image.show()

        if write_to_path is not None:
            # Creating directories if they don't exist
            makedirs(dirname(write_to_path), exist_ok=True)

            # Writing image to disk
            image.save(write_to_path, "PNG")

        return image

    def get_n_atoms(self):
        return self.mol_graph.GetNumAtoms()


class MolGraphBuilder:

    def __init__(self, actionspace_parameters, action_spaces, qumol_graph=None):
        """
        QuMolGraphBuilder constructor
        :param actionspace_parameters: ActionSpaceParameters object containing all the parameters needed by the
        ActionSpace objects
        :param action_spaces: list or dictionary of ActionSpace objects whose union define the action space of the
        molecular graph creation
        """

        self.parameters = actionspace_parameters

        # Molecular graph initialization
        if qumol_graph is None:
            self.qu_mol_graph = MolGraph()
        else:
            self.qu_mol_graph = qumol_graph

        if isinstance(action_spaces, dict):
            self.action_spaces_d = action_spaces
        else:
            # Mapping the ActionSpaces object on their identifier in a dictionary attribute
            self.action_spaces_d = {}
            for action_space in action_spaces:
                self.action_spaces_d[action_space.action_space_type_id()] = action_space

    def copy(self):
        """
        Returning a deep copy of the current instance
        :return:
        """

        new_qumol_builder = MolGraphBuilder(self.parameters, self.action_spaces_d)
        new_qumol_builder.qu_mol_graph = self.qu_mol_graph.copy()
        return new_qumol_builder

    def get_action_spaces_keys(self):
        """
        Returning the keys of the different action spaces
        """
        return list(self.action_spaces_d.keys())

    def get_valid_mask_from_key(self, action_space_k):
        """
        Returning the list of valid actions for the given key
        """

        return self.action_spaces_d[action_space_k].get_valid_actions_mask(self.parameters,
                                                                           self.qu_mol_graph)

    def get_action_spaces_masks(self):
        """
        Returning a dictionary of boolean vectors representing the valid action spaces for the current state of
        the molecular graph. The vectors are mapped on the identifier of the ActionSpace class they are produced by.
        :return:
        """

        action_spaces_masks_d = {}

        for action_space_k, action_space in self.action_spaces_d.items():
            action_spaces_masks_d[action_space_k] = action_space.get_valid_actions_mask(self.parameters,
                                                                                        self.qu_mol_graph)

        return action_spaces_masks_d

    def get_action_expl_str(self, action_coords):
        """
        Returning a short string giving information about the context of the given action
        """

        return self.action_spaces_d[action_coords[0]].get_action_expl(action_coords[1], self.parameters,
                                                                      self.qu_mol_graph)

    def get_action_spaces_sizes(self):
        """
        Returning a dictionary containing the sizes of the vectors representing the action spaces. The sizes are mapped
        on the identifier of the ActionSpace class they are produced by.
        :return:
        """

        action_spaces_sizes_d = {}

        for action_space_k, action_space in self.action_spaces_d.items():
            action_spaces_sizes_d[action_space_k] = action_space.get_action_space_size(self.parameters,
                                                                                       self.qu_mol_graph)

        return action_spaces_sizes_d

    def action_to_str(self, action_coord):
        """
        To string version of action of given coordinates (k, l) if allowed
        k is the identifier of the ActionSpace class that manages the type of action
        l gives the id of the action
        :param action_coord:
        :return:
        """

        return self.action_spaces_d[action_coord[0]].action_to_str(action_coord[1], self.parameters, self.qu_mol_graph)

    def _execute_action(self, action_coord):
        """
        Executing action of given coordinates (k, l) if allowed
        k is the identifier of the ActionSpace class that manages the type of action
        l gives the id of the action
        :param action_coord:
        :return:
        """
        return self.action_spaces_d[action_coord[0]].execute_action(action_coord[1], self.parameters, self.qu_mol_graph)

    def _extract_valid_actions_and_weights(self, action_weights):
        """
        Computing two lists of same sizes containing the set of valid action coordinates according to the current
        state of the molecular graph, and the set of corresponding probabilistic weights
        :param action_weights:
        :return:
        """

        # Lists initialization
        valid_actions_coordinates = []
        valid_actions_weights = []

        # Iterating over all the action spaces
        for curr_action_space_type_id, curr_action_space in self.action_spaces_d.items():
            # Extracting the the size of the current action space
            curr_action_space_size = curr_action_space.get_action_space_size(self.parameters, self.qu_mol_graph)

            # Computing the mask of valid actions in the current action space
            curr_valid_action_space_mask = curr_action_space.get_valid_actions_mask(self.parameters,
                                                                                    self.qu_mol_graph)

            # Generating a vector containing the coordinates of all the actions of the current action space
            curr_action_coordinates = np.array((np.full((curr_action_space_size,), curr_action_space_type_id),
                                                np.arange(curr_action_space_size)), dtype=np.object).T

            # Recording the set of valid action coordinates and their associated probabilistic weights for the
            # current action space
            valid_actions_coordinates.extend(list(np.array(curr_action_coordinates[curr_valid_action_space_mask])))
            valid_actions_weights.extend(
                np.array(action_weights[curr_action_space_type_id][curr_valid_action_space_mask]))

        return valid_actions_coordinates, valid_actions_weights

    def execute_action_coords(self, action_coords):

        # Updating the number of modifications
        self.qu_mol_graph.n_modifications += 1

        # Returning the action result
        return self._execute_action(action_coords)
