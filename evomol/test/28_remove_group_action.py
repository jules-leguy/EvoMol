from networkx import node_connected_component
from rdkit.Chem import MolFromSmiles
import networkx as nx

from evomol import generic_action_space
from evomol.molgraphops.exploration import random_neighbour
from evomol.molgraphops.molgraph import MolGraph, MolGraphBuilder

smi = "C(NCF)=O"

graph = MolGraph(MolFromSmiles(smi))

action_spaces, parameters = generic_action_space(["C", "N", "O", "F"], 9, append_atom=False,
                                                 remove_atom=False, change_bond=False,
                                                 substitution=False, cut_insert=False, move_group=False,
                                                 remove_group=True)

builder = MolGraphBuilder(parameters, action_spaces, graph)

graph.draw()

for i in range(30):
    new_graph, _ = random_neighbour(builder, 1, uniform_action_type=True, return_mol_graph=True)
    new_graph.draw()