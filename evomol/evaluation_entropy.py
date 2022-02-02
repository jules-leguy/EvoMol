import os
import sys
import tempfile

import numpy as np
from rdkit import RDPaths, Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdmolfiles import MolToSmiles, MolFromSmiles

from .evaluation import EvaluationStrategy, EvaluationError

ifg_path = os.path.join(RDPaths.RDContribDir, 'IFG')
sys.path.append(ifg_path)
import ifg


class EntropyContribEvaluationStrategy(EvaluationStrategy):
    """
    Scoring the individuals with their contribution to complete population's entropy. The entropy is computed on the
    extended functional groups or generic scaffolds.
    WARGNING : In this implementation, the entropy contribution of individuals evaluated with the evaluate_individual
    method is approximated using the values computed in the last call to evaluate_population
    """

    def __init__(self, n_max_desc, pop_size_max, descriptor_key):

        super().__init__()
        self.pop_size_max = pop_size_max
        self.curr_pop_size = 0

        # Initialization of population descriptors list
        # self.pop_desc_ids_list = [list() for i in range(self.pop_size_max)]
        self.pop_desc_ids_list = list(np.full((pop_size_max,), None))

        # Initialization of the dictionary associating each descriptor to a unique integer
        self.desc_id_dict = {}

        # Initialization of the array containing the number of occurrences of each descriptor
        self.desc_count = np.zeros((n_max_desc,))

        # Recording descriptor key
        self.descriptor_key = descriptor_key

        # Initialization of the attribute containing the last assigned descriptor id
        self.next_descriptor_id = 0

        # Initialization of attributes containing data for entropy computation
        self.entropy = None
        self.pop_desc_contrib = None
        self.pop_desc_minus_one_contrib = None
        self.pop_desc_plus_one_contrib = None
        self.scores = None

    def keys(self):
        return ["entropy_" + self.descriptor_key]

    def end_step_population(self, pop):
        """
        End of step : computing the entropy on complete population
        """
        self.compute_cache_record_scores()

    def compute_cache_record_scores(self):
        """
        Computing entropy contributions cache and scores of all individuals of the population (loss of entropy if they
        were removed).
        """

        print("N active descriptors : " + str(len(self.desc_count[self.desc_count > 0])))
        print("N descriptors : " + str(len(self.desc_id_dict.keys())))
        print("Curr pop size : " + str(self.curr_pop_size))

        if len(self.desc_count[self.desc_count < 0]) > 0:
            print("ERROR : negative descriptor count")
            raise EvaluationError("Negative descriptor count")

        # Computing entropy of descriptors
        self.pop_desc_contrib = self.pop_entropy_per_desc(self.desc_count, self.pop_size_max)

        # Computing entropy of population
        self.entropy = self.pop_desc_contrib.sum()

        # Computing entropy of descriptors with all the descriptors removed once
        self.pop_desc_minus_one_contrib = self.pop_entropy_per_desc(self.desc_count - 1, self.pop_size_max)

        # Computing entropy of descriptors with all the descriptors added once
        self.pop_desc_plus_one_contrib = self.pop_entropy_per_desc(self.desc_count + 1, self.pop_size_max)

        # Computing scores for all individuals
        self.scores = []
        for idx_ind in range(self.pop_size_max):

            curr_ind_desc = self.pop_desc_ids_list[idx_ind]

            # Computing scores only for defined individuals
            if curr_ind_desc is not None:
                loss_remove = self.loss_remove_faster(curr_ind_desc)
                self.scores.append(loss_remove)

    def compute_record_scores_init_pop(self, population):
        """
        Computing and recording entropy of complete population
        """

        # Computing descriptors
        self.compute_descriptors(population)

        # Computing cache and scores
        self.compute_cache_record_scores()

    def record_ind_score(self, idx, new_total_score, new_scores, new_individual):

        # Updating the descriptors of the given individual
        self.update_descriptors(idx, new_individual)

        # Call to super class to update the score of the given individual. Note : useless as all population scores will
        # be computed at the beginning of the next step
        super(EntropyContribEvaluationStrategy, self).record_ind_score(idx, new_total_score, new_total_score, new_individual)

    def evaluate_individual(self, individual, to_replace_idx=None):
        """
        Estimation of the delta of entropy caused by the given individual.
        """

        super().evaluate_individual(individual, to_replace_idx)

        # Extracting the descriptors of the given individual
        ind_to_add_desc = self.extract_descriptors(individual)

        ind_desc_ids = []

        # Extracting the descriptors ids of evaluated individual
        for curr_ind_desc in ind_to_add_desc:
            ind_desc_ids.append(self.get_desc_id(curr_ind_desc))

        # Extracting the descriptors ids of individual to be replaced
        if self.pop_desc_ids_list[to_replace_idx] is not None:
            ind_to_replace_desc_ids = self.pop_desc_ids_list[to_replace_idx]
        else:
            ind_to_replace_desc_ids = []

        # Computing entropy gain
        entropy_gain = self.gain_replace_faster(ind_desc_ids, ind_to_replace_desc_ids)

        return entropy_gain, [entropy_gain]

    def loss_remove_faster(self, ind_desc):
        """
        Computing the entropy loss that would occur in the population of same size if the given descriptors were removed
        """
        return self.pop_desc_contrib[ind_desc].sum() - self.pop_desc_minus_one_contrib[ind_desc].sum()

    def gain_replace_faster(self, ind_to_add_desc, ind_to_remove_desc):
        """
        Computing the gain of entropy that would occur by adding the descriptors of ind_to_add_desc into the population.
        If there is an intersection between the descriptors that are added and removed (ind_to_remove_desc), it is taken
        in account separately so that the computation is based on actual descriptors count.
        """

        intersect_desc = list(set(ind_to_add_desc) & set(ind_to_remove_desc))
        to_add_minus_to_remove_desc = list(set(ind_to_add_desc) - set(ind_to_remove_desc))

        return self.loss_remove_faster(intersect_desc) + self.pop_desc_plus_one_contrib[
                   to_add_minus_to_remove_desc].sum() - self.pop_desc_contrib[to_add_minus_to_remove_desc].sum()

    def get_active_and_removed_descriptors(self):
        removed_descriptors = []
        active_descriptors = []
        for desc, id in self.desc_id_dict.items():
            if self.desc_count[id] > 0:
                active_descriptors.append(desc)
            else:
                removed_descriptors.append(desc)
        return active_descriptors, removed_descriptors

    def pop_entropy_per_desc(self, desc_count, pop_size):
        """
        Computing the entropy of the descriptors
        :return:
        """
        h = np.zeros(desc_count.shape)

        mask = desc_count > 0
        A = desc_count[mask]
        A = A / float(pop_size)
        h[mask] = -(A * np.log(A))
        return h

    def true_delta_ent(self, idx_remove, ind_to_add):

        start_entropy = self.entropy
        desc_id_to_remove = self.pop_desc_ids_list[idx_remove]
        desc_to_add = self.extract_descriptors(ind_to_add)
        desc_count = self.desc_count.copy()

        for desc_id in desc_id_to_remove:
            desc_count[desc_id] -= 1
        for desc in desc_to_add:
            desc_count[self.get_desc_id(desc)] += 1

        end_entropy = self.pop_entropy_per_desc(desc_count, self.pop_size_max).sum()

        return end_entropy - start_entropy

    def extract_descriptors(self, individual):
        """
        Returning the descriptor(s) extracted from the given individual
        :param individual:
        :return:
        """

        if self.descriptor_key == "gen_scaffolds":
            return [MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MolFromSmiles(individual.to_smiles())))]
        elif self.descriptor_key == "ifg":
            curr_ifgs = ifg.identify_functional_groups(MolFromSmiles(individual.to_smiles()))
            return list(set([curr_ifg[2] for curr_ifg in curr_ifgs]))
        elif self.descriptor_key == "atoms":
            return list(set(individual.get_atom_types()))
        elif self.descriptor_key == "shg_1":
            return list(extract_shingles(individual.to_aromatic_smiles(), 1))
        elif self.descriptor_key == "checkmol":
            return list(set(extract_checkmol(individual)))

    def get_desc_id(self, desc):
        """
        Returning the identifier of the given descriptor.
        If the descriptor is not known, assigning it new unique descriptor
        :param desc:
        :return:
        """

        if desc in self.desc_id_dict:
            return self.desc_id_dict[desc]
        else:
            self.desc_id_dict[desc] = self.next_descriptor_id
            self.next_descriptor_id += 1
            return self.desc_id_dict[desc]

    def compute_descriptors(self, population):
        """
        Computing the internal descriptors for the given new population
        :param population:
        :return:
        """

        # Iterating over the individuals of the given population
        for ind_idx, ind in enumerate(population):
            if ind is not None:
                self.update_descriptors(ind_idx, ind)

    def update_descriptors(self, idx, new_ind):
        """
        Updating the internal descriptors for given new individual with given index
        :param idx:
        :param new_ind:
        :return:
        """

        # Updating population size if the new individual replaces an undefined one
        if self.pop_desc_ids_list[idx] is None and new_ind is not None:
            self.curr_pop_size += 1

        # Updating the population size if the new individual is not defined
        if new_ind is None and self.pop_desc_ids_list[idx] is not None:
            self.curr_pop_size -= 1

        # Removing data for existing former individual
        if self.pop_desc_ids_list[idx] is not None:
            for desc_id in self.pop_desc_ids_list[idx]:
                self.desc_count[desc_id] -= 1
            self.pop_desc_ids_list[idx] = None

        # Recording its descriptors if the new individual is defined
        if new_ind is not None:

            # Extracting individual's descriptors
            curr_ind_desc_ids = []
            for curr_desc in self.extract_descriptors(new_ind):
                curr_desc_id = self.get_desc_id(curr_desc)
                self.desc_count[curr_desc_id] += 1
                curr_ind_desc_ids.append(curr_desc_id)

            # Saving ids of descriptors for given individual
            self.pop_desc_ids_list[idx] = curr_ind_desc_ids

    def get_additional_population_scores(self):
        local_d = {
            self.keys()[0] + "_active_desc": len(self.desc_count[self.desc_count > 0]),
            self.keys()[0] + "_total_desc": self.next_descriptor_id,
            self.keys()[0] + "_pop_entropy": self.entropy
        }

        local_d.update(super().get_additional_population_scores())

        return local_d


def extract_checkmol(molgraph):
    """
    Extracting checkmol descriptors from given molecular graph.
    see https://homepage.univie.ac.at/norbert.haider/cheminf/cmmm.html
    """

    obabel_cmd = "obabel" + " \"-:%s\" -omol -O %s 2>/dev/null"
    checkmol_cmd = os.getenv("CHECKMOL_EXE") + " %s > %s"

    smiles = molgraph.to_aromatic_smiles()

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.mol', delete=True) as mol_fic:
        fic_name = mol_fic.name
        os.system(obabel_cmd % (smiles, fic_name))
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.res', delete=True) as mol_ficout:
            ficout_name = mol_ficout.name
            os.system(checkmol_cmd % (fic_name, ficout_name))
            lines = [l.strip() for l in mol_ficout.readlines()]
        if len(lines) == 0:
            return lines
        elif lines[0] == "unknown query file format!":
            return []
    return lines


def extract_shingles(smiles, level, as_list=False):
    """
    Extracting up to the given level from the given smiles
    see https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0321-8
    """

    if as_list:
        qry_shingles = list()
    else:
        qry_shingles = set()

    radius_constr = level + 1

    # Reloading molecule to make it aromatic
    mol = MolFromSmiles(smiles)

    for atm_idx in range(mol.GetNumAtoms()):
        for N in range(1, radius_constr):
            bonds = AllChem.FindAtomEnvironmentOfRadiusN(mol, N, atm_idx)

            if not bonds:
                break

            # the reportedly faster method
            atoms = set()
            for bond_id in bonds:
                bond = mol.GetBondWithIdx(bond_id)
                atoms.add(bond.GetBeginAtomIdx())
                atoms.add(bond.GetEndAtomIdx())

            # Computed rooted shingle
            new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol, list(atoms), bonds, 0, 0,
                                                              False, False, atm_idx, True, False, False)
            if as_list:
                qry_shingles.append(new_shingle)
            else:
                qry_shingles.add(new_shingle)

    return qry_shingles
