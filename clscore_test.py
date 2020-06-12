from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles, AllChem
import pickle


path_to_db_shingles = "/home/jleguy/Documents/these/prod/lib/GDBChEMBL/shingle_libs/chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl"

radius = 3
rooted = True
weighted = True
cut_off = 0.0

with open(path_to_db_shingles, "rb") as pyc:
    db_shingles = pickle.load(pyc)

def extract_shingles(smi):
    mol = MolFromSmiles(smi)
    qry_shingles = set()

    radius_constr = radius + 1

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

            if rooted:
                new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol, list(atoms), bonds, 0, 0,
                                                                  False, False, atm_idx, True, False, False)
            else:
                new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol, list(atoms), bonds, 0, 0,
                                                                  False, False, -1, True, False, False)

            qry_shingles.add(new_shingle)

    return qry_shingles

def compute_clscore(smi):
    """
    Based on https://github.com/reymond-group/GDBChEMBL
    :param individual:
    :return:
    """

    # Insuring SMILES is aromatic
    smi = MolToSmiles(MolFromSmiles(smi))

    # Extracting shingles
    qry_shingles = extract_shingles(smi)

    # calculate shingle count averaged score
    avg_score = 0
    if qry_shingles:
        sum_scores = 0
        # using log10 of shingle frequency
        if weighted:
            for shingle in qry_shingles:
                # if key not present, add 0 per default
                sum_scores += db_shingles.get(shingle, 0)
        # working binary (i.e. if present -> count++ )
        else:
            for shingle in qry_shingles:
                if shingle in db_shingles:
                    sum_scores += 1
        avg_score = sum_scores / len(qry_shingles)

    if cut_off == 0.0 or cut_off <= avg_score:
        return avg_score
