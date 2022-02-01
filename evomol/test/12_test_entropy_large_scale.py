from random import randint, choice

from evomol import EntropyContribEvaluationStrategy
from evomol.molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles


def gen_rand_ind(smi_size_max):
    smi_size = randint(1, smi_size_max)
    smi = ""
    for i in range(smi_size):
        smi += choice(["C", "O", "N", "P", "S"])

    return smi


def gen_random_pop(pop_size, smi_size_max):
    pop = []

    for i in range(pop_size):
        pop.append(MolGraph(MolFromSmiles(gen_rand_ind(smi_size_max))))

    return pop


def loop_test(pop_size, mol_size_max, n_tests_per_loop, n_loops, max_desc):

    for loop in range(n_loops):

        print("loop : " + str(loop))

        pop = gen_random_pop(pop_size, mol_size_max)
        s = EntropyContribEvaluationStrategy(max_desc, pop_size, descriptor_key="atoms")

        s.compute_record_scores_init_pop(pop)

        for i in range(n_tests_per_loop):


            ind_to_replace_idx = randint(0, pop_size-1)
            ind_to_add = gen_rand_ind(mol_size_max)

            print("desc count : " + str(s.desc_count))

            print("Replacing " + str(pop[ind_to_replace_idx]) + " by " + str(ind_to_add))

            true_delta = s.true_delta_ent(ind_to_replace_idx, MolGraph(MolFromSmiles(ind_to_add)))
            computed_delta = -s.scores[ind_to_replace_idx] + s.evaluate_individual(MolGraph(MolFromSmiles(ind_to_add)), ind_to_replace_idx)[0]


            print("loss replace : " + str(s.scores[ind_to_replace_idx]))
            print("gain add : " + str(s.evaluate_individual(MolGraph(MolFromSmiles(ind_to_add)), ind_to_replace_idx)[0]))

            print("true delta : " + str(true_delta))
            print("computed delta : " + str(computed_delta))

            if abs(true_delta - computed_delta) > 0.0001:
                print("Error !")
                exit(1)


loop_test(pop_size=100000, mol_size_max=10, n_tests_per_loop=100, n_loops=10, max_desc=500000)
