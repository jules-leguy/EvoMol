import numpy as np

MAX_DESC = 500


def pop_to_desc(pop):
    desc = np.zeros(MAX_DESC)

    for ind_desc in pop:
        for d in ind_desc:
            desc[d] += 1

    return desc


def entropy(pop, size):
    pop_desc = pop_to_desc(pop)

    A = pop_desc[pop_desc > 0] / float(size)
    return -((A * np.log(A)).sum())


def loss_remove_idx(pop, idx):
    size = len(pop)

    entropy_before = entropy(pop, size)
    mask = np.arange(size)
    mask = mask[mask != idx]
    pop = np.array(pop)[mask]
    entropy_after = entropy(pop, size)

    return entropy_before - entropy_after


def loss_remove_vector(pop):
    vector = np.zeros((len(pop, )))
    for i in range(len(pop)):
        vector[i] = loss_remove_idx(pop, i)

    return vector


def pop_entropy_per_desc(desc_count, pop_size):
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

#
# def extract_desc_count(pop):
#     desc_count = {}
#
#     for ind in pop:
#         for desc in ind:
#             if desc in desc_count:
#                 desc_count[desc] += 1
#             else:
#                 desc_count[desc] = 1
#
#     desc_count_array = []
#     for key in sorted(desc_count.keys()):
#         desc_count_array.append(desc_count[key])
#
#     return np.array(desc_count_array)


def loss_remove_faster(ind_desc, pop_desc_contrib, pop_desc_minus_one_contrib):
    return pop_desc_contrib[ind_desc].sum() - pop_desc_minus_one_contrib[ind_desc].sum()


def vect_loss_remove_faster(pop, pop_desc_contrib, pop_desc_minus_one_contrib):
    loss_remove_vect = []

    for ind in pop:
        loss_remove_vect.append(
            loss_remove_faster(ind, pop_desc_contrib, pop_desc_minus_one_contrib)
        )

    return loss_remove_vect


def gain_replace_faster(ind_to_add_desc, ind_to_remove_desc, pop_desc_contrib, pop_desc_minus_one_contrib,
                        pop_desc_plus_one_contrib):
    intersect_desc = list(set(ind_to_add_desc) & set(ind_to_remove_desc))
    to_add_minus_to_remove_desc = list(set(ind_to_add_desc) - set(ind_to_remove_desc))

    return loss_remove_faster(intersect_desc, pop_desc_contrib, pop_desc_minus_one_contrib) + pop_desc_plus_one_contrib[
        to_add_minus_to_remove_desc].sum() - pop_desc_contrib[to_add_minus_to_remove_desc].sum()


def test_delta(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, ind_to_remove_idx,
               ind_to_add_desc):
    loss_remove = loss_remove_faster(pop[ind_to_remove_idx], pop_desc_contrib, pop_desc_minus_one_contrib)
    gain_add = gain_replace_faster(ind_to_add_desc, pop[ind_to_remove_idx], pop_desc_contrib,
                                   pop_desc_minus_one_contrib, pop_desc_plus_one_contrib)
    computed_delta = -loss_remove + gain_add

    true_entropy_before = entropy(pop, len(pop))
    pop_after = pop.copy()
    pop_after[ind_to_remove_idx] = ind_to_add_desc
    true_entropy_after = entropy(pop_after, len(pop))
    true_delta = true_entropy_after - true_entropy_before

    print("computed loss : " + str(loss_remove))
    print("computed gain : " + str(gain_add))
    print("computed delta : " + str(computed_delta))
    print("true delta : " + str(true_delta))


pop = [[0], [1], [1], [2], [2], [2]]

pop_desc_count = pop_to_desc(pop)
pop_desc_contrib = pop_entropy_per_desc(pop_desc_count, len(pop))
pop_desc_minus_one_contrib = pop_entropy_per_desc(pop_desc_count - 1, len(pop))
pop_desc_plus_one_contrib = pop_entropy_per_desc(pop_desc_count + 1, len(pop))


def test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, idx_to_remove,
                     desc_to_add):
    print("Replace " + str(idx_to_remove) + " -> " + str(desc_to_add))
    test_delta(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, idx_to_remove, desc_to_add)


test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 0, [0])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 2, [2])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 3, [1])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 3, [3])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 2, [0, 1, 2])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 2, [1, 3])
test_replacement(pop, pop_desc_contrib, pop_desc_minus_one_contrib, pop_desc_plus_one_contrib, 2, [3])
