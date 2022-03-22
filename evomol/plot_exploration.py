import csv
from os.path import join

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from PIL import Image, ImageDraw
from matplotlib.colors import LogNorm
from networkx.drawing.nx_agraph import graphviz_layout
from rdkit.Chem.Draw import MolsToGridImage, MolToImage, DrawingOptions
from rdkit.Chem.rdmolfiles import MolFromSmiles

from .evaluation import EvaluationStrategy
from .molgraphops.molgraph import MolGraph

figsize = (15, 10)
dpi = 500


def extract_history_data(filepath):
    actions_history_scores = {}
    actions_history_smi = {}

    with open(filepath, "r", newline='') as f:
        reader = csv.reader(f)
        keys = []
        for i, row in enumerate(reader):
            if i == 0:
                for key in row:
                    keys.append(key)
            else:
                scores_dict = {}
                curr_row = list(enumerate(row))

                # Checking that current row is defined (= there is a non-empty action history value)
                if curr_row[keys.index("history_data")][1] != "":
                    for j, value in curr_row:
                        if keys[j] == "history_data":
                            action_history = value
                        elif keys[j] == "smiles":
                            smi = value
                        else:
                            scores_dict[keys[j]] = float(value)

                    actions_history_scores[action_history] = scores_dict
                    actions_history_smi[action_history] = smi

    return actions_history_scores, actions_history_smi


def fix_act_desc(act_desc):
    act_desc = act_desc.replace("MoveFG", "Mv")
    act_desc = act_desc.replace("AddA", "Ad")
    act_desc = act_desc.replace("RmA", "Rm")
    act_desc = act_desc.replace("InsA", "In")
    act_desc = act_desc.replace("CutA", "Ct")
    act_desc = act_desc.replace("ChB", "Ch")

    return act_desc


def compute_edge_label(act_at_step):
    edge_label_list = []

    act_at_step = fix_act_desc(act_at_step)

    for act in act_at_step.split("_"):
        edge_label_list.append(act.split("-")[0])

    return "_".join(edge_label_list)


def build_graph(graph, actions_history_list, edge_labels):
    # Iterating over all molecules
    for action_history in actions_history_list:

        # Iterating over all children of the molecule
        actions_list = action_history.split("|")
        for i in range(1, len(actions_list)):
            # Extracting the id of the molecule and the id of its child
            node = "|".join(actions_list[:i])
            child_node = "|".join(actions_list[:i + 1])

            # Recording the edge label
            edge_labels[(node, child_node)] = compute_edge_label(actions_list[i])

            # Adding the molecule to the graph
            graph.add_node(node)

            # Adding the child to the graph and linking it to the parent
            graph.add_node(child_node)
            graph.add_edge(node, child_node)


def crop_image_with_transparency(img):
    """
    Cropping image with a transparent channel.
    :param img:
    :return:
    """

    # Insuring the image has an alpha channel
    img.putalpha(255)

    # Image to numpy array
    image_data = np.array(img)

    # Computing the mask of white pixels
    r, g, b, a = np.rollaxis(image_data, axis=-1)
    white_pixels_mask = np.logical_and(np.logical_and(r == 255, g == 255), b == 255)

    # Replacing all white pixels by transparent pixels
    a[white_pixels_mask] = 0

    # Computing bounding box of non zero pixels
    bbox = Image.fromarray(image_data).getbbox()

    # If the bbox is None then the image is empty then a single pixel is selected
    if bbox is None:
        l, u, r, b = 0, 0, 0, 0
    else:
        l, u, r, b = bbox

    w, h = img.size

    mask = Image.new('L', img.size, color=255)
    epsilon = 10

    # Applying transparency (https://stackoverflow.com/questions/4379978/python-pil-how-to-make-area-transparent-in-png)
    for transparent_zone in [(0, 0, l - epsilon, h), (0, 0, w, u - epsilon), (r + epsilon, 0, w, h),
                             (0, b + epsilon, w, h)]:
        draw = ImageDraw.Draw(mask)
        draw.rectangle(transparent_zone, fill=0)
        img.putalpha(mask)

    return img


def compute_mol_legend(action_history_k, smi, action_history_scores, legend_scores_keys_strat=None):
    legend = ""
    last = 0
    scores_float = []

    if legend_scores_keys_strat is not None:
        for i, key_strat in enumerate(legend_scores_keys_strat):

            score = None

            if isinstance(key_strat, str):
                score = action_history_scores[action_history_k][key_strat]
            elif isinstance(key_strat, EvaluationStrategy):
                score = key_strat.evaluate_individual(MolGraph(MolFromSmiles(smi), sanitize_mol=True))

            scores_float.append(score)
            score_str = "{:.2f}".format(score)

            if i == 1:
                legend += " ["
            elif i > 1:
                legend += ", "

            legend += score_str
            last = i

        if last >= 1:
            legend += "]"

    return legend, scores_float


def compute_mol_attributes(graph, labels_dict, actions_history_smi_pop, actions_history_smi_removed,
                           actions_history_scores_pop, actions_history_scores_removed, legend_scores_keys_strat=None):
    images_attributes = {}
    scores_attributes = {}

    draw_opt = DrawingOptions()
    draw_opt.coordScale = 0.9
    draw_opt.dotsPerAngstrom = 30

    for action_history_k in labels_dict.keys():

        if action_history_k in actions_history_smi_pop:

            smi = actions_history_smi_pop[action_history_k]
            img = MolToImage(MolFromSmiles(smi), size=(800, 800),
                             options=draw_opt)
            images_attributes[action_history_k] = crop_image_with_transparency(img)

            legend, _ = compute_mol_legend(action_history_k, smi, actions_history_scores_pop,
                                           legend_scores_keys_strat)
            scores_attributes[action_history_k] = legend

        else:

            smi = actions_history_smi_removed[action_history_k]
            img = MolToImage(MolFromSmiles(smi), size=(800, 800),
                             options=draw_opt)
            images_attributes[action_history_k] = crop_image_with_transparency(img)

            legend, _ = compute_mol_legend(action_history_k, smi, actions_history_scores_removed,
                                           legend_scores_keys_strat)

            scores_attributes[action_history_k] = legend

    nx.set_node_attributes(graph, images_attributes, "image")
    nx.set_node_attributes(graph, scores_attributes, "score_label")


def draw_mol_labels(labels_dict, actions_history_smi_pop, actions_history_smi_removed,
                    actions_history_scores_pop, actions_history_scores_removed, legend_scores_keys_strat=None,
                    problem_type="max", mols_per_row=4, draw_n_mols=None):
    smi_to_draw = {}
    legends_to_draw = {}
    scores_float = {}

    for action_history_k in labels_dict.keys():

        if labels_dict[action_history_k] != "":

            if action_history_k in actions_history_smi_pop:
                smi = actions_history_smi_pop[action_history_k]
                smi_to_draw[labels_dict[action_history_k]] = smi

                legend, scores = compute_mol_legend(action_history_k, smi, actions_history_scores_pop,
                                                    legend_scores_keys_strat)
                legends_to_draw[labels_dict[action_history_k]] = legend
                scores_float[labels_dict[action_history_k]] = scores
            else:
                smi = actions_history_smi_removed[action_history_k]
                smi_to_draw[labels_dict[action_history_k]] = smi

                legend, scores = compute_mol_legend(action_history_k, smi, actions_history_scores_removed,
                                                    legend_scores_keys_strat)
                legends_to_draw[labels_dict[action_history_k]] = legend
                scores_float[labels_dict[action_history_k]] = scores

    mols = []
    legends = []
    scores_to_sort = []
    for k, smi in smi_to_draw.items():
        mols.append(MolFromSmiles(smi))
        legends.append(legends_to_draw[k])
        scores_to_sort.append(scores_float[k][0])


    mols = np.array(mols)
    legends = np.array(legends)

    # Sorting molecules
    sorted_order = np.argsort(scores_to_sort)
    if problem_type == "max":
        sorted_order = sorted_order[::-1]

    # Filtering molecules if necessary
    if draw_n_mols is not None:
        mols = mols[:draw_n_mols]
        legends = legends[:draw_n_mols]

    legends = list(legends[sorted_order])
    mols = list(mols[sorted_order])

    img = MolsToGridImage(mols, legends=legends, molsPerRow=mols_per_row, subImgSize=(200, 200))
    return img


def best_score_node(actions_history_scores_pop, actions_history_scores_removed, prop_to_study_key):
    # Extracting best score node
    best_score = - float("inf")
    best_score_key = None
    for key, score_dict in actions_history_scores_pop.items():
        if score_dict[prop_to_study_key] > best_score:
            best_score_key = key
            best_score = score_dict[prop_to_study_key]
    for key, score_dict in actions_history_scores_removed.items():
        if score_dict[prop_to_study_key] > best_score:
            best_score_key = key
            best_score = score_dict[prop_to_study_key]

    return best_score_key


def normalize_layout(input_layout):
    x_values = [x for key, (x, y) in input_layout.items()]
    y_values = [y for key, (x, y) in input_layout.items()]

    normalized_layout = {}

    for i, key in enumerate(input_layout.keys()):
        normalized_layout[key] = ((x_values[i] - min(x_values)) / (max(x_values) - min(x_values)),
                                  (y_values[i] - min(y_values)) / (max(y_values) - min(y_values)))

    return normalized_layout


def exploration_graph(model_path, neighbours_threshold=0, root_node="C", plot_images=False,
                      mol_size=0.1, figsize=(15, 10), draw_scores=False, draw_actions=False, plot_labels=False,
                      layout="dot", cmap="inferno", prop_to_study_key="total", dpi=300, legend_offset=(0, 0),
                      legend_scores_keys_strat=None, problem_type="max",
                      mols_per_row=4, draw_n_mols=None, legends_font_size=15):

    # Computing file names file names
    pop_filepath = join(model_path, "pop.csv")
    removed_hist_filepath = join(model_path, "removed_ind_act_history.csv")

    # Extracting actions history of individuals in population
    actions_history_scores_pop, actions_history_smi_pop = extract_history_data(pop_filepath)

    # Extracting actions history of removed individuals
    actions_history_scores_removed, actions_history_smi_removed = extract_history_data(removed_hist_filepath)

    # Initialization of the graph
    graph = nx.Graph()
    edge_labels = {}

    # Building the graph for removed individuals
    build_graph(graph, list(actions_history_scores_removed.keys()), edge_labels)

    # Adding the individuals in population to the graph
    build_graph(graph, list(actions_history_scores_pop.keys()), edge_labels)

    cmap = plt.get_cmap(cmap)

    # Setting color to nodes depending on if they are in the final population
    colors = []
    sizes = []
    labels = {}
    next_label_to_assign = 1
    best_score_key = best_score_node(actions_history_scores_pop, actions_history_scores_removed, prop_to_study_key)
    for node in graph.nodes:

        if plot_images:
            sizes.append(0)
        else:
            if node == root_node:
                sizes.append(40)
            elif node in actions_history_scores_pop:
                sizes.append(1)
            else:
                sizes.append(1)

        curr_node_neighbours = graph.neighbors(node)
        if len(list(curr_node_neighbours)) >= neighbours_threshold or node == best_score_key or node == root_node \
                or node in actions_history_scores_pop:
            labels[node] = str(next_label_to_assign)
            next_label_to_assign += 1
        else:
            labels[node] = ""

        if node in actions_history_scores_pop:
            colors.append(cmap(actions_history_scores_pop[node][prop_to_study_key]))

        elif node in actions_history_scores_removed:
            colors.append(cmap(actions_history_scores_removed[node][prop_to_study_key]))

    # Drawing the graph
    plt.figure(figsize=figsize)

    layout = graphviz_layout(graph, prog=layout, root=root_node, args="-Gnodesep=700 -Gminlen=100")
    print("layout computed")
    plt.clf()

    if plot_labels:
        labels_to_print = labels
    else:
        labels_to_print = {}

    nx.draw_networkx(graph, pos=layout, node_size=sizes, node_color=colors, with_labels=True, labels=labels_to_print,
                     font_size=legends_font_size)

    plt.axis('off')

    if draw_actions:
        nx.draw_networkx_edge_labels(graph, pos=layout, edge_labels=edge_labels,
                                     bbox=dict(facecolor='#ffffff', edgecolor='none', alpha=0.3), font_size=15)

    ax = plt.gca()
    # ax.axis('off')
    fig = plt.gcf()
    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform

    if plot_images or draw_scores:

        compute_mol_attributes(graph, labels, actions_history_smi_pop, actions_history_smi_removed,
                               actions_history_scores_pop, actions_history_scores_removed, legend_scores_keys_strat)

        images = nx.get_node_attributes(graph, "image")
        scores = nx.get_node_attributes(graph, "score_label")

        for n in graph.nodes():

            if labels[n] != "":

                (x, y) = layout[n]
                xx, yy = trans((x, y))  # figure coordinates
                xa, ya = trans2((xx, yy))  # axes coordinates

                if plot_images:
                    a = plt.axes([xa - mol_size / 2, ya - mol_size / 2, mol_size, mol_size])
                    a.imshow(images[n])
                    a.set_aspect('equal')
                    a.patch.set_alpha(0)
                    a.axis('off')

                if draw_scores:
                    fig.text(xa + legend_offset[0], ya + legend_offset[1], scores[n], fontsize=legends_font_size,
                             bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'))

    if plot_images:
        img_labels = draw_mol_labels(labels, actions_history_smi_pop, actions_history_smi_removed,
                                     actions_history_scores_pop, actions_history_scores_removed,
                                     legend_scores_keys_strat=legend_scores_keys_strat, problem_type=problem_type,
                                     mols_per_row=mols_per_row, draw_n_mols=draw_n_mols)
        with open(join(model_path, "mol_table.png"), "wb") as f:
            img_labels.save(f, "png")

    if not plot_images:
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm)

    plt.savefig(join(model_path, "expl_tree.png"), dpi=dpi)
