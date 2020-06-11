import csv
from os import makedirs

import networkx as nx
from PIL import Image, ImageChops, ImageDraw
from evomol.evaluation import EvaluationStrategy
from matplotlib.colors import LogNorm
from evomol.molgraphops.molgraph import MolGraph
from networkx.drawing.nx_agraph import graphviz_layout
from rdkit.Chem.Draw import MolsToGridImage, MolToImage, DrawingOptions
from rdkit.Chem.rdmolfiles import MolFromSmiles
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

figsize = (15, 10)
dpi = 500


def extract_history_data(filepath):
    actions_history_scores = {}
    actions_history_smi = {}

    with open(filepath, "r") as f:
        reader = csv.reader(f)
        keys = []
        for i, row in enumerate(reader):
            if i == 0:
                for key in row:
                    keys.append(key)
            else:
                scores_dict = {}
                for j, value in enumerate(row):
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


def trim(im):
    """
    from https://stackoverflow.com/questions/14211340/automatically-cropping-an-image-with-python-pil
    """
    bg = Image.new(im.mode, im.size, im.getpixel((0, 0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 1, 0)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)


def crop_image_with_transparency(img):
    # Extracting zone to conserve (https://stackoverflow.com/questions/14211340/automatically-cropping-an-image-with-python-pil)
    bg = Image.new(img.mode, img.size, img.getpixel((0, 0)))
    diff = ImageChops.difference(img, bg)
    diff = ImageChops.add(diff, diff, 1, 0)
    bbox = diff.getbbox()

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


def crop_images(images_dict):
    """
    Cropping the images on the largest occupied box
    """

    bboxes = {}

    # Extracting boxes for all images
    for k, img in images_dict.items():
        bg = Image.new(img.mode, img.size, img.getpixel((0, 0)))
        diff = ImageChops.difference(img, bg)
        diff = ImageChops.add(diff, diff, 1, 0)
        bboxes[k] = diff.getbbox()

    # Computing largest box
    l = min([v[0] for v in bboxes.values()])
    u = min([v[1] for v in bboxes.values()])
    r = max([v[1] for v in bboxes.values()])
    b = max([v[1] for v in bboxes.values()])
    bbox = (l, u, r, b)

    # Cropping images
    for k, img in images_dict.items():
        bboxes[k] = img.crop(bbox)


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

    img = MolsToGridImage(mols, legends=legends, molsPerRow=mols_per_row, subImgSize=(200, 200), maxMols=200)
    return img


def best_score_node(actions_history_scores_pop, actions_history_scores_removed):
    # Extracting best score node
    best_score = - float("inf")
    best_score_key = None
    for key, score_dict in actions_history_scores_pop.items():
        if score_dict["total"] > best_score:
            best_score_key = key
            best_score = score_dict["total"]
    for key, score_dict in actions_history_scores_removed.items():
        if score_dict["total"] > best_score:
            best_score_key = key
            best_score = score_dict["total"]

    return best_score_key


def normalize_layout(input_layout):
    x_values = [x for key, (x, y) in input_layout.items()]
    y_values = [y for key, (x, y) in input_layout.items()]

    normalized_layout = {}

    for i, key in enumerate(input_layout.keys()):
        normalized_layout[key] = ((x_values[i] - min(x_values)) / (max(x_values) - min(x_values)),
                                  (y_values[i] - min(y_values)) / (max(y_values) - min(y_values)))

    return normalized_layout


def exploration_graph(pop_filepath, removed_hist_filepath, neighbours_threshold, root_node="C", plot_images=False,
                      mol_size=0.1, figsize=(40, 30), draw_scores=False, draw_actions=False, plot_labels=False,
                      layout="dot", cmap="inferno", output_files_prefix=None, dpi=300, legend_offset=(0, 0),
                      legend_scores_keys_strat=None, problem_type="max",
                      mols_per_row=4, draw_n_mols=None):
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

    cc = list(nx.connected_components(graph))
    print(str(len(cc)) + " connected components")

    cmap = plt.get_cmap(cmap)

    # Setting color to nodes depending on if they are in the final population
    colors = []
    sizes = []
    labels = {}
    n_ind_in_pop = 0
    next_label_to_assign = 1
    best_score_key = best_score_node(actions_history_scores_pop, actions_history_scores_removed)
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
            colors.append(cmap(actions_history_scores_pop[node]["total"]))

        elif node in actions_history_scores_removed:
            colors.append(cmap(actions_history_scores_removed[node]["total"]))

    print(len(colors))
    print(str(n_ind_in_pop) + " in final population")

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
                     font_size=25)

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
                    fig.text(xa + legend_offset[0], ya + legend_offset[1], scores[n], fontsize=15,
                             bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'))

    if plot_images:
        img_labels = draw_mol_labels(labels, actions_history_smi_pop, actions_history_smi_removed,
                                     actions_history_scores_pop, actions_history_scores_removed,
                                     legend_scores_keys_strat=legend_scores_keys_strat, problem_type=problem_type,
                                     mols_per_row=mols_per_row, draw_n_mols=draw_n_mols)
    else:
        img_labels = None

    if not plot_images:
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm)

    if output_files_prefix is not None:

        makedirs(output_files_prefix, exist_ok=True)
        plt.savefig(output_files_prefix + "_expl_tree.png", dpi=dpi)

    plt.show()

    return plt, img_labels
