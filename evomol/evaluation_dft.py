import json
import os
import time
from os.path import join
import numpy as np
import cclib
from .evaluation import EvaluationStrategy, EvaluationError
from .molgraphops.molgraph import MolGraph
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import RemoveHs, RemoveStereochemistry


def delete_file(file):
    """
    Code from https://github.com/Cyril-Grl/AlphaSMILES
    Delete the file if exist

    :param file: the file to delete
    :type  file: str
    :return: None
    """
    if os.path.isfile(file):
        os.remove(file)


def remove_files(files_list):
    for filepath in files_list:
        delete_file(filepath)


def smi_to_filename(smi):
    smi = smi.replace('(', '_po_')
    smi = smi.replace(')', '_pf_')
    smi = smi.replace('/', '_s_')
    smi = smi.replace('\\', '_as_')
    smi = smi.replace('@', '_at_')
    smi = smi.replace('#', '_sh_')
    smi = smi.replace("=", '_eq_')

    return smi


def filename_to_smi(filename):
    filename = filename.replace("_po_", "(")
    filename = filename.replace("_pf_", ")")
    filename = filename.replace("_s_", "/")
    filename = filename.replace("_as_", "\\")
    filename = filename.replace("_at_", "@")
    filename = filename.replace("_sh_", "#")
    filename = filename.replace("_eq_", "=")

    return filename


def load_obabel_smi(smi_path, sanitize_mol):
    # Extracting smiles
    with open(smi_path, "r") as f:
        new_smi = f.readline()

        print("obabel new smi : " + new_smi)

        # Loading converged mol
        new_mol = MolFromSmiles(new_smi)

        # Removing stereo information
        RemoveStereochemistry(new_mol)

        # Removing hydrogens
        new_mol = RemoveHs(new_mol)

        # Converting to SMILES
        smi_rdkit = MolGraph(new_mol, sanitize_mol=sanitize_mol).to_aromatic_smiles()
        print("rdkit new smi : " + smi_rdkit)

        return smi_rdkit


def write_input_file(opt_input_path, xyz_path, smi, n_jobs):
    with open(xyz_path, "r") as xyz:
        position = ""
        for i, l in enumerate(xyz):
            if i >= 2:
                position += l

    with open(opt_input_path, "w") as inp:
        inp.write("%Chk=" + smi + "\n")
        inp.write("%NProcShared=" + str(n_jobs) + "\n")
        inp.write("%mem=512MB\n")
        inp.write("#P B3LYP/3-21G* opt gfprint pop=(full,HirshfeldEE)\n")
        inp.write("\n" + smi + "\n\n")
        inp.write("0 1\n")
        inp.write(position + "\n\n\n")


class OPTEvaluationStrategy(EvaluationStrategy):
    """
    Evaluation strategy running a DFT optimization using Gaussian 09 to assess HOMO or LUMO energies.
    The DFT computation is only ran if the SMILES is identical after a molecular mechanics (MM) optimization using
    OpenBabel.
    The DFT computation is considered a success only if the molecule has kept the same SMILES.
    A cache of already performed DFT computations can be provided. It must be a JSON file containing an entry for each
    already computed aromatic canonical SMILES. Each molecule must be represented as a dictionary containing "homo"
    and/or "lumo" keys with the associated value.

    OpenBabel must be installed in a folder referenced with the $OPT_LIBS environment variable. It must
    be set according to the following path. $OPT_LIBS/obabel/openbabel-2.4.1/bin/obabel

    The $OPT_LIBS environment variable must also contain a script named $OPT_LIBS/dft.sh, launching a Gaussian optimization of
    the input file in parameter.
    """

    def __init__(self, prop, n_jobs=2, working_dir_path="/tmp/", cache_files=None):
        """
        Initialization of the DFT evaluation strategy
        :param prop: key of the property to be assessed. Can be "homo" or "lumo"
        :param n_jobs: number of jobs for gaussian optimization
        :param working_dir_path: directory in which computation files will be stored
        :param cache_files: list of JSON file containing a cache of former computations
        """

        self.prop = prop
        self.n_jobs = n_jobs
        self.working_dir_path = working_dir_path
        self.scores = None

        if cache_files is None:
            self.cache_files = []
        else:
            self.cache_files = cache_files

        # Creating the working directory if does not exist
        os.makedirs(working_dir_path, exist_ok=True)

        # Loading cache
        self.cache = {}
        for cache_file in self.cache_files:
            with open(cache_file, "r") as f:
                cache = json.load(f)
                self.cache.update(cache)

        print(str(len(self.cache.keys())) + " molecules in cache")

    def keys(self):
        return [self.prop]

    def get_population_scores(self):
        return self.scores, np.array([self.scores])

    def is_in_cache(self, smi):
        return smi in self.cache

    def get_cache_value(self, smi):
        return self.cache[smi][self.prop]

    def evaluate_individual(self, individual, file_prefix=""):
        """
        Code from https://github.com/Cyril-Grl/AlphaSMILES
        """

        print()
        print()
        print("evaluating " + str(individual.to_aromatic_smiles()))

        # Extracting SMILES
        smi = individual.to_aromatic_smiles()

        # If already calculated, loading from cache
        if self.is_in_cache(smi):

            score = self.get_cache_value(smi)

            if score is None:
                raise EvaluationError("DFT error in cache for " + smi)
            else:
                return score, [score]

        # If score never computed, starting DFT
        else:

            # Converting the smiles to file name compatible
            filename = file_prefix + smi_to_filename(smi)

            # Computing paths
            smi_path = join(self.working_dir_path, filename + ".smi")
            post_MM_smi_path = join(self.working_dir_path, filename + ".MM.smi")
            post_opt_smi_path = join(self.working_dir_path, filename + ".opt.smi")
            xyz_path = join(self.working_dir_path, filename + ".xyz")
            opt_input_path = join(self.working_dir_path, filename + "_OPT.inp")
            opt_log_path = join(self.working_dir_path, filename + "_OPT.log")

            # Writing SMILES in file
            with open(smi_path, "w") as f:
                f.write(smi)

            # Converting SMILES to XYZ after computing MM
            command_obabel = join(os.getenv("OPT_LIBS"), "obabel/openbabel-2.4.1/bin/obabel") + " -ismi " + smi_path \
                             + " -oxyz -O " + xyz_path + " --gen3d"
            os.system(command_obabel)

            # Converting XYZ to smi
            command_obabel = join(os.getenv("OPT_LIBS"), "obabel/openbabel-2.4.1/bin/obabel") + " -ixyz " + xyz_path \
                             + " -osmi -O " + post_MM_smi_path
            os.system(command_obabel)

            post_MM_smi = load_obabel_smi(post_MM_smi_path, sanitize_mol=True)

            if post_MM_smi == smi:

                # Creating input file for OPT
                write_input_file(opt_input_path, xyz_path, smi, self.n_jobs)

                # Calculate OPT in the working directory
                command_opt = "cd " + self.working_dir_path + "; " + join(os.environ["OPT_LIBS"],
                                                                          "dft.sh") + " " + opt_input_path
                print("Starting OPT")
                start = time.time()
                os.system(command_opt)
                stop = time.time()
                print("Execution time OPT: " + repr(int(stop - start)) + "s")

                # Checking that normal termination occurred
                with open(opt_log_path, "r") as log:
                    last_line = log.readlines()[-1]

                # if the OTP end up well
                if "Normal termination" in last_line:

                    # Extracting the smiles from the log file
                    command_obabel = join(os.getenv("OPT_LIBS"),
                                          "obabel/openbabel-2.4.1/bin/obabel") + " -ilog " + opt_log_path \
                                     + " -ocan -O " + post_opt_smi_path
                    os.system(command_obabel)

                    post_opt_smi_rdkit = load_obabel_smi(post_opt_smi_path, sanitize_mol=True)

                    # If before/after SMILES are identical
                    if smi == post_opt_smi_rdkit:

                        with open(opt_log_path, "r") as log:
                            data = cclib.io.ccread(log, optdone_as_list=True)
                            print("There are %i atoms and %i MOs" % (data.natom, data.nmo))
                            homos = data.homos
                            energies = data.moenergies

                        if len(homos) == 1:

                            if self.prop == "homo":
                                score = energies[0][homos[0]]
                            elif self.prop == "lumo":
                                score = energies[0][homos[0] + 1]
                            else:
                                score = None
                            print("score : " + str(score))

                            # Removing files
                            remove_files([smi_path, post_MM_smi_path, post_opt_smi_path, xyz_path, opt_input_path])

                            # Returning score
                            return score, [score]
                        else:
                            # Removing files
                            remove_files([smi_path, post_MM_smi_path, post_opt_smi_path, xyz_path, opt_input_path])

                            raise EvaluationError("DFT error : |homos| > 1 for " + smi)

                    else:
                        # Removing files
                        remove_files([smi_path, post_MM_smi_path, post_opt_smi_path, xyz_path, opt_input_path])

                        raise EvaluationError("DFT error : Different SMILES : " + smi + " " + post_opt_smi_rdkit)

                else:
                    # Removing files
                    remove_files([smi_path, post_MM_smi_path, post_opt_smi_path, xyz_path, opt_input_path])

                    raise EvaluationError("DFT error : Error during OPT for " + smi)
            else:
                # Removing files
                remove_files([smi_path, post_MM_smi_path, post_opt_smi_path, xyz_path, opt_input_path])

                raise EvaluationError("MM error : Different SMILES : " + smi + " " + post_MM_smi)

    def compute_record_scores(self, population):
        self.scores = []
        for idx, ind in enumerate(population):
            if ind is not None:
                self.scores.append(self.evaluate_individual(ind)[0])

    def record_score(self, idx, total_new_score, new_scores):
        if idx == len(self.scores):
            self.scores.append(None)
        self.scores[idx] = total_new_score
