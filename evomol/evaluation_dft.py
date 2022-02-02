import json
import os
import time
import uuid
from os.path import join, exists

import cclib
import numpy as np
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule, UFFOptimizeMolecule
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToXYZBlock, MolToSmiles
from rdkit.Chem.rdmolops import AddHs
from rdkit.Chem.rdmolops import RemoveHs, RemoveStereochemistry

from .evaluation import EvaluationStrategy, EvaluationError


def check_identical_geometries(xyz_path, pre_optim_smiles, post_optim_smiles_path):
    """
    Checking if the given geometry (XYZ filepath) yields the same SMILES as the pre-optimization smiles
    :param xyz_path: path to the XYZ file
    :param pre_optim_smiles: SMILES before geometrical optimization
    :param post_optim_smiles_path: path where to temporary write the post optimization SMILES
    """

    # Converting XYZ to smi
    command_obabel = "obabel -ixyz " + xyz_path + " -osmi -O " + post_optim_smiles_path
    os.system(command_obabel + " > /dev/null 2> /dev/null")

    # Reading post optim SMI file
    post_optim_smi = load_obabel_smi(post_optim_smiles_path)

    # Assigning the value of success if the SMILES has stayed identical
    return pre_optim_smiles == post_optim_smi


def obabel_mmff94_xyz(smiles, **kwargs):
    """
    Returns the string of the XYZ file obtained performing the MMFF94 molecular mechanics optimization of the given
    SMILES using obabel.
    Writing temporary files in $MM_WORKING_DIR if defined or otherwise in /tmp
    :param smiles : input SMILES
    :return : XYZ string of optimized geometry, success (whether the MM optimization was successful and the smiles has
    stayed identical after optimization)
    """

    working_dir = os.environ["MM_WORKING_DIR"] if "MM_WORKING_DIR" in os.environ else "/tmp"

    # Computing RDKIT canonical SMILES
    smi_canon = MolToSmiles(MolFromSmiles(smiles))
    filename_smiles = str(os.getpid()) + "_" + smi_to_filename(smi_canon)

    # Computing files paths
    smi_path = join(working_dir, filename_smiles + ".smi")
    xyz_path = join(working_dir, filename_smiles + ".xyz")
    post_MM_smi_path = join(working_dir, filename_smiles + ".post_MM.smi")

    try:

        # Writing smiles to file
        with open(smi_path, "w") as f:
            f.write(smi_canon)

        # Converting SMILES to XYZ after computing MM (Obabel MMFF94)
        command_obabel = "obabel -ismi " + smi_path + " -oxyz -O " + xyz_path + " --gen3d"
        os.system(command_obabel + " > /dev/null 2> /dev/null")

        # Reading XYZ string
        with open(xyz_path, "r") as f:
            xyz_str = f.read()

        # Success if the post MM smiles is identical the pre MM smiles
        success = check_identical_geometries(xyz_path, smi_canon, post_MM_smi_path)

    except Exception as e:
        success = False
        xyz_str = None
    finally:
        # Removing files
        remove_files([smi_path, xyz_path, post_MM_smi_path])

    return xyz_str, success


def rdkit_mm_xyz(smiles, ff="MMFF94", **kwargs):
    """
    Returns the string of the XYZ file obtained performing the MMFF94 or UFF molecular mechanics optimization of the
    given SMILES using RDKit.
    Writing temporary files in $MM_WORKING_DIR if defined or otherwise in /tmp
    :param smiles: input_SMILES
    :param ff: whether to use MMFF94 force field ("MMFF94") or UFF force field ("UFF")
    :return : XYZ string of optimized geometry, success (whether the MM optimization was successful and the smiles has
    stayed identical after optimization)
    """

    working_dir = os.environ["MM_WORKING_DIR"] if "MM_WORKING_DIR" in os.environ else "/tmp"

    # Converting the molecule to RDKit object
    mol = MolFromSmiles(smiles)
    smi_canon = MolToSmiles(MolFromSmiles(smiles))

    # Setting paths
    filename_smiles = str(os.getpid()) + "_" + smi_to_filename(smi_canon)
    xyz_path = join(working_dir, filename_smiles + '.xyz')
    post_MM_smi_path = join(working_dir, filename_smiles + '.smi')

    # Computing geometry
    try:

        # Adding implicit hydrogens
        mol = AddHs(mol)

        # MM optimization
        EmbedMolecule(mol)

        if ff == "MMFF94":
            value = MMFFOptimizeMolecule(mol, maxIters=kwargs["max_iterations"])
        elif ff == "UFF":
            value = UFFOptimizeMolecule(mol, maxIters=kwargs["max_iterations"])

        # Success if returned value is null
        success_RDKIT_output = value == 0

        # Computing XYZ from optimized molecule
        xyz_str = MolToXYZBlock(mol)

        # Writing optimized XYZ to file
        with open(xyz_path, "w") as f:
            f.writelines(xyz_str)

        # Success if the optimization has converged and the post MM smiles is identical the pre MM smiles
        success = success_RDKIT_output and check_identical_geometries(xyz_path, smi_canon, post_MM_smi_path)

    except Exception as e:
        success = False
        xyz_str = None
    finally:
        # Removing files
        remove_files([post_MM_smi_path, xyz_path])

    return xyz_str, success


def rdkit_mmff94_xyz(smiles, **kwargs):
    """
    Returns the string of the XYZ file obtained performing the MMFF94 molecular mechanics optimization of the given
    SMILES using RDKit.
    Writing temporary files in $MM_WORKING_DIR if defined or otherwise in /tmp
    :param smiles: input_SMILES
    :param max_iterations: max number of iterations (default 500)
    :return : XYZ string of optimized geometry, success (whether the MM optimization was successful and the smiles has
    stayed identical after optimization)


    NOTE : DEPRECATED FUNCTION. Kept here for backwards compatibility. Now it is better to call rdkit_mm_xyz using
    the ff="MMFF94" parameter.
    """

    working_dir = os.environ["MM_WORKING_DIR"] if "MM_WORKING_DIR" in os.environ else "/tmp"

    # Converting the molecule to RDKit object
    mol = MolFromSmiles(smiles)
    smi_canon = MolToSmiles(MolFromSmiles(smiles))

    # Setting paths
    filename_smiles = str(os.getpid()) + "_" + smi_to_filename(smi_canon)
    xyz_path = join(working_dir, filename_smiles + '.xyz')
    post_MM_smi_path = join(working_dir, filename_smiles + '.smi')

    # Computing geometry
    try:

        # Adding implicit hydrogens
        mol = AddHs(mol)

        # MM optimization
        EmbedMolecule(mol)

        value = MMFFOptimizeMolecule(mol, maxIters=kwargs["max_iterations"])

        # Success if returned value is null
        success_RDKIT_output = value == 0

        # Computing XYZ from optimized molecule
        xyz_str = MolToXYZBlock(mol)

        # Writing optimized XYZ to file
        with open(xyz_path, "w") as f:
            f.writelines(xyz_str)

        # Success if the optimization has converged and the post MM smiles is identical the pre MM smiles
        success = success_RDKIT_output and check_identical_geometries(xyz_path, smi_canon, post_MM_smi_path)

    except Exception as e:
        success = False
        xyz_str = None
    finally:
        # Removing files
        remove_files([post_MM_smi_path, xyz_path])

    return xyz_str, success


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


def load_obabel_smi(smi_path):
    """
    Converting a OpenBabel SMILES into a canonical aromatic RDKit SMILES
    :param smi_path:
    :return:
    """

    # Extracting smiles
    with open(smi_path, "r") as f:
        new_smi = f.readline()

        # Loading converged mol
        new_mol = MolFromSmiles(new_smi)

        # Removing stereo information
        RemoveStereochemistry(new_mol)

        # Removing hydrogens
        new_mol = RemoveHs(new_mol)

        # Converting to SMILES
        smi_rdkit = MolToSmiles(MolFromSmiles(MolToSmiles(new_mol)))

        return smi_rdkit


def write_input_file(opt_input_path, xyz_path, smi, n_jobs, dft_base="3-21G*"):
    with open(xyz_path, "r") as xyz:
        position = ""
        for i, l in enumerate(xyz):
            if i >= 2:
                position += l

    with open(opt_input_path, "w") as inp:
        inp.write("%Chk=" + smi + "\n")
        inp.write("%NProcShared=" + str(n_jobs) + "\n")
        inp.write("%mem=512MB\n")
        inp.write("#P B3LYP/" + dft_base + " opt Symmetry=(NoInt,NoGrad,None) gfprint pop=(full,HirshfeldEE)\n")
        inp.write("\n" + smi + "\n\n")
        inp.write("0 1\n")
        inp.write(position + "\n\n\n")


class SharedLastComputation:
    """
    Object that can be shared by several OPTEvaluationStrategy instances and that contains the values of the last
    DFT computation. It allows to only perform one calculation in case of the evaluation of a combination of
    OPTEvaluationStrategy instances.
    """

    def __init__(self):
        self.smiles = None
        self.homo = None
        self.lumo = None
        self.gap = None
        self.homo_m1 = None


def compress_log_file(log_path):
    """
    Using Gzip to compress the log output file
    :param log_path: path to the log file
    :return:
    """

    cmd = "gzip -f " + log_path
    os.system(cmd)


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

    The $OPT_LIBS environment variable must also contain a script named $OPT_LIBS/dft.sh, starting a Gaussian
    optimization of the input file in parameter.
    """

    def __init__(self, prop, n_jobs=1, working_dir_path="/tmp/", cache_files=None, MM_program="obabel_mmff94",
                 cache_behaviour="retrieve_OPT_data", remove_chk_file=True, shared_last_computation=None,
                 dft_base="3-21G*"):
        """
        Initialization of the DFT evaluation strategy
        :param prop: key of the property to be assessed. Can be "homo", "lumo", "gap" or "homo-1"
        :param n_jobs: number of jobs for gaussian optimization
        :param working_dir_path: directory in which computation files will be stored
        :param cache_files: list of JSON file containing a cache of former computations
        :param MM_program: program used to compute MM. Options are :
            - "obabel" or "obabel_mmff94" for MMFF94 optimization using OpenBabel
            - "rdkit" or "rdkit_mmff94" for MMFF94 optimization using RDKit
            - "rdkit_uff" for UFF optimization using RDKit
        :param cache_behaviour : configuration of the behaviour when cache files are given. "retrieve_OPT_data"
        (default): if the molecule is known in the cache, no DFT computation is made and values are retrieved.
        "compute_again_delete_files": DFT computation are made for all molecules but DFT files are removed for molecules
        that are already in cache.
        :param remove_chk_file: whether the G09 CHK file is removed after DFT computation (default:True)
        :param shared_last_computation: SharedLastComputation instance to share the values of the last computation
        values with several OPTEvaluationStrategy instances
        :param dft_base: base of G09 DFT computation (default : "3-21G*")
        """

        super().__init__()
        self.prop = prop
        self.n_jobs = n_jobs
        self.scores = None
        self.shared_last_computation = shared_last_computation

        self.MM_program = MM_program

        if cache_files is None:
            self.cache_files = []
        else:
            self.cache_files = cache_files

        self.cache = {}
        # Reversing the list so that values of first files are taken primarily if there exists an intersection of
        # SMILES keys
        self.cache_files.reverse()

        # Loading cache
        for cache_file in self.cache_files:
            with open(cache_file, "r") as f:
                cache = json.load(f)
                self.cache.update(cache)

        # Creating the root directory if does not exist
        os.makedirs(working_dir_path, exist_ok=True)

        # Computing a unique identifier for the current instance
        computed_uuid = str(uuid.uuid4())

        # Computing the working directory path by concatenating the working directory root and the uuid
        self.working_dir_path_uuid = join(working_dir_path, computed_uuid)

        self.cache_behaviour = cache_behaviour
        self.remove_chk_file = remove_chk_file
        self.dft_base = dft_base

        print("DFT MM " + str(self.MM_program))
        print(str(len(self.cache.keys())) + " molecules in cache")

    def keys(self):
        return [self.prop]

    def is_in_cache(self, smi):
        return smi in self.cache

    def get_cache_value(self, prop, smi):

        # Computing the success value if the "success" key in in the entry
        if "success" in self.cache[smi]:
            success = self.cache[smi]["success"]
        # Otherwise the success is computed as whether the property is not None
        else:
            success = prop in self.cache[smi] and self.cache[smi][prop] is not None

        if prop in self.cache[smi]:
            value = self.cache[smi][prop]
        else:
            value = None
            success = False

        return value, success

    def remove_evaluation_files(self, post_opt_smi_path, xyz_path, opt_input_path, chk_path,
                                log_path, is_in_cache):
        """
        Removing files created during the MM + DFT computation.
        The existence of the files is checked before removal.
        CHK file is removed iff. self.remove_chk_file is True
        Log file is removed iff. molecule is in cache and self.cache_behaviour is set to "compute_again_delete_files"
        :param post_opt_smi_path: path to the file containing the SMILES after DFT optimization (.smi)
        :param xyz_path: path to the file containing the XYZ data after MM optimization (.xyz)
        :param opt_input_path: path to the input of G09 (.inp)
        :param chk_path: path to the CHK file generated by G09 (.chk)
        :param log_path: path to the G09 LOG path (.log)
        :param is_in_cache: whether the molecule is known in the cache
        """

        remove_files([post_opt_smi_path, xyz_path, opt_input_path])

        # Removing CHK file if self.remove_chk_file is set to True
        if self.remove_chk_file:
            remove_files([chk_path])

        # Removing log path if solution is in cache and self.cache_behaviour is set to "compute_again_delete_files"
        if self.cache_behaviour == "compute_again_delete_files" and is_in_cache:
            remove_files([log_path])

    def evaluate_individual(self, individual, to_replace_idx=None, file_prefix=""):
        """
        Code from https://github.com/Cyril-Grl/AlphaSMILES (Cyril Grelier)

        MIT License

        Copyright (c) 2019 Cyril-Grl

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.
        """

        super().evaluate_individual(individual, to_replace_idx)

        # Extracting SMILES
        smi = individual.to_aromatic_smiles()

        # Converting the smiles to file name compatible
        filename = file_prefix + smi_to_filename(smi)

        # Computing paths
        post_opt_smi_path = join(self.working_dir_path_uuid, filename + ".opt.smi")
        xyz_path = join(self.working_dir_path_uuid, filename + ".xyz")
        opt_input_path = join(self.working_dir_path_uuid, filename + "_OPT.inp")
        opt_log_path = join(self.working_dir_path_uuid, filename + "_OPT.log")
        chk_path = join(self.working_dir_path_uuid, smi + ".chk")

        # Computing whether the solution is known in the cache
        ind_is_in_cache = self.is_in_cache(smi)

        # If already calculated and cache behaviour is set to "retrieve_OPT_data", loading from cache
        if ind_is_in_cache and self.cache_behaviour == "retrieve_OPT_data":

            homo_cache, success_homo_cache = self.get_cache_value("homo", smi)
            homo_m1_cache, success_homo_m1_cache = self.get_cache_value("homo-1", smi)
            lumo_cache, success_lumo_cache = self.get_cache_value("lumo", smi)
            gap_cache, success_gap_cache = self.get_cache_value("gap", smi)

            if self.prop == "gap":
                success = success_homo_cache and success_lumo_cache and success_gap_cache
                score = gap_cache
                scores = [gap_cache]
            elif self.prop == "homo":
                success = success_homo_cache
                score = homo_cache
                scores = [homo_cache]
            elif self.prop == "homo-1":
                success = success_homo_m1_cache
                score = homo_m1_cache
                scores = [homo_m1_cache]
            elif self.prop == "lumo":
                success = success_lumo_cache
                score = lumo_cache
                scores = [lumo_cache]

            if not success:
                raise EvaluationError("DFT failure in cache for " + smi)
            else:
                return score, scores

        # Case in which the computation has just been performed by another OPTEvaluationStrategy instance that
        # shares the same SharedLastComputation instance
        elif self.shared_last_computation is not None and individual.to_aromatic_smiles() == self.shared_last_computation.smiles:

            # Returning score
            if self.prop == "homo":
                return self.shared_last_computation.homo, [self.shared_last_computation.homo]
            elif self.prop == "lumo":
                return self.shared_last_computation.lumo, [self.shared_last_computation.lumo]
            elif self.prop == "gap":
                return self.shared_last_computation.gap, [self.shared_last_computation.gap]
            elif self.prop == "homo-1":
                return self.shared_last_computation.homo_m1, [self.shared_last_computation.homo_m1]

        # If score never computed or cache behaviour is set to "compute_again_delete_files", starting DFT
        else:

            # Creating the working directory if does not exist
            os.makedirs(self.working_dir_path_uuid, exist_ok=True)

            print("computing dft for " + str((individual.to_aromatic_smiles())))

            try:

                # Performing Obabel MM
                if self.MM_program == "obabel" or self.MM_program == "obabel_mmff94":
                    # Converting SMILES to XYZ after computing MM (RDKit MMFF94)
                    xyz_str, success_MM = obabel_mmff94_xyz(smi)

                # Performing RDkit MM
                elif self.MM_program == "rdkit" or self.MM_program == "rdkit_mmff94":
                    # Converting SMILES to XYZ after computing MM (RDKit MMFF94)
                    xyz_str, success_MM = rdkit_mm_xyz(smi, ff="MMFF94", max_iterations=500)

                elif self.MM_program == "rdkit_uff":
                    # Converting SMILES to XYZ after computing MM (RDKit MMFF94)
                    xyz_str, success_MM = rdkit_mm_xyz(smi, ff="UFF", max_iterations=500)

                if success_MM:

                    # Writing optimized XYZ to file
                    with open(xyz_path, "w") as f:
                        f.writelines(xyz_str)

                    # Creating input file for OPT
                    write_input_file(opt_input_path, xyz_path, smi, self.n_jobs, dft_base=self.dft_base)

                    # Calculate OPT in the working directory
                    command_opt = "cd " + self.working_dir_path_uuid + "; " + join(os.environ["OPT_LIBS"],
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
                        command_obabel = "obabel -ilog " + opt_log_path + " -ocan -O " + post_opt_smi_path
                        os.system(command_obabel)

                        post_opt_smi_rdkit = load_obabel_smi(post_opt_smi_path)

                        # If before/after SMILES are identical
                        if smi == post_opt_smi_rdkit:

                            with open(opt_log_path, "r") as log:
                                data = cclib.io.ccread(log, optdone_as_list=True)
                                print("There are %i atoms and %i MOs" % (data.natom, data.nmo))
                                homos = data.homos
                                energies = data.moenergies

                            if len(homos) == 1:

                                homo = energies[0][homos[0]]
                                lumo = energies[0][homos[0] + 1]
                                homo_m1 = energies[0][homos[0] - 1]
                                gap = abs(homo - lumo)

                                # Removing files
                                self.remove_evaluation_files(post_opt_smi_path, xyz_path, opt_input_path, chk_path,
                                                             opt_log_path, is_in_cache=ind_is_in_cache)

                                # Compressing log file
                                compress_log_file(opt_log_path)

                                # Saving values in SharedLastComputation instance if defined
                                if self.shared_last_computation is not None:
                                    self.shared_last_computation.smiles = individual.to_aromatic_smiles()
                                    self.shared_last_computation.homo = homo
                                    self.shared_last_computation.lumo = lumo
                                    self.shared_last_computation.gap = gap
                                    self.shared_last_computation.homo_m1 = homo_m1

                                # Returning score
                                if self.prop == "homo":
                                    return homo, [homo]
                                elif self.prop == "lumo":
                                    return lumo, [lumo]
                                elif self.prop == "gap":
                                    return gap, [gap]
                                elif self.prop == "homo-1":
                                    return homo_m1, [homo_m1]

                            else:
                                raise EvaluationError("DFT error : |homos| > 1 for " + smi)
                        else:
                            raise EvaluationError(
                                "DFT error : Different SMILES : " + smi + " " + post_opt_smi_rdkit)

                    else:
                        raise EvaluationError("DFT error : Error during OPT for " + smi)
                else:
                    raise EvaluationError("MM error")

            except Exception as e:
                print(e)
                # Removing files
                self.remove_evaluation_files(post_opt_smi_path, xyz_path, opt_input_path, chk_path, opt_log_path,
                                             is_in_cache=ind_is_in_cache)

                if exists(opt_log_path):
                    compress_log_file(opt_log_path)

                raise EvaluationError("DFT caused exception " + str(e))

    def compute_record_scores_init_pop(self, population):
        self.scores = []

        for idx, ind in enumerate(population):
            if ind is not None:
                score, _ = self.evaluate_individual(ind)
                self.scores.append(score)

    def record_ind_score(self, idx, new_total_score, new_scores, new_individual):

        if idx == len(self.scores):
            self.scores.append(None)

        self.scores[idx] = new_total_score

    def get_population_scores(self):

        return self.scores, np.array([self.scores])