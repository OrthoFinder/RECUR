from __future__ import absolute_import
from recur import  __version__, helpinfo
import subprocess
import multiprocessing as mp               
import platform                                
import sys
import os
import csv
import shutil
from typing import Optional, Dict, List, Tuple, Union, Set
import threading       
import gc
import time
import datetime
import numpy as np
from collections import Counter
import scipy as spy
import warnings
import psutil
import signal
import traceback
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
import dendropy
from functools import partial
from collections import defaultdict
from recur.run import run_commands
from recur.utils import files, util, process_args
import logging
import warnings
# warnings.filterwarnings("ignore", module='dendropy')

def setup_environment():
    max_int = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int = int(max_int / 10)
    
    # Set maximum recursion limit and OpenBLAS thread limit
    sys.setrecursionlimit(10**6)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    my_env = os.environ.copy()

    if getattr(sys, 'frozen', False):  # For PyInstaller
        base_dir = os.path.dirname(sys.executable)
    else:
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

    local_iqtree2_path = os.path.join(base_dir, 'recur', 'bin', 'iqtree2')
    bin_dir = os.path.dirname(local_iqtree2_path)

    conda_prefix = my_env.get('CONDA_PREFIX')
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, 'Scripts') if os.name == 'nt' else os.path.join(conda_prefix, 'bin')
        my_env['PATH'] = conda_bin + os.pathsep + my_env['PATH']
    
    # Restore original paths if running from a frozen environment
    if getattr(sys, 'frozen', False):
        if os.name == 'nt':  # Windows-specific
            my_env['PATH'] = my_env.get('PATH_ORIG', '')
        else:  # Unix-like systems
            my_env['LD_LIBRARY_PATH'] = my_env.get('LD_LIBRARY_PATH_ORIG', '')
            my_env['DYLD_LIBRARY_PATH'] = my_env.get('DYLD_LIBRARY_PATH_ORIG', '')
    
    return my_env, local_iqtree2_path, bin_dir, conda_prefix 

def CanRunCommand(command: str) -> bool:
    try:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        stdout, stderr = process.communicate()

        # if stdout:
        #     print(stdout.strip())
        # if stderr:
        #     print(stderr.strip())

        returncode = process.returncode
        if returncode == 0:
            return True
        else:
            print(f"Command '{command}' returned a non-zero exit code.")
            return False
    except FileNotFoundError:
        print(f"Command '{command}' not found.")
        return False
    except Exception as e:
        print(f"An error occurred while trying to run '{command}': {e}")
        return False


def initialise_recur(iqtree_version: Optional[str] = None):
    my_env, local_iqtree2_path, bin_dir, conda_prefix = setup_environment()
    system = platform.system()

    try:
        if system == "Windows":
            if mp.get_start_method(allow_none=True) != 'spawn':
                mp.set_start_method('spawn', force=True)
        elif system in ["Linux", "Darwin"]:
            if mp.get_start_method(allow_none=True) != 'fork':
                mp.set_start_method('fork')
    except RuntimeError as e:
        print(f"Multiprocessing context setting error on {system}: {e}")
        pass

    if system == "Linux" and iqtree_version == 'local':
        my_env['PATH'] = bin_dir + os.pathsep + my_env['PATH']
        if CanRunCommand(f"{local_iqtree2_path} --version"):
            return

    if iqtree_version == "conda" and conda_prefix:
        iqtree_path = shutil.which("iqtree2")
        if iqtree_path:
            print("\nConda version of IQ-TREE2 found.")
            print(f"IQ-TREE2 path: {iqtree_path}")
            return

    if iqtree_version == "system" and not conda_prefix:
        iqtree_path = shutil.which("iqtree2")
        if iqtree_path:
            print("\nSystem-wide version of IQ-TREE2 found.")
            print(f"IQ-TREE2 path: {iqtree_path}")
            return

    if conda_prefix and CanRunCommand("iqtree2 --version"):
        print("Local IQ-TREE2 binary failed to run, falling back to the conda version.")
        return

    if CanRunCommand("iqtree2 --version"):
        print("Local IQ-TREE2 binary failed to run, falling back to system-wide binary.")
        return

    print("Cannot proceed. IQ-TREE2 does not exist in either local bin or system-wide PATH.")
    print("Please ensure IQ-TREE2 is properly installed before running RECUR!\n")
    sys.exit(1)


def initialise_terminate_flag() -> threading.Event:
    manager = mp.Manager()
    return manager.Event()

def reset_terminate_flag(event: threading.Event) -> None:
    event.clear()

def initialise_worker(event: threading.Event) -> None:
    setup_signal_handler(event)

def signal_handler(sig: int, frame, event: threading.Event) -> None:
    print(f"Received signal {sig}")
    event.set()
    gc.collect()
    sys.exit(0)

def setup_signal_handler(event: Optional[threading.Event] = None) -> threading.Event:
    if os.name != 'nt':
        signals = [signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGUSR1, signal.SIGUSR2]
    else:
        signals = [signal.SIGTERM, signal.SIGINT]

    if event is None:
        event = initialise_terminate_flag()

    for sig in signals:
        signal.signal(sig, lambda s, f: signal_handler(s, f, event))

    return event

def cleanup() -> None:
    try:
        gc.collect()
    except Exception as e:
        print(f"ERROR during cleanup: {e}")
        print(traceback.format_exc())
    sys.exit()

def kill_child_processes(parent_pid: int, sig: signal.Signals = signal.SIGTERM, retries: int = 3) -> None:
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        os._exit(1)  # Immediate exit without cleanup

    children = parent.children(recursive=True)
    for process in children:
        for attempt in range(retries):
            try:
                process.send_signal(sig)
                process.wait(timeout=3)
                if not process.is_running():
                    break
            except psutil.NoSuchProcess:
                break
            except psutil.TimeoutExpired:
                continue
            except Exception as e:
                print(f"ERROR killing process {process.pid}: {e}")
                break
        else:
            print(f"Failed to kill process {process.pid} after {retries} attempts")
            print(f"Force killing process {process.pid} after {retries} attempts")
            try:
                process.kill()
                process.wait(timeout=3)
                if not process.is_running():
                    print(f"Successfully force killed process {process.pid}")
                    break
            except psutil.NoSuchProcess:
                print(f"Process {process.pid} does not exist")
                break
            except Exception as e:
                print(f"ERROR force killing process {process.pid}: {e}")
                break

    os._exit(0) 

def combine_mut_matrix_and_recurrence_list(mut_matrices: Dict[int, np.typing.NDArray[np.int64]],
                                           recurrence_list: List[List[Union[str, int, float]]],
                                           combined_prot_seqs_dict: Dict[str, str],
                                           species_of_interest: List[str],
                                           residue_dict_flip: Dict[int, str]) -> List[List[Union[str, int, float]]]:
    
    extant_seq = {species: seq for species, seq in combined_prot_seqs_dict.items() if species in species_of_interest}
    ident_dict = {}

    for rec_loc, res in enumerate(zip(*extant_seq.values())):
        ident_dict[rec_loc] = res

    for rec_list in recurrence_list:
        res_loc = int(rec_list[0])
        mut = mut_matrices[res_loc]
        coo = spy.sparse.coo_matrix(mut)
        data, row, col = coo.data, coo.row, coo.col
        sorted_indices = np.argsort(data)[::-1]
        data = data[sorted_indices]
        row = row[sorted_indices]
        col = col[sorted_indices]

        data_str_list = [*(map(str, data))]
        parent_child = []
        for row, col in zip(row, col):
            parent = residue_dict_flip[row]
            child = residue_dict_flip[col]
            parent_child.append(">".join((parent, child)))
        parent_child_data = [*zip(parent_child, data_str_list)]
        parent_child_data_str = ",".join([":".join(pcd) for pcd in parent_child_data])
        rec_list.append(parent_child_data_str)

        res_freq = [*Counter(ident_dict[res_loc]).items()]
        res_freq = sorted(res_freq, reverse=True, key=lambda x: x[1])
        res_freq_str = ",".join([":".join((res, str(freq))) for res, freq in res_freq])
        rec_list.append(res_freq_str)

    return recurrence_list

def get_subtree_species(node):

    species = []
    for leaf in node.leaf_iter():
        species.append(leaf.taxon.label)
    return species


def ParentChildRelation(t, 
                        outgroup_species, 
                        sequence_dict: Dict[str, str], 
                        residue_dict: Dict[str, int], 
                        res_loc: int, 
                        terminate_flag: threading.Event) -> Tuple[int, Optional[str], List[str], np.typing.NDArray[np.int64], str]:
    try:
        t.is_rooted = True
        if len(outgroup_species) == 1:
            outgroup_mrca = t.find_node_with_taxon_label(outgroup_species[0])
        else:
            outgroup_mrca = t.mrca(taxon_labels=outgroup_species)

        t.reroot_at_edge(outgroup_mrca.edge, update_bipartitions=False)
        rt = t.seed_node
        rt_children = rt.child_nodes()
        
        root_of_interest = next((child for child in rt_children if child != outgroup_mrca), None)

        if root_of_interest is None:
 
            return res_loc, None, outgroup_species, np.zeros((20, 20), dtype=np.int64), "No root of interest found; tree structure might be incorrect."

        root_node = root_of_interest.label
        root_node = root_node if "/" not in root_node else root_node.split("/")[0]

        outgroup_subtree_species = get_subtree_species(outgroup_mrca)

        residue_mut = np.zeros((20, 20), dtype=np.int64)
        taxon_count = 0
        for nd in root_of_interest.postorder_iter():
            if terminate_flag.is_set():
                return res_loc, root_node, outgroup_species, np.zeros((20, 20), dtype=np.int64), "Terminated"

            if nd.parent_node is None:
                continue

            if nd.is_internal():
                child = nd.label
            else:
                child = nd.taxon.label

            if child is not None and nd.parent_node.label is not None:
                parent = nd.parent_node.label
                child = child if "/" not in child else child.split("/")[0]
                parent = parent if "/" not in parent else parent.split("/")[0]

                child_ident = sequence_dict[child][res_loc]
                parent_ident = sequence_dict[parent][res_loc]

                child_pos = residue_dict.get(child_ident)  # implicitly avoiding -, *, ., etc.
                parent_pos = residue_dict.get(parent_ident)

                if (child_pos is not None and parent_pos is not None) and (child_pos != parent_pos):
                    residue_mut[parent_pos, child_pos] += 1

                taxon_count += 1

        return taxon_count, root_node, outgroup_subtree_species, residue_mut, ""
    except BrokenPipeError:
        print("Broken pipe error while checking termination flag.")
        return res_loc, root_node, outgroup_species, np.zeros((20, 20), dtype=np.int64), "Broken pipe error"
    except Exception as e:
        error_msg = f"ERROR in ParentChildRelation: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return res_loc, root_node, outgroup_species, np.zeros((20, 20), dtype=np.int64), error_msg

def CountMutations(treefile: str,
                   outgroup_species: List[str],
                   sequence_dict: Dict[str, str], 
                   residue_dict: Dict[str, int], 
                   n_species: int,
                   res_loc: int,
                   terminate_flag: threading.Event) -> Tuple[int, Optional[str], List[str], np.typing.NDArray[np.int64], str]:
    try:
        with open(treefile, 'r') as f:
            t = dendropy.Tree.get(file=f, schema="newick")

        taxon_count, root_node, outgroup_subtree_species, residue_mut, error_msg = ParentChildRelation(t, outgroup_species, sequence_dict, residue_dict, res_loc, terminate_flag)
        
        if error_msg:
            return res_loc, root_node, outgroup_species, np.zeros((20, 20), dtype=np.int64), error_msg
        
        if len(outgroup_subtree_species) != len(outgroup_species):
            outgroup_species = outgroup_subtree_species

        expected_relationships = 2 * (n_species - len(outgroup_species)) - 2
        if taxon_count != expected_relationships:
            return res_loc, root_node, outgroup_species, residue_mut, f"Taxon count {taxon_count} does not match expected relationships {expected_relationships}"
        
        return res_loc, root_node, outgroup_species, residue_mut, ""
    except Exception as e:
        error_msg = f"ERROR in CountMutations for treefile {treefile}, res_loc {res_loc}: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return res_loc, None, outgroup_species, np.zeros((20, 20), dtype=np.int64), error_msg


def count_mutations(treefile: str, 
                    sequence_dict: Dict[str, str], 
                    residue_dict: Dict[str, int], 
                    protein_len: int,
                    n_species: int, 
                    outgroups: List[str], 
                    nthreads: int, 
                    production_logger: logging.Logger,
                    terminate_flag: threading.Event,
                    batch_size: Optional[int] = None) -> Tuple[Optional[str], List[str], Dict[int, np.typing.NDArray[np.int64]]]:
    
    all_res_locs = range(protein_len)
    total_memory = psutil.virtual_memory().total

    if batch_size is None:
        sample_data = list(all_res_locs)[:10]  # Use a sample of locations to estimate
        batch_size = util.determine_optimal_batch_size(sample_data, 'mutation_count', total_memory)
    batches = [all_res_locs[i:i + batch_size] for i in range(0, len(all_res_locs), batch_size)]

    mut_results_dict = {}
    root: Optional[str] = None
    outgroup_species: Optional[List[str]] = None

    try:
        with ProcessPoolExecutor(max_workers=nthreads, initializer=initialise_worker, initargs=(terminate_flag,)) as executor:
            futures = [executor.submit(CountMutations, treefile, outgroups, sequence_dict, residue_dict, n_species, res_loc, terminate_flag) for batch in batches for res_loc in batch]

            for future in as_completed(futures):
                if terminate_flag.is_set():
                    break
                try:
                    result = future.result()
                    res_loc, batch_root, batch_outgroup_species, residue_mut, error_msg = result
                    if error_msg:
                        print(error_msg)
                        production_logger.error(error_msg)
                        terminate_flag.set()
                        break

                    if root is None:
                        root = batch_root
                    mut_results_dict[res_loc] = residue_mut
                    outgroup_species = batch_outgroup_species

                except BrokenPipeError:
                    error_msg = "Broken pipe error during parallel processing."
                    print(error_msg)
                    production_logger.error(error_msg)
                    terminate_flag.set()
                    break
                except Exception as e:
                    error_msg = f"ERROR during parallel processing: {e}"
                    print(error_msg)
                    print(traceback.format_exc())
                    production_logger.error(error_msg)
                    terminate_flag.set()
                    break

    except Exception as e:
        error_msg = f"ERROR during parallel processing: {e}"
        print(error_msg)
        print(traceback.format_exc())
        production_logger.error(error_msg)
    finally:
        if outgroup_species and len(outgroups) != len(outgroup_species):
            warnings.warn(f"Outgroup species provided not monophyletic. Outgroups will be updated. Please find the updated outgroups in the log file.")
            production_logger.info(f"Updated outgroups: {outgroup_species}", extra={'to_file': True, 'to_console': False})
            outgroups = outgroup_species
        return root, outgroups, mut_results_dict

def GetRecurrenceList(rec_loc: int,
                      mut_matrix: np.typing.NDArray[np.int64], 
                      terminate_flag) -> Union[List[List[Union[int, str, float]]], str]:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            recurrence_list: List[List[Union[int, str, float]]] = []
            rows, cols = np.where(mut_matrix > 1)  # Changed condition from > 0 to > 1 as per original code
            for i in range(len(rows)):
                if terminate_flag.is_set():
                    print(f"Terminating GetRecurrenceList for rec_loc {rec_loc} due to termination signal")
                    return "Terminated"
                r = rows[i]
                c = cols[i]
                recurrence = mut_matrix[r, c]
                flipflop = mut_matrix[c, r]
                parent = util.residues[r]
                child = util.residues[c]
                recurrence_list.append([rec_loc, parent, child, recurrence, flipflop])
        
        return recurrence_list
    except Exception as e:
        error_msg = f"ERROR in GetRecurrenceList for rec_loc {rec_loc}: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return error_msg

def get_recurrence_list(mut_matrices: Dict[int, np.typing.NDArray[np.int64]], 
                        extant_seq: Dict[str, str], 
                        outgroups: List[str],
                        nthreads: int,
                        terminate_flag: threading.Event) -> List[List[Union[str, int, float]]]:

    if len(extant_seq) == 0 or len(outgroups) == 0:
        error_msg = "ERROR: extant_seq or outgroups is None. Cannot proceed with get_recurrence_list."
        print(error_msg)
        return []
    
    tasks = [(rec_loc, mut_matrix, terminate_flag) for rec_loc, mut_matrix in mut_matrices.items()]

    recurrence_list = []
    try:
        with ProcessPoolExecutor(max_workers=nthreads, initializer=initialise_worker, initargs=(terminate_flag,)) as executor:
            futures = {executor.submit(GetRecurrenceList, *task): task for task in tasks}
            for future in as_completed(futures):
                if terminate_flag.is_set():
                    print("Terminating get_recurrence_list due to termination signal")
                    break
                try:
                    result = future.result()
                    if isinstance(result, str) and result.startswith("ERROR"):
                        print(result)
                    elif isinstance(result, list) and len(result) > 0:
                        recurrence_list.extend(result)
                except BrokenPipeError:
                    error_msg = "Broken pipe error during recurrence list processing."
                    print(error_msg)
                    terminate_flag.set()
                    break
                except Exception as e:
                    error_msg = f"ERROR during recurrence list processing: {e}"
                    print(error_msg)
                    print(traceback.format_exc())
                    terminate_flag.set()
                    break  # Stop processing further futures
                
            # wait(futures)
    except Exception as e:
        error_msg = f"ERROR during get_recurrence_list: {e}"
        print(error_msg)
        print(traceback.format_exc())

    return recurrence_list


def WorkerProcessAndCount(file: str, 
                          mcs_alnDir: str,
                          mcs_treefile: str, 
                          isnuc_fasta: bool, 
                          residue_dict: Dict[str, int], 
                          n_species: int, 
                          outgroup_mrca: List[str], 
                          res_loc_set: Set[int],
                          production_logger: logging.Logger,
                          terminate_flag: threading.Event) -> Tuple[Optional[Dict[int, np.typing.NDArray[np.int64]]], Optional[str]]:
    try:
        file_path = os.path.join(mcs_alnDir, file)
        mcs_combined_prot_seqs_dict, _ = files.FileReader.ReadAlignment(file_path)

        if isnuc_fasta:
            mcs_combined_prot_seqs_dict, _ = util.GetSeqsDict(mcs_combined_prot_seqs_dict, "CODON1")

        if terminate_flag.is_set():
            msg = f"Terminating processing for file: {file} due to termination signal"
            return None, msg

        mut_matrices = {}
        root_node = None
        for res_loc in res_loc_set:
            if terminate_flag.is_set():
                msg = f"Terminating processing for file: {file} at res_loc {res_loc} due to termination signal"
                return None, msg

            res_loc, current_root_node, _, mut_matrix, error_msg = CountMutations(mcs_treefile,
                                                                    outgroup_mrca,
                                                                    mcs_combined_prot_seqs_dict,
                                                                    residue_dict,
                                                                    n_species,
                                                                    res_loc,
                                                                    terminate_flag)

            if error_msg:
                print(error_msg)
                production_logger.error(error_msg)
                continue

            if current_root_node is None:
                msg = f"CountMutations returned None for root_node at res_loc {res_loc}"
                print(msg)
                continue

            mut_matrices[res_loc] = mut_matrix
            root_node = current_root_node

        if root_node:
            return mut_matrices, None
        else:
            error_msg = f"No valid root_node found for file: {file}, returning None."
            print(error_msg)
            return mut_matrices, error_msg

    except BrokenPipeError:
        error_msg = "Broken pipe error while processing file."
        print(error_msg)
        return None, error_msg
    except Exception as e:
        error_msg = f"ERROR in WorkerProcessAndCount for mcs task {file}: {e}"
        print(error_msg)
        print(traceback.format_exc())
        production_logger.error(error_msg)
        return None, error_msg
    

def process_mcs_files_in_chunks(mcs_alnDir: str, 
                                mcs_treefile: str, 
                                residue_dict: Dict[str, int], 
                                n_species: int, 
                                outgroup_mrca: List[str], 
                                nthreads: int, 
                                isnuc_fasta: bool, 
                                res_loc_set: Set[int],
                                production_logger: logging.Logger,
                                terminate_flag: threading.Event,
                                batch_size: Optional[int] = None) -> List[Dict[int, np.typing.NDArray[np.int64]]]:

    mcs_files = os.listdir(mcs_alnDir)
    total_files_count = len(mcs_files)
    total_memory = psutil.virtual_memory().total

    if batch_size is None:
        sample_file = os.path.join(mcs_alnDir, mcs_files[0])
        batch_size = util.determine_optimal_batch_size(sample_file, 'file_processing', total_memory)
    results = []
    worker = partial(WorkerProcessAndCount, 
                     mcs_alnDir=mcs_alnDir,
                     mcs_treefile=mcs_treefile,
                     isnuc_fasta=isnuc_fasta,
                     residue_dict=residue_dict, 
                     n_species=n_species, 
                     outgroup_mrca=outgroup_mrca, 
                     res_loc_set=res_loc_set,
                     production_logger=production_logger,
                     terminate_flag=terminate_flag)
    
    batches = [mcs_files[i:i + batch_size] for i in range(0, total_files_count, batch_size)]

    try:
        with ProcessPoolExecutor(max_workers=nthreads) as executor:
            futures = [executor.submit(worker, file_data) for batch in batches for file_data in batch]
            
            for future in as_completed(futures):
                if terminate_flag.is_set():
                    print("Terminating processing due to termination signal")
                    break
                try:
                    result, error_msg = future.result()
                    if error_msg:
                        print(error_msg)
                        production_logger.error(error_msg)
                        terminate_flag.set()
                        break

                    if result:
                        results.append(result)

                except Exception as e:
                    error_msg = f"ERROR during processing: {e}"
                    print(error_msg)
                    print(traceback.format_exc())
                    production_logger.error(error_msg)
                    terminate_flag.set()
                    break

    except Exception as e:
        error_msg = f"ERROR during processing files in chunks: {e}"
        print(error_msg)
        print(traceback.format_exc())
        production_logger.error(error_msg)

    return results


def compute_p_value_for_rec_list(rec_list: List[Union[str, int, float]], 
                                 mcs_results: Dict[int, List[np.typing.NDArray[np.int64]]], 
                                 residue_dict: Dict[str, int], 
                                 nalign: int,
                                 terminate_flag: threading.Event) -> Union[List[Union[str, int, float]], None]:

    try:
        if terminate_flag.is_set():
            print(f"Terminating p-value computation for rec_list: {rec_list} due to termination signal")
            return None

        rec_loc = int(rec_list[0])
        parent = str(rec_list[1])
        child = str(rec_list[2])
        recurrence = int(rec_list[3])

        r = residue_dict[parent]
        c = residue_dict[child]
        count_greater = 0
        mcs_result = mcs_results.get(rec_loc)

        if mcs_result is not None:
            for m in mcs_result:
                if np.all(m == 0):
                    continue
                mcs_rec = m[r, c]
                if mcs_rec >= recurrence:
                    count_greater += 1

        p_value = (count_greater + 1) / (nalign + 1)
        rec_list.append(np.round(p_value, len(str(nalign))))
        # print(f"{rec_list} - {count_greater}")
        return rec_list

    except Exception as e:
        error_msg = f"ERROR computing p-value for rec_list {rec_list}: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return None


def compute_p_values(mcs_results: List[Dict[int, np.typing.NDArray[np.int64]]],
                     recurrence_list: List[List[Union[str, int, float]]],
                     residue_dict: Dict[str, int],
                     nalign: int,
                     terminate_flag: threading.Event) -> List[List[Union[str, int, float]]]:

    recurrence_list_pvalue = []

    mcs_results_dict = defaultdict(list)
    
    for mcs_result in mcs_results:
        for res_loc, val in mcs_result.items():
            mcs_results_dict[res_loc].append(val)

    partial_compute = partial(compute_p_value_for_rec_list, 
                            mcs_results=mcs_results_dict, 
                            residue_dict=residue_dict, 
                            nalign=nalign,
                            terminate_flag=terminate_flag)

    try:
        for rec_list in recurrence_list:
            if terminate_flag.is_set():
                print(f"Terminating p-value computation due to termination signal")
                break
            try:
                result = partial_compute(rec_list)
                if result is not None:
                    recurrence_list_pvalue.append(result)
            except Exception as e:
                error_msg = f"ERROR during p-value computation for rec_list {rec_list}: {e}"
                print(error_msg)
                print(traceback.format_exc())
                terminate_flag.set()
                break 

        recurrence_list_pvalue.sort(key=lambda x: (float(x[-1]), -float(x[-3])))
    except Exception as e:
        error_msg = f"ERROR during compute_p_values: {e}"
        print(error_msg)
        print(traceback.format_exc())

    return recurrence_list_pvalue

def main(args: Optional[List[str]] = None):

    d_results = None 
    production_logger = None
    terminate_flag = initialise_terminate_flag()
    setup_signal_handler(terminate_flag)

    if not args:
        args = sys.argv[1:]

    if not args or len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.PrintHelp()
        sys.exit()

    elif args[0] == "-v" or args[0] == "--version":
        print(f"RECUR:v{__version__}")
        sys.exit()
 
    # util.get_system_info()
    start_main = time.perf_counter()

    try:
        input_command = "RECUR command: " + " ".join(sys.argv)  + "\n"
        
        options, alnDir, alnPath, resultsDir_nonDefault = process_args.ProcessArgs(args)
        iqtree_version = options.iqtree_version if options.iqtree_version else "system"
        if not os.getenv('ENV_SETUP_DONE'):
            initialise_recur(iqtree_version)
            os.environ['ENV_SETUP_DONE'] = '1'


        aln_path_dict = files.FileHandler.ProcessesNewAln(alnDir, alnPath)

        aln_len = len(aln_path_dict)
        if isinstance(options.gene_tree, dict):
            if len(options.gene_tree) <= aln_len and len(options.gene_tree) == 1:
                options.gene_tree = {gene: [*options.gene_tree.values()][0] for gene in aln_path_dict}
        
        if isinstance(options.outgroups, dict):
            if len(options.outgroups) <= aln_len and len(options.outgroups) == 1:
                options.outgroups = {gene: [*options.outgroups.values()][0] for gene in aln_path_dict}

        count = 0
        asr = True

        for gene, aln_path in aln_path_dict.items():

            if terminate_flag.is_set():
                break
            
            try:    
                start = time.perf_counter()

                if count > 0:
                    try:
                        width = os.get_terminal_size().columns
                    except OSError as e:
                        width = 80
                    util.print_centered_text(width, f"Processing gene {gene}")
                count += 1

                if isinstance(options.outgroups, dict):
                    outgroup_mrca = options.outgroups[gene]
                else:
                    outgroup_mrca = options.outgroups
                
                filehandler = files.FileHandler()
                filereader = files.FileReader()
                filewriter = files.FileWriter()

                alnFN = os.path.basename(aln_path)
                filehandler.gene_of_interest = alnFN
                base_dir = os.path.join(resultsDir_nonDefault, f"{alnFN}.recur") if resultsDir_nonDefault else os.path.join(alnDir, f"{alnFN}.recur")
                if not os.path.exists(base_dir): 
                    os.mkdir(base_dir)
            
                filehandler.CreateOutputDirectories(options, base_dir)

                results_dir = filehandler.GetResultsDirectory()
                production_logger = util.setup_logging(results_dir, "w", "brief")
                
                production_logger.info(f"{input_command}", extra={'to_file': True, 'to_console': False})
                print()
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Starting RECUR v%s" % __version__, extra={'to_file': True, 'to_console': True})

                alignment_dict, alignment_len = filereader.ReadAlignment(aln_path)
                n_species = len(alignment_dict)

                production_logger.info(f"Analysing: {gene}", extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Number of extant species found: {n_species}", extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Length of the {gene} alignment: {alignment_len}", extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Results Directory: {results_dir}\n", extra={'to_file': True, 'to_console': True})

                species_of_interest = [*alignment_dict.keys()]
                isnuc = util.CheckSequenceType([*alignment_dict.values()])

                if isnuc and options.sequence_type == "AA":
                    print("\nWARNING: The input sequence type does not match the existing alignments. RECUR will convert the codon alignments into protein alignments with the default CODON model\n")
                    options.sequence_type = "CODON1"

                check_exist = 0
                real_phyDir = filehandler.GetRealPhylogenyDir()
                for file in os.listdir(real_phyDir):
                    if file.rsplit(".", 1)[1] in ["state", "treefile", "iqtree"]:
                        check_exist += 1
                
                restart_step1 = False
                gene_tree = options.gene_tree[gene] if isinstance(options.gene_tree, dict) else None
                if not gene_tree:
                    step1_info = "Step1: Inferring phylogenetic tree and model of evolution"
                else:
                    step1_info = "Step1: Inferring ancestral sequences, phylogenetic tree and model of evolution"
                production_logger.info(step1_info, extra={'to_file': True, 'to_console': True})
                production_logger.info("="*len(step1_info), extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Results Directory: {real_phyDir}\n", extra={'to_file': True, 'to_console': True})


                if (options.usr_state and options.usr_tree and options.usr_iqtree):
                    statefile = options.usr_state
                    treefile = options.usr_tree
                    iqtreefile = options.usr_iqtree
                    production_logger.info("NOTE: with the provided statefile and treefile, RECUR will skip step1.\n", extra={'to_file': True, 'to_console': False})
                elif check_exist == 3 and options.restart_from != 1:
                    production_logger.info("NOTE: with the existing statefile and treefile, RECUR will skip step1.\n", extra={'to_file': True, 'to_console': False})
                    statefile = filehandler.GetStateFileFN() 
                    treefile = filehandler.GetTreeFileFN() 
                    iqtreefile = filehandler.GetIQTreeFileFN()
                else:
                    if options.restart_from == 1:
                        production_logger.info("Restart RECUR from step1\n", extra={'to_file': True, 'to_console': False})
                    
                    if check_exist > 0 and check_exist < 3:
                        production_logger.info("NOTE: RECUR is forced to restart from step1 due some missing files\n", extra={'to_file': True, 'to_console': False})
                    restart_step1 = True
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Ran IQ-TREE to build the gene trees for the real phylogeny", extra={'to_file': True, 'to_console': True})
                    production_logger.info("Using %d RECUR thread(s), %d IQ-TREE2 thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads), extra={'to_file': True, 'to_console': True})

                    if not gene_tree:
                        asr = False

                    commands = run_commands.GetGeneTreeBuildCommands([aln_path], 
                                                        real_phyDir, 
                                                        options.evolution_model,
                                                        options.iqtree_nthreads,
                                                        phy_seed=options.seed,
                                                        asr=asr,
                                                        sequence_type=options.sequence_type,
                                                        gene_tree=gene_tree,
                                                        bootstrap=options.bootstrap,
                                                        sh_alrt=options.bootstrap)
                    
                    if not gene_tree:
                        production_logger.info(f"step1 iqtree2 gene tree building command: ", extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    else:
                        production_logger.info(f"step1 iqtree2 ancestral state reconstruction command: ", extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    run_commands.RunCommand(commands, 
                                            real_phyDir, 
                                            nthreads=options.recur_nthreads,
                                            delete_files=True,
                                            files_to_keep=["state", "treefile", "iqtree"],
                                            fd_limit=options.fd_limit)
                    
                    treefile = filehandler.GetTreeFileFN() 
                    iqtreefile = filehandler.GetIQTreeFileFN()
                    

                best_evolution_model = filereader.ReadIQTreeFile(iqtreefile)

                if not gene_tree:
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + f"Tree inference complete, best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file':True, 'to_console': True})
            
                else:
                    statefile = filehandler.GetStateFileFN()
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + f"Ancestral state reconstruction complete, best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file': True, 'to_console': True})

                if not gene_tree and restart_step1:

                    step2_info = f"Step2: Inferring ancestral sequences"
                    production_logger.info(step2_info, extra={'to_file': True, 'to_console': True})
                    production_logger.info("="*len(step2_info), extra={'to_file': True, 'to_console': True})
                    step2_results_info = f"Results Directory: {filehandler.GetRealPhylogenyDir()}\n"
                    production_logger.info(step2_results_info, extra={'to_file': True, 'to_console': True})
                
                    filehandler.UpdateTreeFile(treefile)

                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + f"Starting ancestral state reconstruction.", extra={'to_file': True, 'to_console': True})

                    commands = run_commands.GetGeneTreeBuildCommands([aln_path], 
                                                        real_phyDir, 
                                                        best_evolution_model,
                                                        options.iqtree_nthreads,
                                                        phy_seed=options.seed,
                                                        sequence_type=options.sequence_type,
                                                        gene_tree=treefile,
                                                        asr=True)
                    
                    production_logger.info(f"step2 iqtree2 ancestral state reconstruction command: ",  extra={'to_file': True, 'to_console': False})
                    production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})
    
                    run_commands.RunCommand(commands, 
                                            real_phyDir,
                                            nthreads=options.recur_nthreads, 
                                            delete_files=True,
                                            files_to_keep=["state", "treefile", "iqtree"],
                                            fd_limit=options.fd_limit)

                    statefile = filehandler.GetStateFileFN()
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + f"Ancestral state reconstruction complete\n", extra={'to_file': True, 'to_console': True})

                filehandler.CheckFileCorrespondance(gene, statefile, treefile)
                residue_dict, residue_dict_flip = util.residue_table()
                
                node_seq_dict = filereader.ReadStateFile(statefile)

                combined_seq_dict = {k: v for d in (node_seq_dict, alignment_dict) for k, v in d.items()}
                alignment_dict.clear()

                if options.sequence_type == "AA":
                    node_prot_seqs_fn = filehandler.GetNodeProtSeqsFN()
                    filewriter.WriteSeqsToAln(node_seq_dict, node_prot_seqs_fn)
                    combined_prot_seqs_dict = combined_seq_dict.copy()
                    protein_len = len([*node_seq_dict.values()][0])
                else:
                    node_dna_seqs_fn = filehandler.GetNodeDNASeqsFN()
                    filewriter.WriteSeqsToAln(node_seq_dict, node_dna_seqs_fn)
                    combined_prot_seqs_dict, protein_len = util.GetSeqsDict(combined_seq_dict, options.sequence_type)
                    node_prot_seqs_dict, _ = util.GetSeqsDict(node_seq_dict, options.sequence_type)
                    node_prot_seqs_fn = filehandler.GetNodeProtSeqsFN()
                    filewriter.WriteSeqsToAln(node_prot_seqs_dict, node_prot_seqs_fn)
                
                node_seq_dict.clear()
                combined_seq_dict.clear()

                try:
                    root_node, outgroup_mrca, mut_matrices = count_mutations(treefile, 
                                                                            combined_prot_seqs_dict,
                                                                            residue_dict, 
                                                                            protein_len, 
                                                                            n_species,
                                                                            outgroup_mrca, 
                                                                            options.nthreads,
                                                                            production_logger,
                                                                            terminate_flag,
                                                                            batch_size=options.batch_size)
                finally:
                    reset_terminate_flag(terminate_flag)

                production_logger.info(f"Root of species of interest: {root_node}", extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Substitution matrix output: {filehandler.GetMutMatrixDir()}\n", extra={'to_file': True, 'to_console': True})

                filewriter.WriteMutMatrix(mut_matrices, 
                                          residue_dict_flip, 
                                          filehandler.GetMutCountMatricesFN(),
                                          filehandler.GetAccumMutCountMatricesFN())
                
                try:
                    recurrence_list = get_recurrence_list(mut_matrices, 
                                                        combined_prot_seqs_dict,
                                                        outgroup_mrca, 
                                                        options.nthreads,
                                                        terminate_flag)
                finally:
                    reset_terminate_flag(terminate_flag)

                if len(recurrence_list) == 0:
                    production_logger.info(f"ATTENTION: No recurrence has identified for gene {gene}! Monte-Carlo Simiatlion will be SKIPPED!\n", 
                                           extra={'to_file': True, 'to_console': True})
                    continue 


                res_loc_set = set(int(res_list[0]) for res_list in recurrence_list)

                if not gene_tree:
                    step2_info = f"Step3: Simulating Sequence Evolution with {options.nalign} replicates"
                else: 
                    step2_info = f"Step2: Simulating Sequence Evolution with {options.nalign} replicates"
                    
                production_logger.info(step2_info, extra={'to_file': True, 'to_console': True})
                production_logger.info("="*len(step2_info), extra={'to_file': True, 'to_console': True})
                step2_results_info = f"Results Directory: {filehandler.GetMCSimulationDir()}\n"
                production_logger.info(step2_results_info, extra={'to_file': True, 'to_console': True})

                mcs_faDir = filehandler.GetMCSimulationDir()
                if len(os.listdir(mcs_faDir)) != options.nalign or \
                    (options.restart_from == 2 or options.restart_from == 1):

                    if len(os.listdir(mcs_faDir)) != options.nalign and len(os.listdir(mcs_faDir)) > 0:
                        util.delete_files_in_directory(mcs_faDir)

                    identifier = "rooted_" + gene + "_alisim"
                    output_prefix = os.path.join(mcs_faDir, identifier)
                    if root_node is not None:
                        if options.sequence_type == "AA":
                            fn_root_node = ",".join((node_prot_seqs_fn, root_node))
                        else:
                            fn_root_node = ",".join((node_dna_seqs_fn, root_node))
                    else:
                        raise ValueError("Root node of interest is None.")
                    
                    mcs_commands = run_commands.GetMCsimulationCommand(output_prefix,
                                                                    options.iqtree_nthreads,
                                                                    options.mcs_seed, 
                                                                    best_evolution_model,
                                                                    treefile, 
                                                                    fn_root_node,
                                                                    options.nalign)
                    
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Starting Monte-Carlo Simulation.", extra={'to_file': True, 'to_console': True})
                    production_logger.info("Using %d RECUR thread(s), %d IQ-TREE2 thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads), 
                                        extra={'to_file': True, 'to_console': False})
                    production_logger.info("step2 iqtrees command: ", extra={'to_file': True, 'to_console': False})
                    production_logger.info(f"{mcs_commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    run_commands.RunCommand(mcs_commands, 
                                            mcs_faDir, 
                                            nthreads=options.recur_nthreads,
                                            delete_files=True,
                                            files_to_keep=["fasta", "fa"],
                                            fd_limit=options.fd_limit)
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Monte-Carlo Simulation complete.\n", extra={'to_file': True, 'to_console': True})
                else:
                    production_logger.info("NOTE: With the existing Monte-Carlo simulated *.fa files, RECUR will skip step2.\n", 
                                        extra={'to_file': True, 'to_console': True})
                    
                for file in os.listdir(real_phyDir):
                    file_path = os.path.join(real_phyDir, file)
                    if os.path.exists(file_path):
                        if file.endswith(".treefile.txt") or file.endswith(".treefile.log"):
                            os.remove(file_path)

                afasta = random.choice(os.listdir(mcs_faDir))
                fasta_dict, fasta_len = filereader.ReadAlignment(os.path.join(mcs_faDir, afasta))
                isnuc_fasta = util.CheckSequenceType([*fasta_dict.values()])

                if options.usr_mcs_alnDir:
                    mcs_alnDir = options.usr_mcs_alnDir
                else:
                    mcs_alnDir = mcs_faDir

                if not options.recDir:
                    recurrenceDir = alnDir
                else:
                    recurrenceDir = options.recDir
                if not gene_tree:
                    step3_info = f"Step4: Analysing recurrent substitutions"
                else:
                    step3_info = f"Step3: Analysing recurrent substitutions"
                production_logger.info(step3_info, extra={'to_file': True, 'to_console': True})
                production_logger.info("="*len(step3_info), extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Results Directory: {recurrenceDir}\n", extra={'to_file': True, 'to_console': True})

                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Starting create substitution matrices for simulated phylogeny.", extra={'to_file': True, 'to_console': True})
                production_logger.info("Using %d thread(s) for RECUR analysis" % options.nthreads, extra={'to_file': True, 'to_console': True})
                
                try:
                    mcs_results = process_mcs_files_in_chunks(mcs_alnDir, 
                                                        treefile, 
                                                        residue_dict, 
                                                        n_species, 
                                                        outgroup_mrca, 
                                                        options.nthreads,
                                                        isnuc_fasta,
                                                        res_loc_set, 
                                                        production_logger,
                                                        terminate_flag,
                                                        batch_size=options.batch_size)
                finally:
                    reset_terminate_flag(terminate_flag)

                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Substitution matrices creation complete.\n")
                
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Starting compute p values.")

                try:
                    recurrence_list_pvalue = compute_p_values(mcs_results, 
                                                            recurrence_list,
                                                            residue_dict, 
                                                            options.nalign,
                                                            terminate_flag)
                
                finally:
                    reset_terminate_flag(terminate_flag)

                recurrence_list_updated = combine_mut_matrix_and_recurrence_list(mut_matrices,
                                                                                 recurrence_list_pvalue,
                                                                                 combined_prot_seqs_dict,
                                                                                 species_of_interest,
                                                                                 residue_dict_flip)
                
                filewriter.WriteRecurrenceList(recurrence_list_updated, filehandler.GetRecurrenceListFN(recurrenceDir))
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "p values computing complete.", extra={'to_file': True, 'to_console': True})

                d_results = os.path.normpath(filehandler.GetResultsDirectory()) + os.path.sep
                rec_results = os.path.normpath(filehandler.GetRecurrenceListFN(recurrenceDir))
                production_logger.info("\nResults:\n    %s\n" % rec_results, extra={'to_file': True, 'to_console': True})
                
                del combined_prot_seqs_dict, alignment_dict, mut_matrices, recurrence_list
                del recurrence_list_pvalue, recurrence_list_updated
                gc.collect()

                util.log_memory_usage(f"after processing gene {gene}", production_logger)

                end = time.perf_counter()
                duration = end - start
                if production_logger:
                    production_logger.info(f"Finished analysis of {gene} in {duration:.2f} seconds", extra={'to_file': True, 'to_console': True})
                else:
                    print(f"Finished analysis of {gene} in {duration:.2f} seconds")

            except Exception as e:
                print(f"\nERROR occurred during analysis of {gene}: {e}")
                print(traceback.format_exc())
                cleanup()

    except KeyboardInterrupt:
        print("\nReceived KeyboardInterrupt, performing cleanup.")
        cleanup()
    except Exception as e:
        print(f"\nERROR occurred in main: {e}")
        print(traceback.format_exc())
        helpinfo.PrintHelp()
        cleanup()
    finally:
        if d_results:
            util.PrintCitation(d_results)

        util.log_memory_usage("after final cleanup")

        end_main = time.perf_counter()
        duration_main = end_main - start_main
       
        if production_logger:
            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
            production_logger.info(prepend + "RECUR run completed\n", extra={'to_file': True, 'to_console': True})
            production_logger.info(f"*** RECUR finishes in {duration_main:.2f} seconds ***", extra={'to_file': True, 'to_console': True})
            # Flush and close the handlers to ensure all logs are written
            for handler in production_logger.handlers:
                handler.flush()
                handler.close()
        else:
            print("RECUR run completed\n")
            print(f"*** RECUR finishes in {duration_main:.2f} seconds ***")

        cleanup()
        kill_child_processes(os.getpid())
        
if __name__ == "__main__":
    main() 
