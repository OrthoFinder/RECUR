from __future__ import absolute_import

import datetime
import gc
import logging
import multiprocessing as mp
import os
import platform
import random
import shutil
import signal
import subprocess
import sys
import time
import traceback
import warnings
from collections import Counter
# from concurrent.futures import ProcessPoolExecutor, as_completed
# from functools import partial
from typing import Dict, List, Optional, Tuple, Union

import dendropy
import numpy as np
import psutil
from rich import print, progress

from recur import __version__, __location__, helpinfo
from recur.run import run_commands
from recur.utils import files, process_args, util, parallel_task_manager
from recur.utils import analytic_tools as at


# warnings.filterwarnings("ignore", module='dendropy')


def CanRunCommand(
        command: str, 
        env: Optional[Dict[str, str]] = None, 
        print_info: bool = False,
        iqtree_version: str = "iqtree2"
    ) -> bool:
    try:
        process = subprocess.Popen(
            command, 
            env=env, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True
        )
        stdout, stderr = process.communicate()
        if print_info:
            if stdout:
                if command in [
                    f"command -v {iqtree_version}", 
                    f"which {iqtree_version}", 
                    f"where {iqtree_version}", 
                    f"Get-Command {iqtree_version}"
                    ]:
                    print(f"Using {iqtree_version} from: ")
                print(stdout.strip())
            if stderr:
                print(stderr.strip())

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


def setup_environment() -> Dict[str, str]:
    os.environ["OPENBLAS_NUM_THREADS"] = "1"  # Single-threaded numpy/openblas

    my_env = os.environ.copy()
    local_bin_dir = os.path.join(__location__, 'bin')  # Make sure __location__ is defined

    if os.name == "nt":
        # Windows-specific directories
        bin_dirs = [
            os.path.expanduser(r"~/bin"),
            os.path.expanduser(r"~\AppData\Local\Programs\Python\Scripts"),
            r"C:\Program Files\SomeExecutableDir",
        ]
    else:
        # Unix-like directories
        bin_dirs = [
            local_bin_dir,
            "/opt/bin",
            "/usr/bin",
            "/usr/local/bin",
            os.path.expanduser("~/bin"),
            os.path.expanduser("~/.local/bin"),
            os.path.expanduser("~/local/bin"),
        ]
        
    for bin_dir in bin_dirs:
        my_env["PATH"] = bin_dir + os.pathsep + my_env["PATH"]

    conda_prefix = my_env.get("CONDA_PREFIX")
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, "Scripts") if os.name == "nt" else os.path.join(conda_prefix, "bin")
        my_env["PATH"] = conda_bin + os.pathsep + my_env["PATH"]
    
    return my_env

def initialise_recur(show_iqtree_path: bool = False, iqtree_version: str = "iqtree2") -> Dict[str, str]:
    my_env = setup_environment()
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

    iqtree_path = shutil.which(iqtree_version, path=my_env["PATH"])
    iqtree_version_cmd = f'"{iqtree_path}" --version'
    if CanRunCommand(iqtree_version_cmd, env=my_env, iqtree_version=iqtree_version):
        if show_iqtree_path:
            print(f"\ncan run {iqtree_version} - [bold green]ok[/bold green]")
            if system == "Windows":
                CanRunCommand(f"where {iqtree_version}", env=my_env, print_info=True, iqtree_version=iqtree_version)
            else:
                CanRunCommand(f"which {iqtree_version}", env=my_env, print_info=True, iqtree_version=iqtree_version)
        return my_env

    print("Cannot proceed. IQ-TREE does not exist in either local bin or system-wide PATH.")
    print("Please ensure IQ-TREE is properly installed before running RECUR!\n")
    sys.exit(1)

def signal_handler(sig: int, frame) -> None:
    print(f"Received signal {sig}")
    gc.collect()
    sys.exit(0)

def setup_signal_handler() -> None:
    if os.name != 'nt':
        signals = [signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGUSR1, signal.SIGUSR2]
    else:
        signals = [signal.SIGTERM, signal.SIGINT]

    for sig in signals:
        signal.signal(sig, signal_handler)

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
        os._exit(1)

    children = parent.children(recursive=True)
    for process in children:
        for _ in range(retries):
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

def ParentChildRelation(treefile: str,
                        outgroup_squences: List[str],
                        n_species: int,
                        preserve_underscores: bool
                        ) -> Tuple[Optional[str], List[str], List[str], List[str], str]:
    try:

        with open(treefile, 'r') as f:
            t = dendropy.Tree.get(file=f,
                                  schema="newick",
                                  preserve_underscores=preserve_underscores,
                                  case_sensitive_taxon_labels=False)

        t.is_rooted = True
        if len(outgroup_squences) == 1:
            outgroup_mrca = t.find_node_with_taxon_label(outgroup_squences[0])
        else:
            outgroup_mrca = t.mrca(taxon_labels=outgroup_squences)

        t.reroot_at_edge(outgroup_mrca.edge, update_bipartitions=False)
        rt = t.seed_node
        rt_children = rt.child_nodes()

        root_of_interest = next((child for child in rt_children if child != outgroup_mrca), None)

        if root_of_interest is None:
            return None, outgroup_squences, [], [], "No root of interest found; tree structure might be incorrect."

        root_node = root_of_interest.label
        root_node = root_node if "/" not in root_node else root_node.split("/")[0]

        outgroup_subtree_species = [leaf.taxon.label for leaf in outgroup_mrca.leaf_iter()]

        branch_count = 0
        parent_list: List[str] = []
        child_list: List[str] = []
        for nd in root_of_interest.postorder_iter():

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

                parent_list.append(parent)
                child_list.append(child)

                branch_count += 1

        if len(outgroup_subtree_species) != len(outgroup_squences):
            outgroup_squences = outgroup_subtree_species
        
        # Bifurcating check
        expected_relationships = 2 * (n_species - len(outgroup_squences)) - 2

        if branch_count != expected_relationships:
            return root_node, outgroup_squences, parent_list, child_list, \
                f"Taxon count {branch_count} does not match expected relationships {expected_relationships}."

        return root_node, outgroup_squences, parent_list, child_list, ""
    except BrokenPipeError:
        print("Broken pipe error.")
        return root_node, outgroup_squences, parent_list, child_list, "Broken pipe error"
    except Exception as e:
        error_msg = f"ERROR in ParentChildRelation: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return root_node, outgroup_squences, parent_list, child_list, error_msg


def count_mutations(parent_list: List[str],
                    child_list: List[str],
                    sequence_dict: Dict[str, str],
                    residue_dict: Dict[str, int],
                    dash_exist: bool = False,
                    binary_sequence_dict: Optional[Dict[str, str]] = None
                    ) -> Counter[Tuple[int, int, int]]:

    parent_num_list = [
        [residue_dict.get(res, 0) for res in sequence_dict[parent]]
        for parent in parent_list
    ]
    child_num_list = [
        [residue_dict.get(res, 0) for res in sequence_dict[child]]
        for child in child_list
    ]
    parent_array = np.array(parent_num_list)
    child_array = np.array(child_num_list)

    del sequence_dict, parent_num_list, child_num_list

    if dash_exist:
        binary_parent_num_list = [
            [1 if res == "1" else 0 for res in binary_sequence_dict[parent]]
            for parent in parent_list
        ]
        binary_child_num_list = [
            [1 if res == "1" else 0 for res in binary_sequence_dict[child]]
            for child in child_list
        ]

        binary_parent_array = np.array(binary_parent_num_list)
        binary_child_array = np.array(binary_child_num_list)

        parent_array = parent_array * binary_parent_array
        child_array = child_array * binary_child_array

        del binary_parent_num_list, binary_child_num_list, \
            binary_parent_array, binary_child_array

    parent_child_diff = parent_array != child_array
    row_indices, col_indices = np.where(parent_child_diff)

    parent_res_id = parent_array[row_indices, col_indices]
    child_res_id = child_array[row_indices, col_indices]

    # parent_mask = np.isin(parent_res_id, util.reserved_chars_index, invert=True)
    parent_mask = np.where(parent_res_id != 0)
    parent_res_id = parent_res_id[parent_mask]
    child_res_id = child_res_id[parent_mask]
    col_indices = col_indices[parent_mask]

    # child_mask = np.isin(child_res_id, util.reserved_chars_index, invert=True)
    child_mask = np.where(child_res_id != 0)

    parent_res_id = parent_res_id[child_mask]
    child_res_id = child_res_id[child_mask]
    col_indices = col_indices[child_mask]

    parent_child_tuples = [*zip(col_indices, parent_res_id, child_res_id)]
    rec_loc_count_dict = Counter(parent_child_tuples)

    del parent_child_diff, row_indices, col_indices, child_mask, \
        parent_mask, parent_res_id, child_res_id, parent_child_tuples

    return rec_loc_count_dict, parent_array, child_array

def get_recurrence_list(rec_loc_count_dict: Counter[Tuple[int, int, int]],
                        residue_dict_flip: Dict[int, str],
                        ) -> List[List[Union[int, str, float]]]:

    rec_loc_count_dict2 = {}
    flipflop_dict = {}
    for key, val in rec_loc_count_dict.most_common():
        if val < 2:
            break
        rec_loc_count_dict2[key] = val
        res_loc, parent_id, child_id = key
        flipflop_key = (res_loc, child_id, parent_id)
        flipflop = rec_loc_count_dict.get(flipflop_key)
        if flipflop is not None:
            flipflop_dict[key] = flipflop

    recurrence_list: List[List[Union[int, str, float]]] = [
        [
            int(res_loc),
            str(residue_dict_flip[parent_id]),
            str(residue_dict_flip[child_id]),
            int(recurrence),
            int(flipflop_dict.get((res_loc, parent_id, child_id), 0))
         ]
        for (res_loc, parent_id, child_id), recurrence in rec_loc_count_dict2.items()
    ]

    return recurrence_list

# def WorkerProcessAndCount(file: str,
#                           mcs_alnDir: str,
#                           parent_list: List[str],
#                           child_list: List[str],
#                           isnuc_fasta: bool,
#                           sequence_type: str,
#                           residue_dict: Dict[str, int],
#                           res_loc_list: List[int],
#                           production_logger: logging.Logger,
#                           dash_exist: bool = False,
#                           binary_sequence_dict: Dict[str, str] = {},
#                           ) -> Tuple[Dict[Tuple[int, int, int], int], str]:

#     rec_loc_count_dict: Dict[Tuple[int, int, int], int] = {}
#     error_msg = ""

#     try:
#         file_path = os.path.join(mcs_alnDir, file)
#         mcs_combined_prot_seqs_dict, _, _ = files.FileReader.ReadAlignment(file_path)

#         if isnuc_fasta:
#             mcs_combined_prot_seqs_dict, _ = util.GetSeqsDict(mcs_combined_prot_seqs_dict, sequence_type)

#         parent_num_list = [
#             [residue_dict.get(res, 0) for res in mcs_combined_prot_seqs_dict[parent]]
#             for parent in parent_list
#         ]
#         child_num_list = [
#             [residue_dict.get(res, 0) for res in mcs_combined_prot_seqs_dict[child]]
#             for child in child_list
#         ]

#         parent_array = np.array(parent_num_list)
#         child_array = np.array(child_num_list)

#         del mcs_combined_prot_seqs_dict, parent_num_list, child_num_list

#         if dash_exist:
#             binary_parent_num_list = [
#                 [1 if res == "1" else 0 for res in binary_sequence_dict[parent]]
#                 for parent in parent_list
#             ]
#             binary_child_num_list = [
#                 [1 if res == "1" else 0 for res in binary_sequence_dict[child]]
#                 for child in child_list
#             ]

#             binary_parent_array = np.array(binary_parent_num_list)
#             binary_child_array = np.array(binary_child_num_list)

#             parent_array = parent_array * binary_parent_array
#             child_array = child_array * binary_child_array

#             del binary_parent_num_list, binary_child_num_list, \
#                 binary_parent_array, binary_child_array

#         parent_child_diff = parent_array != child_array
#         row_indices, col_indices = np.where(parent_child_diff)

#         mask = np.isin(col_indices, res_loc_list)
#         col_idx = col_indices[mask]
#         row_idx = row_indices[mask]

#         parent_res_id = parent_array[row_idx, col_idx]
#         child_res_id = child_array[row_idx, col_idx]

#         # parent_mask = np.isin(parent_res_id, util.reserved_chars_index, invert=True)
#         parent_mask = np.where(parent_res_id != 0)
#         parent_res_id = parent_res_id[parent_mask]
#         child_res_id = child_res_id[parent_mask]
#         col_idx = col_idx[parent_mask]

#         # child_mask = np.isin(child_res_id, util.reserved_chars_index, invert=True)
#         child_mask = np.where(child_res_id != 0)
#         parent_res_id = parent_res_id[child_mask]
#         child_res_id = child_res_id[child_mask]
#         col_idx = col_idx[child_mask]

#         parent_child_tuples = [*zip(col_idx, parent_res_id, child_res_id)]
#         rec_loc_count_dict = Counter(parent_child_tuples)
#         del parent_child_diff, mask, row_idx, col_idx, row_indices, col_indices, \
#             parent_res_id, child_res_id, parent_child_tuples, parent_array, child_array

#     except BrokenPipeError:
#         error_msg = "Broken pipe error while processing file."
#         print(error_msg)

#     except Exception as e:
#         error_msg = f"ERROR in WorkerProcessAndCount for mcs task {file}: {e}"
#         print(error_msg)
#         print(traceback.format_exc())
#         production_logger.error(error_msg)

#     return rec_loc_count_dict, error_msg


# def process_mcs_files_in_chunks(mcs_alnDir: str,
#                                 parent_list: List[str],
#                                 child_list: List[str],
#                                 residue_dict: Dict[str, int],
#                                 nthreads: int,
#                                 isnuc_fasta: bool,
#                                 sequence_type: str,
#                                 res_loc_list: List[int],
#                                 production_logger: logging.Logger,
#                                 window_width: int,
#                                 dash_exist: bool = False,
#                                 binary_sequence_dict: Optional[Dict[str, str]] = None,
#                                 update_cycle: Optional[int] = None,
#                                 mcs_batch_size: Optional[int] = None
#                                 ) -> List[Dict[Tuple[int, int, int], int]]:

#     mcs_files = os.listdir(mcs_alnDir)
#     total_file_count = len(mcs_files)

#     results = []
#     worker = partial(WorkerProcessAndCount,
#                      mcs_alnDir=mcs_alnDir,
#                      parent_list=parent_list,
#                      child_list=child_list,
#                      isnuc_fasta=isnuc_fasta,
#                      sequence_type=sequence_type,
#                      residue_dict=residue_dict,
#                      res_loc_list=res_loc_list,
#                      production_logger=production_logger,
#                      dash_exist=dash_exist,
#                      binary_sequence_dict=binary_sequence_dict,
#                      )
#     if mcs_batch_size is not None:
#         batches = [mcs_files[i:i + mcs_batch_size] for i in range(0, total_file_count, mcs_batch_size)]

#     if update_cycle is not None:
#         mcs_progress = progress.Progress(
#         progress.TextColumn("[progress.description]{task.description}"),
#         progress.BarColumn(bar_width=window_width // 2),
#         progress.SpinnerColumn(),
#         progress.MofNCompleteColumn(),
#         progress.TimeElapsedColumn(),
#         transient=False,
#         # progress.TextColumn("{task.completed}/{task.total}")
#         )
#         task = mcs_progress.add_task("[magenta]Processing...", total=total_file_count)

#         mcs_progress.start()
#     try:
#         with ProcessPoolExecutor(max_workers=nthreads) as executor:

#             if mcs_batch_size is not None:
#                 futures = [executor.submit(worker, file_data) for batch in batches for file_data in batch]
#             else:
#                 futures = [executor.submit(worker, file_data) for file_data in mcs_files]


#             for i, future in enumerate(as_completed(futures)):
#                 try:
#                     result, error_msg = future.result()
#                     if error_msg:
#                         print(error_msg)
#                         production_logger.error(error_msg)
#                         break

#                     if result:
#                         results.append(result)

#                 except Exception as e:
#                     error_msg = f"ERROR during processing: {e}"
#                     print(error_msg)
#                     print(traceback.format_exc())
#                     production_logger.error(error_msg)
#                     break

#                 finally:
#                     if update_cycle is not None:
#                         if (i + 1) % update_cycle == 0:
#                             mcs_progress.update(task, advance=update_cycle)

#             if update_cycle is not None:
#                 mcs_progress.stop()

#     except Exception as e:
#         error_msg = f"ERROR during processing files: {e}"
#         print(error_msg)
#         print(traceback.format_exc())
#         production_logger.error(error_msg)

#     return results

def mcs_count_greater(
        mcs_results: List[Dict[Tuple[int, int, int], int]],
        recurrence_list: List[List[Union[str, int, float]]],
        residue_dict: Dict[str, int],
    ) -> List[int]:
    try:
        count_greater_list = []
        for rec_list in recurrence_list:
            rec_loc = int(rec_list[0])
            parent = str(rec_list[1])
            child = str(rec_list[2])
            recurrence = int(rec_list[3])

            count_greater = 0
            for mcs_count_dict in mcs_results:
                mcs_rec = mcs_count_dict.get((rec_loc, residue_dict[parent], residue_dict[child]))
                if mcs_rec is not None:
                    if mcs_rec >= recurrence:
                        count_greater += 1
            
            count_greater_list.append(count_greater)

    except Exception as e:
        error_msg = f"ERROR during compute_p_values: {e}"
        print(error_msg)
        print(traceback.format_exc())
    finally:
        return count_greater_list
    

def update_recurrence_list(
        R: List[int],
        B: int,
        res_loc_count_dict: Dict[Tuple[int, int, int], int],
        recurrence_list: List[List[Union[str, int, float]]],
        combined_prot_seqs_dict: Dict[str, str],
        species_of_interest: List[str],
        residue_dict_flip: Dict[int, str],
        protein_len: int,
        alpha: float = 0.05,
        q: float = 0.05,
        method: Optional[str] = "fdr_bh",
        pval_stats: bool = False,
    ) -> List[List[Union[str, int, float]]]:

    extant_seq = {species: seq for species, seq in combined_prot_seqs_dict.items() if species in species_of_interest}
    ident_dict = {}
    
    p_hat, p_adj, ci_lo_adj, ci_hi_adj, decision_sig, robust_sig = \
        at.sitewise_decision(R, B, alpha=alpha, q=q, method=method)

    res_loc_info_dict= util.get_sorted_res_loc_info(res_loc_count_dict, protein_len)
    for rec_loc, res in enumerate(zip(*extant_seq.values())):
        ident_dict[rec_loc] = res
    
    precsion = len(str(B))
    print()
    for i, rec_list in enumerate(recurrence_list):
        res_loc = int(rec_list[0])
        parent_child = []
        counts = []

        rec_list.append(np.round(p_hat[i], precsion))
        rec_list.append(np.round(p_adj[i], precsion))

        if pval_stats:
            lower_ci = np.round(ci_lo_adj[i], precsion)
            upper_ci = np.round(ci_hi_adj[i], precsion)
            rec_list.append(lower_ci)
            rec_list.append(upper_ci)
            decision = "Reject" if decision_sig[i] == 1 else "Accept"
            rec_list.append(decision)
            robust = "Yes" if robust_sig[i] == 1 else "No"
            rec_list.append(robust)

        for parent_id, child_id, recurrence in res_loc_info_dict[res_loc]:
            parent = residue_dict_flip[parent_id]
            child = residue_dict_flip[child_id]
            parent_child.append(">".join((parent, child)))
            counts.append(recurrence)

        data_str_list = [*(map(str, counts))]
        parent_child_data = [*zip(parent_child, data_str_list)]
        parent_child_data_str = ",".join([":".join(pcd) for pcd in parent_child_data])
        rec_list.append(parent_child_data_str)

        res_freq = [*Counter(ident_dict[res_loc]).items()]
        res_freq = sorted(res_freq, reverse=True, key=lambda x: x[1])
        res_freq_str = ",".join([":".join((res, str(freq))) for res, freq in res_freq])
        rec_list.append(res_freq_str)

    recurrence_list.sort(key=lambda x: (-float(x[6]), float(x[3])), reverse=True)

    return recurrence_list

def main(args: Optional[List[str]] = None):

    d_results = None
    production_logger = None
    setup_signal_handler()

    if not args:
        args = sys.argv[1:]

    if not args or len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.PrintHelp()
        sys.exit()

    elif args[0] == "-v" or args[0] == "--version":
        print(f"[dark_goldenrod]RECUR[/dark_goldenrod]:v[deep_sky_blue2]{__version__}[/deep_sky_blue2]")
        sys.exit()

    start_main = time.perf_counter()

    try:
        input_command = "RECUR command: " + " ".join(sys.argv)  + "\n"

        options, alnDir, alnPath, resultsDir_nonDefault = process_args.ProcessArgs(args)
        
        # iqtree_version = options.iqtree_version if options.iqtree_version else "system"
        my_env = initialise_recur(options.show_iqtree_path, iqtree_version=options.iqtree_version)

        if options.system_info:
            util.get_system_info()

        aln_path_dict = files.FileHandler.ProcessesNewAln(alnDir, alnPath)

        aln_len = len(aln_path_dict)
        if isinstance(options.gene_tree, dict):
            if len(options.gene_tree) <= aln_len and len(options.gene_tree) == 1:
                options.gene_tree = {gene: [*options.gene_tree.values()][0] for gene in aln_path_dict}

        if isinstance(options.outgroups, dict):
            if len(options.outgroups) <= aln_len and len(options.outgroups) == 1:
                options.outgroups = {gene: [*options.outgroups.values()][0] for gene in aln_path_dict}

        residue_dict, residue_dict_flip = util.residue_table()

        if options.compute_recurrence:
            count = 0
            for gene, aln_path in aln_path_dict.items():

                asr = True
                fix_branch_length = True
                try:
                    width = os.get_terminal_size().columns
                except OSError as e:
                    width = 80

                if count > 0:
                    util.print_centered_text(width, f"Processing gene {gene}")
                count += 1
                alnFN = os.path.basename(aln_path)
                base_dir = os.path.join(options.project_dir, f"{alnFN}.recur")

                filehandler = files.FileHandler(base_dir)
                filereader = files.FileReader()
                filewriter = files.FileWriter()
                filehandler.gene_of_interest = alnFN
                mcs_results, B = filereader.ReadMCSRecurrenceCount(base_dir)
                recurrence_count_file = filehandler.GetRecurrenceCountPhylogenyFN(options.project_dir)
                rec_loc_count_dict = filereader.ReadRecurrenceCount(recurrence_count_file)
                combined_prot_seqs_fn = os.path.join(base_dir, filehandler.GetCombinedProtSeqsFN())

                combined_prot_seqs_dict, protein_len, _ = filereader.ReadAlignment(combined_prot_seqs_fn)
                alignment_dict, _, _ = filereader.ReadAlignment(aln_path)
                species_of_interest = [*alignment_dict.keys()]
                
                recurrence_list_fn = filehandler.GetRecurrenceListRealPhylogenyFN(options.project_dir)
                recurrence_list = filereader.ReadRecurrenceList(recurrence_list_fn)
                M = len(recurrence_list)
                print()
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                print(prepend + "Starting compute p values.")
                R = mcs_count_greater(
                    mcs_results,
                    recurrence_list,
                    residue_dict,
                )
                if options.pval_adjust_method is None:
                    options.pval_adjust_method = at.method_selection(M, options.site_dependence)

                recurrence_list_updated = update_recurrence_list(
                    R, 
                    B,
                    rec_loc_count_dict,
                    recurrence_list,
                    combined_prot_seqs_dict,
                    species_of_interest,
                    residue_dict_flip,
                    protein_len,
                    alpha=options.significance_level, # needs be provided by the user, if they don't use the default setting
                    q=options.fdr_level, # needs be provided by the user, if they don't use the default setting
                    method=options.pval_adjust_method, # needs be provided by the user, if they don't use the default setting
                    pval_stats=options.pval_stats
                )

                filewriter.WriteRecurrenceList(
                    recurrence_list_updated, 
                    filehandler.GetRecurrenceListFN(options.project_dir),
                    options
                    )
    
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                print(prepend + "p values computing complete.")

                os.remove(combined_prot_seqs_fn)
                os.remove(recurrence_count_file)
                os.remove(recurrence_list_fn)

                d_results = os.path.normpath(base_dir) + os.path.sep
                rec_results = os.path.normpath(filehandler.GetRecurrenceListFN(options.project_dir))
                print("\nResults:\n    %s\n" % rec_results)

        else:
                
            count = 0
            for gene, aln_path in aln_path_dict.items():

                asr = True
                fix_branch_length = True

                try:
                    start = time.perf_counter()

                    try:
                        width = os.get_terminal_size().columns
                    except OSError as e:
                        width = 80

                    if count > 0:
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
                    base_dir = os.path.join(resultsDir_nonDefault, f"{alnFN}.recur") \
                        if resultsDir_nonDefault else os.path.join(alnDir, f"{alnFN}.recur")
                    # if not os.path.exists(base_dir):
                    #     os.mkdir(base_dir)

                    filehandler.CreateOutputDirectories(options, base_dir)
                    filehandler.CreateMCSDirectories(options)
                    results_dir = filehandler.GetResultsDirectory()
                    production_logger = util.setup_logging(results_dir, "w", "brief")

                    production_logger.info(f"{input_command}", extra={'to_file': True, 'to_console': False})
                    print()
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Starting RECUR v%s" % __version__, extra={'to_file': True, 'to_console': True})

                    alignment_dict, alignment_len, dash_exist = filereader.ReadAlignment(aln_path)
                    n_species = len(alignment_dict)

                    production_logger.info(f"Analysing: {os.path.basename(aln_path)}", extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Number of sequences found: {n_species}", extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Length of alignment: {alignment_len}", extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Results Directory: {results_dir}\n", extra={'to_file': True, 'to_console': True})

                    species_of_interest = [*alignment_dict.keys()]
                    isnuc = util.CheckSequenceType([*alignment_dict.values()])

                    if isnuc and options.sequence_type == "AA":
                        print("\nWARNING: The input sequence type does not match the existing alignments. RECUR will convert the codon alignments into protein alignments with the default CODON model\n")
                        options.sequence_type = "CODON1"

                    check_exist = 0
                    real_phyDir = filehandler.GetRealPhylogenyDir()
                    exist_treefile = False
                    exist_iqtree = False
                    exist_state = False

                    for file in util.iter_dir(real_phyDir):
                        file_extension = file.rsplit(".", 1)[1]
                        if file_extension == "state":
                            exist_state = True
                            check_exist += 1

                        elif file_extension == "treefile":
                            exist_treefile = True
                            check_exist += 1

                        elif file_extension == "iqtree":
                            exist_iqtree = True
                            check_exist += 1

                    restart_step1 = False
                    restart_step2 = False
                    restart_step3 = False
                    override = options.override

                    gene_tree = options.gene_tree[gene] if isinstance(options.gene_tree, dict) else None
                    if gene_tree is not None:
                        step1_info = "Step1: Inferring phylogenetic tree and model of evolution"
                        if options.restart_from == 1:
                            restart_step1 = True
                            restart_step3 = True
                            override = False
                        elif options.restart_from == 2:
                            restart_step3 = True
                            override = False
                    else:
                        step1_info = "Step1: Inferring ancestral sequences, phylogenetic tree and model of evolution"
                        if options.restart_from == 1:
                            restart_step1 = True
                            restart_step2 = True
                            restart_step3 = True
                            override = False
                        elif options.restart_from == 2:
                            restart_step2 = True
                            restart_step3 = True
                            override = False
                        elif options.restart_from == 3:
                            restart_step3 = True
                            override = False

                    production_logger.info(step1_info, extra={'to_file': True, 'to_console': True})
                    production_logger.info("="*len(step1_info), extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Results Directory: {real_phyDir}\n", extra={'to_file': True, 'to_console': True})

                    if (options.usr_state and options.usr_tree and options.usr_iqtree):
                        statefile = options.usr_state
                        treefile = options.usr_tree
                        iqtreefile = options.usr_iqtree
                        if gene_tree is None:
                            production_logger.info("NOTE: with the provided treefile, RECUR will skip Step1.\n", extra={'to_file': True, 'to_console': True})
                        else:
                            production_logger.info("NOTE: with the provided statefile and treefile, RECUR will skip Step1.\n", extra={'to_file': True, 'to_console': True})

                    elif (check_exist == 3 and not override) and not restart_step1:
                        if gene_tree is not None:
                            statefile = filehandler.GetStateFileFN()
                            treefile = filehandler.GetTreeFileFN()
                            iqtreefile = filehandler.GetIQTreeFileFN()
                            production_logger.info("NOTE: with the existing statefile, treefile and iqtreefile, RECUR will skip Step1.\n", extra={'to_file': True, 'to_console': True})
                        else:
                            treefile = filehandler.GetTreeFileFN()
                            iqtreefile = filehandler.GetIQTreeFileFN()
                            production_logger.info("NOTE: with the existing treefile and iqtreefile, RECUR will skip Step1.\n", extra={'to_file': True, 'to_console': True})
                    else:
                        if check_exist == 2:
                            if exist_iqtree and exist_treefile and not exist_state:
                                if gene_tree is not None:
                                    production_logger.info("NOTE: RECUR is forced to restart from Step1 due to the missing statefile\n", extra={'to_file': True, 'to_console': True})
                                else:
                                    restart_step2 = True

                            elif exist_iqtree and not exist_treefile and exist_state:
                                production_logger.info("NOTE: RECUR is forced to restart from Step1 due to the missing treefile\n", extra={'to_file': True, 'to_console': True})
                                if gene_tree is None:
                                    restart_step2 = True

                            elif not exist_iqtree and exist_treefile and exist_state:
                                production_logger.info("NOTE: RECUR is forced to restart from Step1 due to the missing iqtreefile\n", extra={'to_file': True, 'to_console': True})
                                if gene_tree is None:
                                    restart_step2 = True

                        elif check_exist > 0 and check_exist < 2:
                            production_logger.info("NOTE: RECUR is forced to restart from Step1 due to some missing files\n", extra={'to_file': True, 'to_console': False})
                            if gene_tree is None:
                                restart_step2 = True

                        if restart_step1 or not override:
                            production_logger.info("### Restart RECUR from Step1 ###\n", extra={'to_file': True, 'to_console': True})

                        if override or restart_step1:
                            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                            production_logger.info(prepend + "Running IQ-TREE to build the real phylogeny", extra={'to_file': True, 'to_console': True})
                            production_logger.info("Using %d RECUR thread(s), %d IQ-TREE thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads), extra={'to_file': True, 'to_console': True})

                        if gene_tree is None:
                            asr = False
                            fix_branch_length = False
                        else:
                            fix_branch_length = options.fix_branch_length

                        commands = run_commands.GetGeneTreeBuildCommands(
                            [aln_path],
                            real_phyDir,
                            options.evolution_model,
                            options.iqtree_nthreads,
                            phy_seed=options.seed,
                            asr=asr,
                            sequence_type=options.sequence_type,
                            gene_tree=gene_tree,
                            bootstrap=options.bootstrap,
                            sh_alrt=options.bootstrap,
                            branch_test=options.branch_test,
                            fix_branch_length=fix_branch_length,
                            iqtree_version=options.iqtree_version,
                            iqtree_cmd_dict=options.iqtree_cmd_dict
                        )

                        if gene_tree is None:
                            production_logger.info(f"step1 {options.iqtree_version} gene tree building command: ", extra={'to_file': True, 'to_console': False})
                            production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                        else:
                            production_logger.info(f"step1 {options.iqtree_version} ancestral state reconstruction command: ", extra={'to_file': True, 'to_console': False})
                            production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                        run_commands.RunCommand(
                            commands,
                            real_phyDir,
                            env=my_env,
                            nthreads=options.recur_nthreads,
                            delete_files=True,
                            files_to_keep=["state", "treefile", "iqtree"],
                            fd_limit=options.fd_limit
                        )

                        treefile = filehandler.GetTreeFileFN()
                        iqtreefile = filehandler.GetIQTreeFileFN()

                    best_evolution_model = filereader.ReadIQTreeFile(iqtreefile)

                    if gene_tree is None:
                        if override or restart_step1:
                            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                            production_logger.info(prepend + f"Tree inference complete, best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file':True, 'to_console': True})
                        else:
                            production_logger.info(f"Best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file':True, 'to_console': True})

                        step2_info = f"Step2: Inferring ancestral sequences"
                        production_logger.info(step2_info, extra={'to_file': True, 'to_console': True})
                        production_logger.info("="*len(step2_info), extra={'to_file': True, 'to_console': True})
                        step2_results_info = f"Results Directory: {filehandler.GetRealPhylogenyDir()}\n"
                        production_logger.info(step2_results_info, extra={'to_file': True, 'to_console': True})

                        if not restart_step2 and exist_state and not override:
                            statefile = filehandler.GetStateFileFN()
                            production_logger.info("NOTE: with the existing statefile, RECUR will skip Step2.\n", extra={'to_file': True, 'to_console': True})

                        else:
                            if not restart_step2 and (exist_iqtree and exist_treefile and not exist_state):
                                production_logger.info("NOTE: RECUR is forced to restart from Step2 due to the missing statefile\n", extra={'to_file': True, 'to_console': True})

                            elif not restart_step1 and restart_step2:
                                production_logger.info("### Restart RECUR from Step2 ###\n", extra={'to_file': True, 'to_console': True})

                            filehandler.UpdateTreeFile(treefile)

                            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                            production_logger.info(prepend + f"Starting ancestral state reconstruction.", extra={'to_file': True, 'to_console': True})

                            commands = run_commands.GetGeneTreeBuildCommands(
                                [aln_path],
                                real_phyDir,
                                best_evolution_model,
                                options.iqtree_nthreads,
                                phy_seed=options.seed,
                                sequence_type=options.sequence_type,
                                gene_tree=treefile,
                                asr=True,
                                fix_branch_length=options.fix_branch_length,
                                iqtree_version=options.iqtree_version,
                                iqtree_cmd_dict=options.iqtree_cmd_dict
                            )

                            production_logger.info(f"step2 {options.iqtree_version} ancestral state reconstruction command: ",  extra={'to_file': True, 'to_console': False})
                            production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                            run_commands.RunCommand(
                                commands,
                                real_phyDir,
                                env=my_env,
                                nthreads=options.recur_nthreads,
                                delete_files=True,
                                files_to_keep=["state", "treefile", "iqtree"],
                                fd_limit=options.fd_limit
                            )

                            statefile = filehandler.GetStateFileFN()
                            if not restart_step2 and restart_step3:
                                production_logger.info("Ancestral state reconstruction complete\n", extra={'to_file': True, 'to_console': True})
                            else:
                                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                                production_logger.info(prepend + "Ancestral state reconstruction complete\n", extra={'to_file': True, 'to_console': True})
                    else:
                        statefile = filehandler.GetStateFileFN()
                        if not restart_step1 and restart_step3:
                            production_logger.info(f"Best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file': True, 'to_console': True})
                        else:
                            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                            production_logger.info(prepend + f"Ancestral state reconstruction complete, best fitting model of sequence evolution: {best_evolution_model}\n", extra={'to_file': True, 'to_console': True})

                    filehandler.CheckFileCorrespondance(gene, statefile, treefile)
                    node_seq_dict = filereader.ReadStateFile(statefile)

                    combined_seq_dict = {k: v for d in (node_seq_dict, alignment_dict) for k, v in d.items()}

                    outgroup_mrca, preserve_underscores = util.CheckOutgroups(outgroup_mrca, alignment_dict)

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
                        alignment_dict, _ = util.GetSeqsDict(alignment_dict, options.sequence_type)
                    
                    if options.disk_save and options.multi_stage:
                        combined_prot_seqs_fn = filehandler.GetCombinedProtSeqsFN()
                        filewriter.WriteSeqsToAln(combined_prot_seqs_dict, combined_prot_seqs_fn)


                    # # ----------------- Binary phylogeny analysis --------------------------
                    if not dash_exist:
                        binary_combined_seq_dict: Dict[str, str] = {}
                    else:
                        prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                        production_logger.info(prepend + "Starting ancestral indel estimation.", extra={'to_file': True, 'to_console': True})

                        binary_alignment_dict = util.ConvertToBinary(alignment_dict)
                        binary_aln_path = filehandler.GetBinarySeqsFN()
                        filewriter.WriteSeqsToAln(binary_alignment_dict, binary_aln_path)

                        binary_tree_commands = run_commands.GetGeneTreeBuildCommands(
                            [binary_aln_path],
                            filehandler.GetBinaryPhylogenyDir(),
                            options.binary_model,
                            options.iqtree_nthreads,
                            phy_seed=options.seed,
                            gene_tree=treefile,
                            asr=True,
                            fix_branch_length=options.binary_blfix,
                            iqtree_version=options.iqtree_version,
                            iqtree_cmd_dict=options.iqtree_cmd_dict
                        )

                        production_logger.info(f"{options.iqtree_version} ancestral indel estimation command: ",  extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{binary_tree_commands[0]}\n", extra={'to_file': True, 'to_console': False})

                        run_commands.RunCommand(
                            binary_tree_commands,
                            filehandler.GetBinaryPhylogenyDir(),
                            env=my_env,
                            nthreads=options.recur_nthreads,
                            delete_files=True,
                            files_to_keep=["state", "treefile", "aln"],
                            fd_limit=options.fd_limit
                        )

                        binary_statefile = filehandler.GetBinaryStateFileFN()

                        binary_node_seq_dict = filereader.ReadStateFile(binary_statefile)
                        binary_combined_seq_dict = {k: v for d in (binary_node_seq_dict, binary_alignment_dict) for k, v in d.items()}

                        binary_node_seqs_fn = filehandler.GetBinaryNodeSeqsFN()
                        filewriter.WriteSeqsToAln(binary_node_seq_dict, binary_node_seqs_fn)

                        prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                        production_logger.info(prepend + "Ancestral indel estimation complete.\n", extra={'to_file': True, 'to_console': True})
                    # # ------------------------------------------------------------

                    alignment_dict.clear()
                    node_seq_dict.clear()
                    combined_seq_dict.clear()

                    root_node, outgroup_squences, parent_list, child_list, error_msg = \
                        ParentChildRelation(
                            treefile,
                            outgroup_mrca,
                            n_species,
                            preserve_underscores,
                        )

                    if error_msg:
                        production_logger.error(error_msg)
                        if not options.continue_on_error:
                            continue

                    if outgroup_mrca and len(outgroup_mrca) != len(outgroup_squences):
                        warnings.warn(f"Outgroup sequences provided not monophyletic. Outgroups will be updated. Please find the updated outgroups in the log file.")
                        production_logger.info("Updated outgroups: {outgroup_mrca}", extra={'to_file': True, 'to_console': False})
                        outgroup_mrca = outgroup_squences

                    rec_loc_count_dict, parent_arr, child_arr = count_mutations(
                        parent_list,
                        child_list,
                        combined_prot_seqs_dict,
                        residue_dict,
                        dash_exist=dash_exist,
                        binary_sequence_dict=binary_combined_seq_dict,
                    )

                    if dash_exist:
                        binary_modified_seq = util.ConvertToSequence(
                            parent_list,
                            child_list,
                            parent_arr,
                            child_arr,
                            residue_dict_flip
                        )
                        binary_modified_seq_fn = filehandler.GetBinaryCombinedProtSeqsFN()
                        filewriter.WriteSeqsToAln(binary_modified_seq, binary_modified_seq_fn)

                    del parent_arr, child_arr
                    production_logger.info(f"Root of subtree of interest (excluding outgroup sequences): {root_node}", extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Substitution matrix output: {filehandler.GetMutMatrixDir()}\n", extra={'to_file': True, 'to_console': True})
                    
                    filewriter.WriteMutMatrix(
                        rec_loc_count_dict,
                        residue_dict_flip,
                        protein_len,
                        filehandler.GetMutCountMatricesFN(),
                        filehandler.GetAccumMutCountMatricesFN()
                    )
                
                    recurrence_list = get_recurrence_list(
                        rec_loc_count_dict, 
                        residue_dict_flip
                    )
                    
                    # T_obs = np.asarray([item[3] for item in recurrence_list])
                    M = len(recurrence_list)

                    if not options.recDir:
                        recurrenceDir = alnDir
                    else:
                        recurrenceDir = options.recDir
                    
                    # if (options.disk_save and options.multi_stage) or options.just_recurrence:
                    filewriter.WriteRecurrenceListRealPhylogeny(
                        recurrence_list, 
                        filehandler.GetRecurrenceListRealPhylogenyFN(recurrenceDir)
                    )

                    filewriter.WriteRecurrenceCountToFile(
                        rec_loc_count_dict,
                        filehandler.GetRecurrenceCountPhylogenyFN(recurrenceDir)
                    )

                    if (options.disk_save and options.multi_stage) :
                        print("NOTE: You are running one the disk saving mode!")
                        print(
                            "If you don't wish to use `multi-stage-recur.sh`, "
                            "you need to run `recur -f project_dir -cr -pam pval_adjust_method` "
                            "to obtain the recurrence list.\n"
                        )

                    if len(recurrence_list) == 0:
                        production_logger.info(f"ATTENTION: No recurrence has identified for gene {gene}! Monte-Carlo Simiatlion will be SKIPPED!\n",
                                            extra={'to_file': True, 'to_console': True})
                        continue

                    res_loc_list = [int(res_list[0]) for res_list in recurrence_list]
                    
                    if not restart_step3:
                        if options.nalign is None:
                            mcp_method, B = at.min_mcs(
                                    M,
                                    method=options.pval_adjust_method,
                                    alpha=options.significance_level,
                                    q=options.fdr_level,
                                    rel_tol=options.relative_tolerance,
                                    extra_grid_cushion=options.grid_cushion,
                                    suspect_dependence=options.site_dependence,
                                    mc_error_control=options.mc_error_control
                                )
                            
                            options.nalign = B
                            options.pval_adjust_method = mcp_method
                        else:
                            options.pval_adjust_method = at.method_selection(M, options.site_dependence)
                        
                        # mcs_info = (
                        #     f"With {M} recurrence test, "
                        #     f"{options.nalign} Monte Carlo Simulations will be conducted based on {options.pval_adjust_method} "
                        #     f"(i.e., {at.METHODS_AVAILABLE[options.pval_adjust_method]}) adjustment method to adjust the p-values in multipletests.\n"
                        # )
                        # production_logger.info(mcs_info, extra={'to_file': False, 'to_console': True})
                        

                        plural_rec  = "test" if M == 1 else "tests"
                        plural_sims = "simulation" if options.nalign == 1 else "simulations"

                        mcs_info = (
                            "Monte-Carlo Recurrence Analysis Parameters\n"
                            "- Recurrence {label:<10} : {rec:,d} {plural_rec}\n"
                            "- Monte-Carlo           : {sims:,d} {plural_sims}\n"
                            "- P-value adjust method : {method_key} (i.e., {method_readable})\n"
                        ).format(
                            label="tests",
                            rec=M,
                            plural_rec=plural_rec,
                            sims=options.nalign,
                            plural_sims=plural_sims,
                            method_readable=at.METHODS_AVAILABLE[options.pval_adjust_method],
                            method_key=options.pval_adjust_method,
                        )
                        production_logger.info(
                            mcs_info,
                            extra={"to_file": False, "to_console": True},
                        )

                        msg = (
                            f"MCS Parameters | num_rec_tests={M} | num_mc_sims={options.nalign} | "
                            f"pval_adj_method={options.pval_adjust_method}"
                            "\n"
                        )
                        production_logger.info(
                            msg,
                            extra={
                                "to_file": True,
                                "to_console": False,
                                "num_rec_tests":   M,
                                "num_mc_sims":     options.nalign,
                                "pval_adj_method": options.pval_adjust_method,
                            },
                        )

                        if options.nalign > options.recur_limit and not options.just_recurrence:
                            warning_msg = (
                                f"\n"
                                f"  WARNING: Requested {options.nalign} alignments exceeds the current limit ({options.recur_limit}).\n"
                                f"\n"
                                f"   To save disk space and avoid performance issues, please use the multi-stage batch runner:\n"
                                f"     ./multi-stage-recur.sh -n {options.nalign} -d {aln_path} -b <batch_size> [optional: -m <pval_adjust_method>]\n"
                                f"\n"
                                f"   This will perform {options.nalign} Monte Carlo simulations "
                                f"using the p-value adjustment method: {options.pval_adjust_method}\n"
                                f"\n"
                                f"   If you still want to run everything in one stage, you can override the limit with:\n"
                                f"     --recur-limit {options.nalign} or -rl {options.nalign}\n"
                                f"\n"
                                f"  RECUR will now exit safely without running the simulation.\n"
                            )

                            production_logger.warning(
                                warning_msg,
                                extra={
                                    "to_file": True,
                                    "to_console": True,
                                    "num_rec_tests":   M,
                                    "num_mc_sims":     options.nalign,
                                    "pval_adj_method": options.pval_adjust_method,
                                },
                            )
                            sys.exit(0)

                    
                    if options.just_recurrence:
                        print("Done Computing Recurrence for Real Phylogeny.")
                        util.PrintCitation()
                        sys.exit(0)

                    if gene_tree is None:
                        step2_info = f"Step3: Simulating Sequence Evolution with {options.nalign} replicates"
                    else:
                        step2_info = f"Step2: Simulating Sequence Evolution with {options.nalign} replicates"

                    production_logger.info(step2_info, extra={'to_file': True, 'to_console': True})
                    production_logger.info("="*len(step2_info), extra={'to_file': True, 'to_console': True})
                    step2_results_info = f"Results Directory: {filehandler.mcs_dir}\n"
                    production_logger.info(step2_results_info, extra={'to_file': True, 'to_console': True})

                    mcs_faDir = filehandler.mcs_dir
                    if len(os.listdir(mcs_faDir)) != options.nalign or \
                        (restart_step1 or restart_step2 or restart_step3) or override:

                        if len(os.listdir(mcs_faDir)) != options.nalign and len(os.listdir(mcs_faDir)) > 0:
                            util.delete_files_in_directory(mcs_faDir)


                        msg = (
                            f"MCS Parameters | num_rec_tests={M} | num_mc_sims={options.nalign} | "
                            f"pval_adj_method={options.pval_adjust_method}"
                            "\n"
                        )
                        production_logger.info(
                            msg,
                            extra={
                                "to_file": True,
                                "to_console": False,
                                "num_rec_tests":   M,
                                "num_mc_sims":     options.nalign,
                                "pval_adj_method": options.pval_adjust_method,
                            },
                        )

                        identifier = "rooted_" + gene + "_alisim"
                        output_prefix = os.path.join(mcs_faDir, identifier)
                        if root_node is not None:
                            if options.sequence_type == "AA":
                                fn_root_node = ",".join((node_prot_seqs_fn, root_node))
                            else:
                                fn_root_node = ",".join((node_dna_seqs_fn, root_node))
                        else:
                            raise ValueError("Root node of interest is None.")
                        

                        mcs_commands = run_commands.GetMCsimulationCommand(
                            output_prefix,
                            options.iqtree_nthreads,
                            options.mcs_seed,
                            best_evolution_model,
                            treefile,
                            fn_root_node,
                            options.nalign,
                            iqtree_version=options.iqtree_version,
                            iqtree_cmd_dict=options.iqtree_cmd_dict
                        )
                        
                        # mcs_seed_loc = int(mcs_commands[0].split().index("--seed")) + 1
                        # mcs_nalign_loc = int(mcs_commands[0].split().index("--num-alignments")) + 1
                        # if options.multi_stage:
                        #     mcs_commands[0].replace(mcs_commands[0][mcs_nalign_loc], str(options.nalign_batch))
                            
                        #     base_seed = int(mcs_commands[0].split()[mcs_seed_loc])

                        #     nbatch = options.nalign // options.nalign_batch
                        #     res_nbatch = options.nalign - nbatch * options.nalign_batch

                        #     for i in range(nbatch):
                        #         cmd_copy = mcs_commands[0][:]
                        #         cmd_copy.replace(cmd_copy.split()[mcs_seed_loc], str(base_seed + i + 1))
                        #         mcs_commands.append(cmd_copy)

                        #     if res_nbatch != 0:
                        #         cmd_copy = mcs_commands[0][:]
                        #         cmd_copy.replace(cmd_copy.split()[mcs_nalign_loc], str(res_nbatch))
                        #         cmd_copy.replace(cmd_copy.split()[mcs_seed_loc], str(base_seed + nbatch + 1))
                        #         mcs_commands.append(cmd_copy)

                        if restart_step3:
                            if gene_tree is None and not restart_step2:
                                production_logger.info("### Restart RECUR from Step3 ###\n", extra={'to_file': True, 'to_console': True})
                            elif gene_tree is not None and not restart_step1:
                                production_logger.info("### Restart RECUR from Step2 ###\n", extra={'to_file': True, 'to_console': True})

                        prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                        production_logger.info(prepend + "Starting Monte-Carlo Simulation.", extra={'to_file': True, 'to_console': True})
                        production_logger.info("Using %d RECUR thread(s), %d IQ-TREE thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads),
                                            extra={'to_file': True, 'to_console': False})
                        production_logger.info("step2 iqtrees command: ", extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{mcs_commands[0]}\n", extra={'to_file': True, 'to_console': False})
                        
                        run_commands.RunCommand(
                            mcs_commands,
                            mcs_faDir,
                            env=my_env,
                            nthreads=options.recur_nthreads,
                            delete_files=True,
                            files_to_keep=["fasta", "fa"],
                            fd_limit=options.fd_limit,
                        )
                        prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                        production_logger.info(prepend + "Monte-Carlo Simulation complete.\n", extra={'to_file': True, 'to_console': True})
                    else:
                        if gene_tree is None:
                            production_logger.info("NOTE: With the existing Monte-Carlo simulated *.fa files, RECUR will skip Step3.\n",
                                                extra={'to_file': True, 'to_console': True})
                        else:
                            production_logger.info("NOTE: With the existing Monte-Carlo simulated *.fa files, RECUR will skip Step2.\n",
                                                extra={'to_file': True, 'to_console': True})

                    for file in util.iter_dir(real_phyDir):
                        file_path = os.path.join(real_phyDir, file)
                        if os.path.exists(file_path):
                            if file.endswith(".treefile.txt") or file.endswith(".treefile.log"):
                                os.remove(file_path)

                    afasta = random.choice(os.listdir(mcs_faDir))
                    fasta_dict, _, _ = filereader.ReadAlignment(os.path.join(mcs_faDir, afasta))
                    isnuc_fasta = util.CheckSequenceType([*fasta_dict.values()])

                    if options.usr_mcs_alnDir:
                        mcs_alnDir = options.usr_mcs_alnDir
                    else:
                        mcs_alnDir = mcs_faDir

                    if not gene_tree:
                        step3_info = f"Step4: Analysing recurrent substitutions"
                    else:
                        step3_info = f"Step3: Analysing recurrent substitutions"
                    production_logger.info(step3_info, extra={'to_file': True, 'to_console': True})
                    production_logger.info("="*len(step3_info), extra={'to_file': True, 'to_console': True})
                    production_logger.info(f"Results Directory: {recurrenceDir}\n", extra={'to_file': True, 'to_console': True})

                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Starting substitution matrices calculation for simulated phylogeny.", extra={'to_file': True, 'to_console': True})
                    production_logger.info("Using %d thread(s) for RECUR analysis" % options.nthreads, extra={'to_file': True, 'to_console': True})

                    mcs_results = parallel_task_manager.process_mcs_files_in_chunks(
                        mcs_alnDir,
                        parent_list,
                        child_list,
                        residue_dict,
                        options.nthreads,
                        isnuc_fasta,
                        options.sequence_type,
                        res_loc_list,
                        production_logger,
                        width,
                        dash_exist=dash_exist,
                        binary_sequence_dict=binary_combined_seq_dict,
                        update_cycle=options.update_cycle,
                        mcs_batch_size=options.mcs_batch_size
                    )

                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Substitution matrices complete.")
                    
                    if options.disk_save:
                        for file in util.iter_dir(mcs_alnDir):
                            file_path = os.path.join(mcs_alnDir, file)
                            if os.path.exists(file_path):
                                if file.endswith(".fasta") or file.endswith(".fa"):
                                    os.remove(file_path)
                    
                        filewriter.WriteMCSRecurrenceCountToFile(
                            mcs_results,
                            filehandler.mcs_dir
                        )
                        print("NOTE: You are running one the disk saving mode!")
                        
                        if options.multi_stage:
                            print(
                                "If you don't wish to use `multi-stage-recur.sh`, "
                                "you need to run `recur -f project_dir -cr -pam pval_adjust_method` "
                                "to obtain the recurrence list.\n"
                            )
                            sys.exit(0)

                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Starting to compute p values.")

                    R = mcs_count_greater(
                        mcs_results,
                        recurrence_list,
                        residue_dict,
                    )

                    recurrence_list_updated = update_recurrence_list(
                        R, 
                        options.nalign,
                        rec_loc_count_dict,
                        recurrence_list,
                        combined_prot_seqs_dict,
                        species_of_interest,
                        residue_dict_flip,
                        protein_len,
                        alpha=options.significance_level,
                        q=options.fdr_level,
                        method=options.pval_adjust_method,
                        pval_stats=options.pval_stats
                    )

                    filewriter.WriteRecurrenceList(
                        recurrence_list_updated, 
                        filehandler.GetRecurrenceListFN(recurrenceDir),
                        options
                    )
                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "p values computing complete.", extra={'to_file': True, 'to_console': True})

                    d_results = os.path.normpath(filehandler.GetResultsDirectory()) + os.path.sep
                    rec_results = os.path.normpath(filehandler.GetRecurrenceListFN(recurrenceDir))
                    production_logger.info("\nResults:\n    %s\n" % rec_results, extra={'to_file': True, 'to_console': True})

                    del parent_list, child_list, combined_prot_seqs_dict, alignment_dict, rec_loc_count_dict, recurrence_list
                    del recurrence_list_updated
                    gc.collect()

                    # util.log_memory_usage(f"after processing gene {gene}", production_logger)

                    end = time.perf_counter()
                    duration = end - start
                    if production_logger:
                        production_logger.info(f"Finished analysis of {os.path.basename(aln_path)} in {duration:.2f} seconds", extra={'to_file': True, 'to_console': True})
                    else:
                        print(f"Finished analysis of {os.path.basename(aln_path)} in {duration:.2f} seconds")

                except Exception as e:
                    print(f"\nERROR occurred during analysis of {os.path.basename(aln_path)}: {e}")
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

        # util.log_memory_usage("after final cleanup")

        end_main = time.perf_counter()
        duration_main = end_main - start_main

        if production_logger:
            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
            production_logger.info(prepend + "RECUR run completed\n", extra={'to_file': True, 'to_console': True})
            production_logger.info(f"*** RECUR finished in {duration_main:.2f} seconds ***", extra={'to_file': True, 'to_console': True})
            # Flush and close the handlers to ensure all logs are written
            for handler in production_logger.handlers:
                handler.flush()
                handler.close()
        else:
            print("RECUR run completed\n")
            print(f"*** RECUR finished in {duration_main:.2f} seconds ***")

        cleanup()
        kill_child_processes(os.getpid())

if __name__ == "__main__":
    main()
