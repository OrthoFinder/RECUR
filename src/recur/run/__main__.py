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
import gc
import time
import datetime
from rich import progress
import numpy as np
from collections import Counter
import copy
import warnings
import psutil
import signal
import traceback
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
import dendropy
from functools import partial
from recur.run import run_commands
from recur.utils import files, util, process_args
import logging
import warnings
# warnings.filterwarnings("ignore", module='dendropy')

def setup_environment() -> Tuple[Dict[str, str], str, str, Optional[str]]:
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

def CanRunCommand(command: str, env: Optional[Dict[str, str]] = None) -> bool:
    try:
        process = subprocess.Popen(command, env=env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
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


def initialise_recur(iqtree_version: Optional[str] = None) -> Dict[str, str]:
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
        if CanRunCommand(f"{local_iqtree2_path} --version", env=my_env):
            return my_env

    if iqtree_version == "conda" and conda_prefix:
        iqtree_path = shutil.which("iqtree2")
        if iqtree_path:
            print("\nConda version of IQ-TREE2 found.")
            print(f"IQ-TREE2 path: {iqtree_path}")
            return my_env

    if iqtree_version == "system" and not conda_prefix:
        iqtree_path = shutil.which("iqtree2")
        if iqtree_path:
            print("\nSystem-wide version of IQ-TREE2 found.")
            print(f"IQ-TREE2 path: {iqtree_path}")
            return my_env

    if conda_prefix and CanRunCommand("iqtree2 --version", env=my_env):
        print("Local IQ-TREE2 binary failed to run, falling back to the conda version.")
        return my_env

    if CanRunCommand("iqtree2 --version", env=my_env):
        print("Local IQ-TREE2 binary failed to run, falling back to system-wide binary.")
        return my_env

    print("Cannot proceed. IQ-TREE2 does not exist in either local bin or system-wide PATH.")
    print("Please ensure IQ-TREE2 is properly installed before running RECUR!\n")
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

def ParentChildRelation(treefile: str, 
                        outgroup_species: List[str],
                        n_species: int,
                        ) -> Tuple[Optional[str], List[str], List[str], List[str], str]:
    try:

        with open(treefile, 'r') as f:
            t = dendropy.Tree.get(file=f, schema="newick")

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
 
            return None, outgroup_species, [], [], "No root of interest found; tree structure might be incorrect."

        root_node = root_of_interest.label
        root_node = root_node if "/" not in root_node else root_node.split("/")[0]

        outgroup_subtree_species = [leaf.taxon.label for leaf in outgroup_mrca.leaf_iter()]

        taxon_count = 0
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

                taxon_count += 1

        if len(outgroup_subtree_species) != len(outgroup_species):
            outgroup_species = outgroup_subtree_species

        expected_relationships = 2 * (n_species - len(outgroup_species)) - 2

        if taxon_count != expected_relationships:
            return root_node, outgroup_species, parent_list, child_list, f"Taxon count {taxon_count} does not match expected relationships {expected_relationships}."
        
        return root_node, outgroup_species, parent_list, child_list, ""
    except BrokenPipeError:
        print("Broken pipe error.")
        return root_node, outgroup_species, parent_list, child_list, "Broken pipe error"
    except Exception as e:
        error_msg = f"ERROR in ParentChildRelation: {e}"
        print(error_msg)
        print(traceback.format_exc())
        return root_node, outgroup_species, parent_list, child_list, error_msg 

def count_mutations(parent_list: List[str],
                    child_list: List[str],
                    sequence_dict: Dict[str, str],
                    residue_dict: Dict[str, int]
                    ) -> Counter[Tuple[int, int, int]]:

    parent_num_list = [
        [residue_dict[res] for res in sequence_dict[parent]] 
        for parent in parent_list
    ]
    child_num_list = [
        [residue_dict[res] for res in sequence_dict[child]] 
        for child in child_list
    ]

    parent_array = np.array(parent_num_list)
    child_array = np.array(child_num_list)

    parent_child_diff = parent_array != child_array
    row_indices, col_indices = np.where(parent_child_diff)

    parent_res_id = parent_array[row_indices, col_indices]
    child_res_id = child_array[row_indices, col_indices]
    
    parent_child_tuples = [*zip(col_indices, parent_res_id, child_res_id)]
    rec_loc_count_dict = Counter(parent_child_tuples)

    return rec_loc_count_dict

def get_recurrence_list(rec_loc_count_dict: Counter[Tuple[int, int, int]], 
                        residue_dict_flip: Dict[int, str],
                        ) -> List[List[Union[int, str, float]]]:
    
    rec_loc_count_dict2 = copy.deepcopy(dict(rec_loc_count_dict))
    flipflop_dict = {}
    for key, _ in rec_loc_count_dict.most_common():
        res_loc, parent_id, child_id = key
        flipflop_key = (res_loc, child_id, parent_id)
        flipflop = rec_loc_count_dict.get(flipflop_key)

        if flipflop is not None:
            rec_loc_count_dict2.pop(flipflop_key, None)
            flipflop_dict[key] = flipflop

    recurrence_list: List[List[Union[int, str, float]]] = [
        [
            int(res_loc), 
            str(residue_dict_flip[parent_id]), 
            str(residue_dict_flip[child_id]), 
            int(recurrence), 
            int(flipflop_dict.get(key, 0))
         ] 
        for (res_loc, parent_id, child_id), recurrence in rec_loc_count_dict2.items() if int(recurrence) > 1
    ]
    return recurrence_list

def WorkerProcessAndCount(file: str, 
                          mcs_alnDir: str,
                          parent_list: List[str],
                          child_list: List[str],
                          isnuc_fasta: bool,
                          sequence_type: str, 
                          residue_dict: Dict[str, int], 
                          res_loc_list: List[int],
                          production_logger: logging.Logger,
                          ) -> Tuple[Dict[Tuple[int, int, int], int], str]:
    
    rec_loc_count_dict: Dict[Tuple[int, int, int], int] = {}
    error_msg = ""

    try:
        file_path = os.path.join(mcs_alnDir, file)
        mcs_combined_prot_seqs_dict, _ = files.FileReader.ReadAlignment(file_path)

        if isnuc_fasta:
            mcs_combined_prot_seqs_dict, _ = util.GetSeqsDict(mcs_combined_prot_seqs_dict, sequence_type)

        parent_num_list = [
            [residue_dict[res] for res in mcs_combined_prot_seqs_dict[parent]] 
            for parent in parent_list
        ]
        child_num_list = [
            [residue_dict[res] for res in mcs_combined_prot_seqs_dict[child]] 
            for child in child_list
        ]
        parent_array = np.array(parent_num_list)
        child_array = np.array(child_num_list)

        parent_child_diff = parent_array != child_array
        row_indices, col_indices = np.where(parent_child_diff)

        mask = np.isin(col_indices, res_loc_list)
        col_idx = col_indices[mask]
        row_idx = row_indices[mask]

        parent_res_id = parent_array[row_idx, col_idx]
        child_res_id = child_array[row_idx, col_idx]

        parent_child_tuples = [*zip(col_idx, parent_res_id, child_res_id)]
        rec_loc_count_dict = Counter(parent_child_tuples)

    except BrokenPipeError:
        error_msg = "Broken pipe error while processing file."
        print(error_msg)

    except Exception as e:
        error_msg = f"ERROR in WorkerProcessAndCount for mcs task {file}: {e}"
        print(error_msg)
        print(traceback.format_exc())
        production_logger.error(error_msg)
    
    return rec_loc_count_dict, error_msg


def process_mcs_files_in_chunks(mcs_alnDir: str, 
                                parent_list: List[str],
                                child_list: List[str], 
                                residue_dict: Dict[str, int], 
                                nthreads: int, 
                                isnuc_fasta: bool, 
                                sequence_type: str,
                                res_loc_list: List[int],
                                production_logger: logging.Logger,
                                window_width: int,
                                update_cycle: Optional[int] = None, 
                                mcs_batch_size: Optional[int] = None
                                ) -> List[Dict[Tuple[int, int, int], int]]:

    mcs_files = os.listdir(mcs_alnDir)
    total_file_count = len(mcs_files)

    results = []
    worker = partial(WorkerProcessAndCount, 
                     mcs_alnDir=mcs_alnDir,
                     parent_list=parent_list,
                     child_list=child_list,
                     isnuc_fasta=isnuc_fasta,
                     sequence_type=sequence_type,
                     residue_dict=residue_dict,
                     res_loc_list=res_loc_list,
                     production_logger=production_logger,
                     )
    if mcs_batch_size is not None:    
        batches = [mcs_files[i:i + mcs_batch_size] for i in range(0, total_file_count, mcs_batch_size)]
    
    if update_cycle is not None:
        mcs_progress = progress.Progress(
        progress.TextColumn("[progress.description]{task.description}"),
        progress.BarColumn(bar_width=window_width // 2),
        progress.SpinnerColumn(),
        progress.MofNCompleteColumn(),
        progress.TimeElapsedColumn(),
        transient=False,
        # progress.TextColumn("{task.completed}/{task.total}")
        )
        task = mcs_progress.add_task("[magenta]Processing...", total=total_file_count)

        mcs_progress.start()
    try:
        with ProcessPoolExecutor(max_workers=nthreads) as executor:
            
            if mcs_batch_size is not None:   
                futures = [executor.submit(worker, file_data) for batch in batches for file_data in batch]
            else:
                futures = [executor.submit(worker, file_data) for file_data in mcs_files]
            
                
            for i, future in enumerate(as_completed(futures)):
                try:
                    result, error_msg = future.result()
                    if error_msg:
                        print(error_msg)
                        production_logger.error(error_msg)
                        break

                    if result:
                        results.append(result)
                            
                except Exception as e:
                    error_msg = f"ERROR during processing: {e}"
                    print(error_msg)
                    print(traceback.format_exc())
                    production_logger.error(error_msg)
                    break

                finally:
                    if update_cycle is not None:
                        if (i + 1) % update_cycle == 0:
                            mcs_progress.update(task, advance=update_cycle)
            
            if update_cycle is not None:
                mcs_progress.stop()

    except Exception as e:
        error_msg = f"ERROR during processing files: {e}"
        print(error_msg)
        print(traceback.format_exc())
        production_logger.error(error_msg)

    return results

def compute_p_values(mcs_results: List[Dict[Tuple[int, int, int], int]],
                     recurrence_list: List[List[Union[str, int, float]]],
                     residue_dict: Dict[str, int],
                     nalign: int) -> List[List[Union[str, int, float]]]:
    try:
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
        
            p_value = (count_greater + 1) / (nalign + 1)
            rec_list.append(np.round(p_value, len(str(nalign))))
            # print(f"{rec_list} - {count_greater}")

    except Exception as e:
        error_msg = f"ERROR during compute_p_values: {e}"
        print(error_msg)
        print(traceback.format_exc())
    finally:
        return recurrence_list

def update_recurrence_list(res_loc_count_dict: Dict[Tuple[int, int, int], int],
                            recurrence_list: List[List[Union[str, int, float]]],
                            combined_prot_seqs_dict: Dict[str, str],
                            species_of_interest: List[str],
                            residue_dict_flip: Dict[int, str],
                            protein_len: int) -> List[List[Union[str, int, float]]]:
    
    extant_seq = {species: seq for species, seq in combined_prot_seqs_dict.items() if species in species_of_interest}
    ident_dict = {}

    res_loc_info_dict= util.get_sorted_res_loc_info(res_loc_count_dict, protein_len)
    for rec_loc, res in enumerate(zip(*extant_seq.values())):
        ident_dict[rec_loc] = res

    for rec_list in recurrence_list:
        res_loc = int(rec_list[0])
        parent_child = []
        counts = []
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
        print(f"RECUR:v{__version__}")
        sys.exit()
 
    start_main = time.perf_counter()

    try:
        input_command = "RECUR command: " + " ".join(sys.argv)  + "\n"
        
        options, alnDir, alnPath, resultsDir_nonDefault = process_args.ProcessArgs(args)
        iqtree_version = options.iqtree_version if options.iqtree_version else "system"
        my_env = initialise_recur(iqtree_version)

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
                override = True

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
                        production_logger.info(prepend + "Ran IQ-TREE to build the gene trees for the real phylogeny", extra={'to_file': True, 'to_console': True})
                        production_logger.info("Using %d RECUR thread(s), %d IQ-TREE2 thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads), extra={'to_file': True, 'to_console': True})
                    
                    if gene_tree is None:
                        asr = False
                        fix_branch_length = False
                    else:
                        fix_branch_length = options.fix_branch_length
                    
                    commands = run_commands.GetGeneTreeBuildCommands([aln_path], 
                                                        real_phyDir, 
                                                        options.evolution_model,
                                                        options.iqtree_nthreads,
                                                        phy_seed=options.seed,
                                                        asr=asr,
                                                        sequence_type=options.sequence_type,
                                                        gene_tree=gene_tree,
                                                        bootstrap=options.bootstrap,
                                                        sh_alrt=options.bootstrap,
                                                        fix_branch_length=fix_branch_length,
                                                        )
                    
                    if gene_tree is None:
                        production_logger.info(f"step1 iqtree2 gene tree building command: ", extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    else:
                        production_logger.info(f"step1 iqtree2 ancestral state reconstruction command: ", extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    run_commands.RunCommand(commands, 
                                            real_phyDir,
                                            env=my_env, 
                                            nthreads=options.recur_nthreads,
                                            delete_files=True,
                                            files_to_keep=["state", "treefile", "iqtree"],
                                            fd_limit=options.fd_limit)
                    
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
                        if not restart_step2 and not exist_state:             
                            production_logger.info("NOTE: RECUR is forced to restart from Step2 due to the missing statefile\n", extra={'to_file': True, 'to_console': True})
                        
                        elif not restart_step1 and restart_step2:
                            production_logger.info("### Restart RECUR from Step2 ###\n", extra={'to_file': True, 'to_console': True})
                            
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
                                                            asr=True,
                                                            fix_branch_length=options.fix_branch_length)
                        
                        production_logger.info(f"step2 iqtree2 ancestral state reconstruction command: ",  extra={'to_file': True, 'to_console': False})
                        production_logger.info(f"{commands[0]}\n", extra={'to_file': True, 'to_console': False})
        
                        run_commands.RunCommand(commands, 
                                                real_phyDir,
                                                env=my_env,
                                                nthreads=options.recur_nthreads, 
                                                delete_files=True,
                                                files_to_keep=["state", "treefile", "iqtree"],
                                                fd_limit=options.fd_limit)

                        statefile = filehandler.GetStateFileFN()
                        if not restart_step2 and restart_step3:
                            production_logger.info(f"Ancestral state reconstruction complete\n", extra={'to_file': True, 'to_console': True})
                        else:
                            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                            production_logger.info(prepend + f"Ancestral state reconstruction complete\n", extra={'to_file': True, 'to_console': True})
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
              
                root_node, outgroup_species, parent_list, child_list, error_msg = ParentChildRelation(treefile, 
                                                                                                      outgroup_mrca, 
                                                                                                      n_species, 
                                                                                                      )
                
                if error_msg:
                    production_logger.error(error_msg)
                    continue

                if outgroup_mrca and len(outgroup_mrca) != len(outgroup_species):
                    warnings.warn(f"Outgroup species provided not monophyletic. Outgroups will be updated. Please find the updated outgroups in the log file.")
                    production_logger.info(f"Updated outgroups: {outgroup_mrca}", extra={'to_file': True, 'to_console': False})
                    outgroup_mrca = outgroup_species


                rec_loc_count_dict = count_mutations(parent_list, 
                                                    child_list,
                                                    combined_prot_seqs_dict,
                                                    residue_dict)
                    
                production_logger.info(f"Root of species of interest: {root_node}", extra={'to_file': True, 'to_console': True})
                production_logger.info(f"Substitution matrix output: {filehandler.GetMutMatrixDir()}\n", extra={'to_file': True, 'to_console': True})

                filewriter.WriteMutMatrix(rec_loc_count_dict,
                                          residue_dict_flip, 
                                          protein_len,
                                          filehandler.GetMutCountMatricesFN(),
                                          filehandler.GetAccumMutCountMatricesFN())
                
                recurrence_list = get_recurrence_list(rec_loc_count_dict, residue_dict_flip)

                if len(recurrence_list) == 0:
                    production_logger.info(f"ATTENTION: No recurrence has identified for gene {gene}! Monte-Carlo Simiatlion will be SKIPPED!\n", 
                                           extra={'to_file': True, 'to_console': True})
                    continue 

                res_loc_list = [int(res_list[0]) for res_list in recurrence_list]

                if gene_tree is None:
                    step2_info = f"Step3: Simulating Sequence Evolution with {options.nalign} replicates"
                else: 
                    step2_info = f"Step2: Simulating Sequence Evolution with {options.nalign} replicates"
                    
                production_logger.info(step2_info, extra={'to_file': True, 'to_console': True})
                production_logger.info("="*len(step2_info), extra={'to_file': True, 'to_console': True})
                step2_results_info = f"Results Directory: {filehandler.GetMCSimulationDir()}\n"
                production_logger.info(step2_results_info, extra={'to_file': True, 'to_console': True})

                mcs_faDir = filehandler.GetMCSimulationDir()
                if len(os.listdir(mcs_faDir)) != options.nalign or \
                    (restart_step1 or restart_step2 or restart_step3) or override:

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
                    if restart_step3:
                        if gene_tree is None and not restart_step2:
                            production_logger.info("### Restart RECUR from Step3 ###\n", extra={'to_file': True, 'to_console': True})
                        elif gene_tree is not None and not restart_step1:
                            production_logger.info("### Restart RECUR from Step2 ###\n", extra={'to_file': True, 'to_console': True})

                    prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                    production_logger.info(prepend + "Starting Monte-Carlo Simulation.", extra={'to_file': True, 'to_console': True})
                    production_logger.info("Using %d RECUR thread(s), %d IQ-TREE2 thread(s)" % ( options.recur_nthreads, options.iqtree_nthreads), 
                                        extra={'to_file': True, 'to_console': False})
                    production_logger.info("step2 iqtrees command: ", extra={'to_file': True, 'to_console': False})
                    production_logger.info(f"{mcs_commands[0]}\n", extra={'to_file': True, 'to_console': False})

                    run_commands.RunCommand(mcs_commands, 
                                            mcs_faDir,
                                            env=my_env, 
                                            nthreads=options.recur_nthreads,
                                            delete_files=True,
                                            files_to_keep=["fasta", "fa"],
                                            fd_limit=options.fd_limit)
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
                fasta_dict, _ = filereader.ReadAlignment(os.path.join(mcs_faDir, afasta))
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
                
                mcs_results = process_mcs_files_in_chunks(mcs_alnDir, 
                                                parent_list,
                                                child_list, 
                                                residue_dict, 
                                                options.nthreads,
                                                isnuc_fasta,
                                                options.sequence_type,
                                                res_loc_list, 
                                                production_logger,
                                                width, 
                                                update_cycle=options.update_cycle,
                                                mcs_batch_size=options.mcs_batch_size)

                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Substitution matrices creation complete.\n")
                
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "Starting compute p values.")


                recurrence_list_pvalue = compute_p_values(mcs_results, 
                                                        recurrence_list,
                                                        residue_dict,
                                                        options.nalign)
                
                recurrence_list_updated = update_recurrence_list(rec_loc_count_dict,
                                                                recurrence_list_pvalue,
                                                                combined_prot_seqs_dict,
                                                                species_of_interest,
                                                                residue_dict_flip,
                                                                protein_len)
                
                filewriter.WriteRecurrenceList(recurrence_list_updated, filehandler.GetRecurrenceListFN(recurrenceDir))
                prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + ": "
                production_logger.info(prepend + "p values computing complete.", extra={'to_file': True, 'to_console': True})

                d_results = os.path.normpath(filehandler.GetResultsDirectory()) + os.path.sep
                rec_results = os.path.normpath(filehandler.GetRecurrenceListFN(recurrenceDir))
                production_logger.info("\nResults:\n    %s\n" % rec_results, extra={'to_file': True, 'to_console': True})
                
                del parent_list, child_list, combined_prot_seqs_dict, alignment_dict, rec_loc_count_dict, recurrence_list
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
