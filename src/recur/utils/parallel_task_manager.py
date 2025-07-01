# -*- coding: utf-8 -*-

import concurrent.futures
import contextlib
import os
import subprocess
import sys
import traceback
import threading
import warnings
import logging
from collections import Counter
from functools import partial
import platform
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from rich import print, progress

from recur.utils import util, files

# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)


lock = threading.RLock()

def drop_privileges():
    if platform.system() == "Windows":
        return
    else:
        os.setgid(os.getgid())
        os.setuid(os.getuid())

def RunCommand(
        command: str,
        env: Optional[Dict[str, str]] = None,
        qPrintOnError: bool = False,
        qPrintStderr: bool = True
    ) -> int:
    
    kwargs = {}
    if platform.system() != "Windows":
        kwargs["preexec_fn"] = drop_privileges

    try:
        popen = subprocess.Popen(
            command,
            env=env,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            **kwargs 
        )
        
        stdout, stderr = popen.communicate()

        if qPrintOnError and popen.returncode != 0:
            with lock:
                print(f"\nERROR: external program returned an error code: {popen.returncode}")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout.decode()}")
                print(f"stderr:\n{stderr.decode()}")

        elif qPrintStderr and len(stderr) > 0:
            with lock:
                print("\nWARNING: program produced output to stderr")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout.decode()}")
                print(f"stderr:\n{stderr.decode()}")

        return popen.returncode

    except Exception as e:
        with lock:
            print(f"Exception occurred while running command: {command}")
            print(f"Exception: {e}")
        return -1

def clean_up_files(fileDir: str,
                   processed_files: Set[str],
                   already_deleted_files: Set[str],
                   files_to_keep: Optional[List[str]] = None,
                   files_to_remove: Optional[List[str]] = None) -> None:
    with lock:
        for file in os.listdir(fileDir):
            file_path = os.path.join(fileDir, file)
            if os.path.isfile(file_path) and os.path.exists(file_path):
                if file_path in already_deleted_files:
                    continue
                if files_to_keep:
                    if file_path in processed_files or \
                       file.rsplit(".", 1)[-1].lower() not in files_to_keep:
                        with contextlib.suppress(FileNotFoundError):
                            os.remove(file_path)
                            already_deleted_files.add(file_path)

                elif files_to_remove:
                    if file_path in processed_files or file.rsplit(".", 1)[-1].lower() in files_to_remove:
                        with contextlib.suppress(FileNotFoundError):
                            os.remove(file_path)
                            already_deleted_files.add(file_path)

def set_file_descriptor_limit(new_limit: int) -> None:
    try:
        if os.name != 'nt':
            import resource
            soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
            print(f"Current file descriptor limits: soft={soft_limit}, hard={hard_limit}")

            if new_limit > soft_limit:
                resource.setrlimit(resource.RLIMIT_NOFILE, (new_limit, hard_limit))
                print(f"New file descriptor limits: soft={new_limit}, hard={hard_limit}")
            else:
                print("New limit is not higher than the current soft limit.")
        else:
            print("File descriptor limit functions are not available on Windows.")
    except AttributeError:
        print("File descriptor limit functions not available on this platform.")
    except ValueError as e:
        print(f"Invalid limit value: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

def RunParallelCommands(nProcesses: int,
                        commands: List[str],
                        fileDir: str,
                        env: Optional[Dict[str, str]] = None,
                        delete_files: bool = False,
                        files_to_keep: Optional[List[str]] = None,
                        files_to_remove: Optional[List[str]] = None,
                        q_print_on_error: bool = False,
                        q_always_print_stderr: bool = False,
                        fd_limit: Optional[int] = None,
                        ):

    processed_files: Set[str] = set()
    already_deleted_files: Set[str] = set()

    try:
        if fd_limit is not None:
            if sys.platform.startswith('linux') or sys.platform == 'darwin':
                set_file_descriptor_limit(fd_limit)
            else:
                warnings.warn(f"File descriptor limit adjustment is not supported on {sys.platform}.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=nProcesses) as executor:
            futures = [executor.submit(RunCommand,
                                       cmd,
                                       env,
                                       q_print_on_error,
                                       q_always_print_stderr)
                       for cmd in commands]
            # concurrent.futures.wait(futures)
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    if result != 0:
                        if q_print_on_error:
                            print(f"ERROR occurred with command: {future}")
                except Exception as exc:
                    if q_print_on_error:
                        print(f"Exception occurred with command: {future}, Error: {exc}")

        if delete_files:
            # util.PrintTime("Cleaning up the directory.")
            clean_up_files(fileDir, processed_files, already_deleted_files, files_to_keep, files_to_remove)

    except KeyboardInterrupt:
        for future in futures:
            future.cancel()

def WorkerProcessAndCount(file: str,
                          mcs_alnDir: str,
                          parent_list: List[str],
                          child_list: List[str],
                          isnuc_fasta: bool,
                          sequence_type: str,
                          residue_dict: Dict[str, int],
                          res_loc_list: List[int],
                          production_logger: logging.Logger,
                          dash_exist: bool = False,
                          binary_sequence_dict: Dict[str, str] = {},
                          ) -> Tuple[Dict[Tuple[int, int, int], int], str]:

    rec_loc_count_dict: Dict[Tuple[int, int, int], int] = {}
    error_msg = ""

    try:
        file_path = os.path.join(mcs_alnDir, file)
        mcs_combined_prot_seqs_dict, _, _ = files.FileReader.ReadAlignment(file_path)

        if isnuc_fasta:
            mcs_combined_prot_seqs_dict, _ = util.GetSeqsDict(mcs_combined_prot_seqs_dict, sequence_type)

        parent_num_list = [
            [residue_dict.get(res, 0) for res in mcs_combined_prot_seqs_dict[parent]]
            for parent in parent_list
        ]
        child_num_list = [
            [residue_dict.get(res, 0) for res in mcs_combined_prot_seqs_dict[child]]
            for child in child_list
        ]

        parent_array = np.array(parent_num_list)
        child_array = np.array(child_num_list)

        del mcs_combined_prot_seqs_dict, parent_num_list, child_num_list

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

        mask = np.isin(col_indices, res_loc_list)
        col_idx = col_indices[mask]
        row_idx = row_indices[mask]

        parent_res_id = parent_array[row_idx, col_idx]
        child_res_id = child_array[row_idx, col_idx]

        # parent_mask = np.isin(parent_res_id, util.reserved_chars_index, invert=True)
        parent_mask = np.where(parent_res_id != 0)
        parent_res_id = parent_res_id[parent_mask]
        child_res_id = child_res_id[parent_mask]
        col_idx = col_idx[parent_mask]

        # child_mask = np.isin(child_res_id, util.reserved_chars_index, invert=True)
        child_mask = np.where(child_res_id != 0)
        parent_res_id = parent_res_id[child_mask]
        child_res_id = child_res_id[child_mask]
        col_idx = col_idx[child_mask]

        parent_child_tuples = [*zip(col_idx, parent_res_id, child_res_id)]
        rec_loc_count_dict = Counter(parent_child_tuples)
        del parent_child_diff, mask, row_idx, col_idx, row_indices, col_indices, \
            parent_res_id, child_res_id, parent_child_tuples, parent_array, child_array

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
                                dash_exist: bool = False,
                                binary_sequence_dict: Optional[Dict[str, str]] = None,
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
                     dash_exist=dash_exist,
                     binary_sequence_dict=binary_sequence_dict,
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
        with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:

            if mcs_batch_size is not None:
                futures = [executor.submit(worker, file_data) for batch in batches for file_data in batch]
            else:
                futures = [executor.submit(worker, file_data) for file_data in mcs_files]


            for i, future in enumerate(concurrent.futures.as_completed(futures)):
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