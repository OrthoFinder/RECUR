# -*- coding: utf-8 -*-

import concurrent.futures
import contextlib
import os
import subprocess
import sys
import threading
import warnings
import platform
from typing import Dict, List, Optional, Set
from rich import print


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
                        fd_limit: Optional[int] = None):

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
