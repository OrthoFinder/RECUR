# -*- coding: utf-8 -*-
#
import os
import sys
import datetime
from recur.citation import citation, print_citation
from recur.utils import files
from typing import List, Dict, Tuple, Optional, Any, Callable
from recur import genetic_codes
from importlib import resources as impresources
import psutil
import logging
import logging.config
import yaml
import tracemalloc
import traceback

residues = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

class ConsoleOnlyFilter(logging.Filter):
    def filter(self, record):
        return getattr(record, 'to_console', True)

class FileOnlyFilter(logging.Filter):
    def filter(self, record):
        return getattr(record, 'to_file', True)
    
def setup_logging(log_folder: str, 
                  mode: Optional[str] = "w", 
                  formatter: Optional[str] = "brief") -> logging.Logger:
    try:
        base_path = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir))
        log_config_path = os.path.join(base_path, 'logging_config.yaml')
        
        if not os.path.exists(log_config_path):
            print(f"Contents of the base path ({base_path}): {os.listdir(base_path)}")  # Debugging line
            raise FileNotFoundError(f"Log configuration file not found at {log_config_path}")
    except Exception as e:
        print(f"Error resolving log config path manually: {str(e)}")  # Debugging line
        log_config_path = None

    production_logger = logging.getLogger("production")

    if production_logger.hasHandlers():
        production_logger.handlers.clear()

    if log_config_path and os.path.exists(log_config_path):
        try:
            with open(log_config_path, 'rt') as f:
                config = yaml.safe_load(f.read())

            log_path = os.path.join(log_folder, 'Log.txt')
            if not os.path.exists(log_folder):
                os.makedirs(log_folder)

            config['handlers']['file']['filename'] = log_path
            config['handlers']['file']['formatter'] = formatter
            config['handlers']['file']['mode'] = mode

            # Dynamically replace the placeholder with the correct module path
            config_str = yaml.dump(config)
            config_str = config_str.replace('placeholder.ConsoleOnlyFilter', f'{__name__}.ConsoleOnlyFilter')
            config_str = config_str.replace('placeholder.FileOnlyFilter', f'{__name__}.FileOnlyFilter')
            config = yaml.safe_load(config_str)

            logging.config.dictConfig(config)
            production_logger = logging.getLogger("production")

        except Exception as e:
            print(f"Error loading logging configuration: {str(e)}")
            logging.basicConfig(level=logging.DEBUG)  # Fallback to basic configuration
            production_logger = logging.getLogger("production")
            production_logger.error("Logging configuration failed. Fallback to basic config.", exc_info=True)
    else:
        print(f"Log configuration file does not exist at: {log_config_path}")  # Debugging line
        logging.basicConfig(level=logging.DEBUG)
        production_logger = logging.getLogger("production")
        production_logger.error("Logging configuration file not found. Using basic config.")

    return production_logger

def log_memory_usage(tag="", production_logger=None):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 3)
    if production_logger:
        production_logger.info(f"Memory usage {tag}: {mem:.2f} GB\n")
    else:
        print(f"Memory usage {tag}: {mem:.2f} GB\n")

def delete_files_in_directory(dirpath: str) -> None:
   try:
     with os.scandir(dirpath) as entries:
       for entry in entries:
         if entry.is_file():
            os.unlink(entry.path)
   except OSError:
     print("ERROR occurred while deleting files.")
     raise

def measure_memory_usage(func: Callable[..., Any], *args: Any, **kwargs: Any) -> int:
    tracemalloc.start()
    func(*args, **kwargs)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak

def adjust_batch_size(current_batch_size: int, max_memory_usage: int, total_memory: int) -> int:
    target_memory_usage = total_memory * 0.75  # Target to use up to 75% of total memory
    if max_memory_usage > target_memory_usage:
        return max(1, current_batch_size // 2)
    else:
        return min(current_batch_size + 1, 1000)  # Ensure batch size doesn't grow indefinitely

def determine_optimal_batch_size(data: Any, 
                                 task_type: str, 
                                 total_memory: int) -> int:

    if task_type == 'file_processing':
        memory_per_task = measure_memory_usage(files.FileReader.ReadAlignment, data)
    elif task_type == 'mutation_count':
        memory_per_task = measure_memory_usage(len, data) * 500  # Measure memory usage for mutation counting
    else:
        raise ValueError("Unknown task type")

    memory_usage = memory_per_task if memory_per_task != 0 else 1

    optimal_batch_size = max(1, (total_memory // 4) // memory_usage)
    
    return int(min(optimal_batch_size, 100))

def print_centered_text(width, text):
    text_length = len(text)
    if text_length >= width:
        print(text)
    else:
        dashes_each_side = (width - text_length - 2) // 2
        dashes = '-' * dashes_each_side
        print()
        print(f"{dashes} {text} {dashes}")
        print()

def get_system_info():
    num_cpus = psutil.cpu_count(logical=False) 
    num_logical_cpus = psutil.cpu_count(logical=True)
    print("Machine Information:")
    print(f"Number of Physical CPUs: {num_cpus}, Number of Logical CPUs: {num_logical_cpus}")
    print()
    print("RECUR Start with Running Processes:")
    recur_proc = 0
    for proc in psutil.process_iter(['pid', 'name', 'num_threads', 'cmdline', 'ppid']):
        try:
            if "recur" in proc.info['name']:
                recur_proc += 1
                print(f"Process ID: {proc.info['pid']}, Name: {proc.info['name']}, "
                      f"Threads: {proc.info['num_threads']}, "
                      f"Parent Process ID: {proc.info['ppid']}")
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    
    print(f"Total RECUR processes found: {recur_proc}")
    if recur_proc > 2:
        print("More than two processes found, please cleanup the RECUR processes and rerun RECUR!")
        print("If you are using WSL2, you can run `wsl --shutdown` to reboot the subsystem.")

def residue_table() -> Tuple[Dict[str, int], Dict[int, str]]:
    residue_pos = [*range(len(residues))]
    residue_tuples = [*zip(residues, residue_pos)]
    residue_dict = {k: int(v) for k, v in residue_tuples}
    residue_dict_flip = {int(v): k for k, v in residue_tuples}
    return residue_dict, residue_dict_flip

def CheckSequenceType(alignments: List[str]) -> bool:
    nuc = {"A", "C", "T", "G", "-"}
    diff_aa_nuc = set(residues) - nuc
    isnuc = True
    for aln in alignments:
        unique_res = set(aln)
        diff_aln = unique_res - nuc 
        if any(item in diff_aa_nuc for item in diff_aln):
            isnuc = False
            break
    return isnuc

def ImportCodon(code_name: str) -> Dict[str, str]:
    
    code_fn = code_name.lower() + ".txt"
    codon_table = {}
    try:
        with impresources.open_text(genetic_codes, code_fn) as reader:
            for line in reader:
                codon, letter = line.strip().split()
                codon_table[codon] = letter
            
    except FileNotFoundError as e:
            print(f"File not found: {e.filename}")
    
    return codon_table

def Translate(seq: str, table: Dict[str, str]) -> str:

    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in table:
                protein += '-'
            else:
                protein += table[codon]

    return protein

def GetSeqsDict(dna_seq_dict: Dict[str, str], sequence_type: str) -> Tuple[Dict[str, str], int]:

    codon_table = ImportCodon(sequence_type)     #NCBI's genetic code 11 (for plant plastids)

    prot_sequence_dict = {}
    for node, seq in dna_seq_dict.items():
        prot_seq = Translate(seq, codon_table)
        prot_sequence_dict[node] = prot_seq
    protein_len = len(prot_seq)

    return prot_sequence_dict, protein_len

def PrintTime(message: str) -> None:
    print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message))
    sys.stdout.flush()

def Fail():
    sys.stderr.flush()
    print(traceback.format_exc())
    sys.exit(1)
               
def GetDirectoryName(baseDirName: str, 
                     i: int, 
                     sequence_type: str, 
                     extended_filename: bool) -> str:
    
    if not extended_filename:
        if i == 0:
            return baseDirName + os.sep
        else:
            return baseDirName + ("_%d" % i) + os.sep
    else:
        sequence_type = sequence_type if sequence_type else ""
        extension = sequence_type

        if i == 0:
            return baseDirName + "_" + extension + os.sep 
        else:
            return baseDirName + ("_%d" % i) + "_" + extension + os.sep

def CreateNewWorkingDirectory(baseDirectoryName: str, 
                              sequence_type: str, 
                              qDate: bool = True,
                              extended_filename: bool = False,
                              keepprev: bool = False) -> str:
    iAppend = 0
    newDirectoryName = GetDirectoryName(baseDirectoryName,
                                        iAppend,
                                        sequence_type,
                                        extended_filename)
    if keepprev:
        dateStr = datetime.date.today().strftime("%b%d") if qDate else ""
        baseDirectoryName = baseDirectoryName  + dateStr
        while os.path.exists(newDirectoryName):
            iAppend += 1
            newDirectoryName = GetDirectoryName(baseDirectoryName, 
                                                iAppend, 
                                                sequence_type,
                                                extended_filename)

    os.makedirs(newDirectoryName, exist_ok=True)

    return newDirectoryName
   
def WriteCitation(d: str) -> None:
    with open(d + "Citation.txt", 'w') as outfile:
        outfile.write(citation)

def PrintCitation(d: Optional[str] = None) -> None:
    if d is not None: 
        WriteCitation(d)
    print()
    print(print_citation)  
