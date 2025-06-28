# -*- coding: utf-8 -*-
#
import datetime
import logging
import logging.config
import os
import re
import sys
import traceback
from importlib import resources as impresources
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import psutil
import yaml
from rich import print

from recur import genetic_codes, helpinfo
from recur.citation import citation, print_citation

reserved_chars = ['-', 'B', 'O', 'J', 'Z', 'U', 'X'] # B O J X Z U ':', '*', '.', '+', ' '
# reserved_chars_index = [*range(20, 20 + len(reserved_chars))]

residues = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', \
            'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

special_chars = reserved_chars + [':', '*', '.', '+', ' ']
legal_chars = residues + reserved_chars

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
            print(f"Contents of the base path ({base_path}): {os.listdir(base_path)}")
            raise FileNotFoundError(f"Log configuration file not found at {log_config_path}")
    except Exception as e:
        print(f"Error resolving log config path manually: {str(e)}")
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

            config_str = yaml.dump(config)
            config_str = config_str.replace('placeholder.ConsoleOnlyFilter', f'{__name__}.ConsoleOnlyFilter')
            config_str = config_str.replace('placeholder.FileOnlyFilter', f'{__name__}.FileOnlyFilter')
            config = yaml.safe_load(config_str)

            logging.config.dictConfig(config)
            production_logger = logging.getLogger("production")

        except Exception as e:
            print(f"Error loading logging configuration: {str(e)}")
            logging.basicConfig(level=logging.DEBUG)
            production_logger = logging.getLogger("production")
            production_logger.error("Logging configuration failed. Fallback to basic config.", exc_info=True)
    else:
        print(f"Log configuration file does not exist at: {log_config_path}")
        logging.basicConfig(level=logging.DEBUG)
        production_logger = logging.getLogger("production")
        production_logger.error("Logging configuration file not found. Using basic config.")

    return production_logger

def log_memory_usage(tag: str = "", production_logger: Optional[logging.Logger] = None) -> None:
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

def print_centered_text(width: int, text: str) -> None:
    text_length = len(text)
    if text_length >= width:
        print(text)
    else:
        dashes_each_side = (width - text_length - 2) // 2
        dashes = '-' * dashes_each_side
        print()
        print(f"{dashes} {text} {dashes}")
        print()

def get_system_info() -> None:
    num_cpus = psutil.cpu_count(logical=False)
    num_logical_cpus = psutil.cpu_count(logical=True)
    print("\nMachine Information:")
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
    residue_dict = dict(zip(residues, range(1, len(residues)+1)))
    residue_dict_flip = dict(zip(range(1, len(residues)+1), residues))
    return residue_dict, residue_dict_flip

def CheckSequenceType(alignments: List[str]) -> bool:
    nuc = {"A", "C", "T", "G"}
    diff_aa_nuc = set(residues) - nuc
    isnuc = True
    for aln in alignments:
        unique_res = set(aln)
        unique_res.difference_update(set(special_chars))
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

def GetSeqsDict(dna_seq_dict: Dict[str, str],
                sequence_type: str) -> Tuple[Dict[str, str], int]:

    codon_table = ImportCodon(sequence_type)     #NCBI's genetic code 11 (for plant plastids)

    prot_sequence_dict = {
        node: Translate(seq, codon_table)
        for node, seq in dna_seq_dict.items()
    }
    protein_len = len(next(iter(prot_sequence_dict.values())))

    return prot_sequence_dict, protein_len

def PrintTime(message: str) -> None:
    print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message))
    sys.stdout.flush()

def Fail():
    sys.stderr.flush()
    print(traceback.format_exc())
    helpinfo.PrintHelp()
    sys.exit(1)

def GetFileName(baseFileName: str,
                i: int,
                sequence_type: str,
                extended_filename: bool) -> str:
    
    if not extended_filename:
        if i == 0:
            return baseFileName + ".tsv"
        else:
            return baseFileName + ("_%d" % i) + ".tsv"
    else:
        sequence_type = sequence_type if sequence_type else ""
        extension = sequence_type

        if i == 0:
            return baseFileName + "." + extension + ".tsv"
        else:
            return baseFileName + ("_%d" % i) + "." + extension + ".tsv"

def CreateNewFileName(baseFileName: str,
                    sequence_type: str,
                    qDate: bool = True,
                    extended_filename: bool = False,
                    keepprev: bool = False) -> str:
    iAppend = 0
    newFileName = GetFileName(baseFileName,
                                    iAppend,
                                    sequence_type,
                                    extended_filename)
    if keepprev:
        dateStr = datetime.date.today().strftime("%b%d") if qDate else ""
        baseFileName = baseFileName  + "." + dateStr
        while os.path.exists(newFileName):
            iAppend += 1.
            newFileName = GetFileName(baseFileName,
                                            iAppend,
                                            sequence_type,
                                            extended_filename)

    return newFileName


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
        dateStr = "." + datetime.date.today().strftime("%b%d") if qDate else ""
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

def iter_dir(d: Optional[str] = None) -> Iterator[str]:
    if d is None:
        Fail()
    with os.scandir(d) as entries:
        for entry in entries:
            if entry.is_file():
                yield entry.name

def get_sorted_res_loc_info(res_loc_count_dict: Dict[Tuple[int, int, int], int],
                            protein_len: int) -> Dict[int, List[Tuple[int, int, int]]]:


    res_loc_info_dict: Dict[int, List[Tuple[int, int, int]]] = {
        res_loc: [] for res_loc in range(protein_len)
    }

    for (res_loc, parent_id, child_id), recurrence in res_loc_count_dict.items():
        res_loc_info_dict[res_loc].append((parent_id, child_id, recurrence))

    res_loc_info_dict_sorted: Dict[int, List[Tuple[int, int, int]]] = {
        res_loc: sorted(val, key=lambda x: x[-1], reverse=True) if val else []
        for res_loc, val in res_loc_info_dict.items()
    }

    return res_loc_info_dict_sorted

def CheckOutgroups(outgroups_mrca: List[str], alignment_dict: Dict[str, str]) -> Tuple[List[str], bool]:

    outgroups = []
    preserve_underscores = False
    alignment_names = {}
    COMPILE = re.compile(r"\s+")
    for key in alignment_dict:
        if "_" in key:
            preserve_underscores = True
        new_key = COMPILE.sub("_", key.strip().lower())
        alignment_names[new_key] = key

    for outgroup in outgroups_mrca:
        new_outgroup = COMPILE.sub("_",  outgroup.lower().strip())
        if new_outgroup in alignment_names:
            outgroups.append(alignment_names[new_outgroup])
        else:
            print(f"Outgroup sequence {outgroup} is not found in the multiple sequence alignment.")
            print(traceback.format_exc())
            sys.exit(1)

    return outgroups, preserve_underscores

def ConvertToBinary(alignment_dict: Dict[str, str]) -> Dict[str, str]:
    msa_binary = {}
    for label, seq in alignment_dict.items():
        binary_seq = [
            1 if res in residues else 0 if res in reserved_chars else 2 for res in seq
        ]
        msa_binary[label] = ''.join(map(str, binary_seq))

    return msa_binary


def ConvertToSequence(parent_list: List[str],
                      child_list: List[str],
                      parent_arr: np.ndarray,
                      child_arr: np.ndarray,
                      residue_dict_flip: Dict[int, str],
                      ) -> Dict[str, str]:

    new_seq = {}
    for i in range(len(parent_list)):
        parent, child = parent_list[i], child_list[i]

        new_seq[parent] = "".join([residue_dict_flip[num] if num != 0 else "-" for num in parent_arr[i, :]])
        new_seq[child] = "".join([residue_dict_flip[num] if num != 0 else "-" for num in child_arr[i, :]])

    return new_seq

