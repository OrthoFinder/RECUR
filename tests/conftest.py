import os 
import sys
import pytest
import shutil
import pathlib
from typing import Tuple 
import multiprocessing as mp
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

@pytest.fixture()
def test_data():
    cwd = os.getcwd()
    test_data_path = os.path.join(cwd, "tests", "benchmark", "pruned_random_tree_15.aln.recfinder.tsv")
    if not os.path.exists(test_data_path):
        print("Test file does not exist!")
        sys.exit()
    return test_data_path

@pytest.fixture()
def get_recurrence_list(test_data):
    recurrance_list = []
    with open(test_data) as reader:
        for line in reader:
            if "Site" not in line:
                line = line.strip().split("\t")[:5] + line.strip().split("\t")[6:]
                recurrance_list.append(" ".join(line))
                print(" ".join(line))
    return recurrance_list







