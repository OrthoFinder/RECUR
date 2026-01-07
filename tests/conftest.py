import os
import pytest

@pytest.fixture()
def test_data():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    test_data_path = os.path.join(
        repo_root, "ExampleData", "example_alignments.aln.recur.tsv"
    )

    if not os.path.isfile(test_data_path):
        pytest.skip(f"Test data file not found: {test_data_path}")

    return test_data_path


@pytest.fixture()
def get_recurrence_list(test_data):
    recurrence_list = []
    with open(test_data, "r", encoding="utf-8") as reader:
        for line in reader:
            if line.startswith("Site"):
                continue
            fields = line.rstrip("\n").split("\t")
            fields = fields[:5] + fields[7:]
            recurrence_list.append(" ".join(fields))
    return recurrence_list
