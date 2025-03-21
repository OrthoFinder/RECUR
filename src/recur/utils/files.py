#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import traceback
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple, Union

import dendropy
import numpy as np
from rich import print

from recur import __version__
from recur.utils import process_args, util


class FileHandler(object):

    def __init__(self, results_dir: str = ""):
        self.wd_current = ""
        self.rd1 = results_dir
        self.gene_of_interest = None

    def CreateOutputDirectories(self, options: process_args.Options, base_dir: str) -> None:

        if options.name == "":
            self.rd1 = util.CreateNewWorkingDirectory(base_dir,
                                                      options.sequence_type,
                                                      extended_filename=options.extended_filename,
                                                      keepprev=options.keepprev)
        else:
            self.rd1 = util.CreateNewWorkingDirectory(base_dir + "_" + options.name,
                                                      options.sequence_type,
                                                      qDate=False,
                                                      extended_filename=options.extended_filename,
                                                      keepprev=options.keepprev)
        self.wd_current = self.rd1 + os.sep
        if not os.path.exists(self.wd_current):
            os.makedirs(self.wd_current, exist_ok=True)

    def CreateMCSDirectories(self, options: process_args.Options) -> None:
        mcs_dir = self.rd1 + "Monte_Carlo_Similation"
        
        if options.disk_save:
            mcs_dir = util.CreateNewWorkingDirectory(mcs_dir,
                                                    options.sequence_type,
                                                    qDate=False,
                                                    extended_filename=False,
                                                    keepprev=True)
        self.mcs_dir = mcs_dir 
        if not os.path.exists(self.mcs_dir):
            os.makedirs(self.mcs_dir, exist_ok=True)


    def CheckFileCorrespondance(self, gene: str, state_file: str, tree_file:str) -> None:
        gene_statfile = os.path.basename(state_file).rsplit(".", 1)[0]
        gene_treefile = os.path.basename(tree_file).rsplit(".", 1)[0]
        if gene != gene_statfile and gene !=  gene_treefile:
            print('ERROR: Gene files have not been correctly processed')
            util.Fail()

    @staticmethod
    def ProcessesNewAln(alnDir: str, alnPath: Optional[str] = None) -> Dict[str, str]:
        originalALNFilenames = {}

        if alnPath is None:
            alnExtensions = {"aln", "fasta", "fa", "faa"} # need to add more potential file extensions
            if not os.path.exists(alnDir):
                print("\nDirectory does not exist: %s" % alnDir)
                util.Fail()
            files_in_directory = sorted([f for f in os.listdir(alnDir) if os.path.isfile(os.path.join(alnDir, f))])
            excludedFiles = []
            for f in files_in_directory:
                if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in alnExtensions and not f.startswith("._"):
                    gene_of_interest = f.rsplit(".", 1)[0]
                    originalALNFilenames[gene_of_interest] = os.path.join(alnDir, f)
                else:
                    excludedFiles.append(os.path.join(alnDir, f))

            if len(excludedFiles) != 0:
                print("\nWARNING: Files have been ignored as they don't appear to be ALN files:")
                for f in excludedFiles:
                    print(f"    {f}")
                print("RECUR expects ALN files to have one of the following extensions: %s" % (", ".join(alnExtensions)))

            if len(originalALNFilenames) == 0:
                print("\nNo fasta files found in supplied directory: %s" % alnDir)
                util.Fail()

        else:
            f = os.path.basename(alnPath)
            if "." in f:
                if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower():
                    gene_of_interest = f.rsplit(".", 1)[0]
                    originalALNFilenames[gene_of_interest] = alnPath
            else:
                originalALNFilenames[f] = alnPath

        return originalALNFilenames

    @staticmethod
    def UpdateTreeFile(treefile: str) -> None:
        try:
            with open(treefile, 'r') as f:
                t = dendropy.Tree.get(file=f, schema="newick")

            for idx, node in enumerate(t.preorder_node_iter()):
                if node.is_leaf():
                    if node.taxon.label and node.taxon.label.isspace():
                        node.taxon.label = "_".join(node.taxon.label.split()).strip("_")
                else:
                    node.label = f"Node{idx}"

            with open(treefile, "w") as f:
                f.write(t.as_string(schema="newick"))

        except Exception as e:
            print(f"ERROR in updating labels to the treefile {treefile}: {e}")

    def GetBinaryStateFileFN(self) -> str:
        if self.GetBinaryPhylogenyDir() is None:
            raise Exception("No Indel_Estimation directory.")
        statefile = [os.path.join(self.GetBinaryPhylogenyDir(), file) for file in os.listdir(self.GetBinaryPhylogenyDir()) if file.endswith(".state")][0]
        return statefile

    def GetBinaryTreeFileFN(self) -> str:
        if self.GetBinaryPhylogenyDir() is None:
            raise Exception("No Indel_Estimation directory.")
        treefile = [os.path.join(self.GetBinaryPhylogenyDir(), file) for file in os.listdir(self.GetBinaryPhylogenyDir()) if file.endswith(".treefile")][0]
        return treefile

    def GetBinarySeqsFN(self) -> str:
        if self.GetBinaryPhylogenyDir() is None:
            raise Exception("No Indel_Estimation directory.")
        binary_seqs_file = os.path.join(self.GetBinaryPhylogenyDir(), f"{self.gene_of_interest}.binary_sequences.aln")
        return binary_seqs_file

    def GetBinaryNodeSeqsFN(self) -> str:
        if self.GetBinaryPhylogenyDir() is None:
            raise Exception("No Indel_Estimation directory.")
        node_seqs_file = os.path.join(self.GetBinaryPhylogenyDir(), f"{self.gene_of_interest}.inferred_binary_sequences.aln")
        return node_seqs_file

    def GetStateFileFN(self) -> str:
        if self.GetRealPhylogenyDir() is None:
            raise Exception("No Real_Phylogeny directory.")
        statefile = [os.path.join(self.GetRealPhylogenyDir(), file) for file in os.listdir(self.GetRealPhylogenyDir()) if file.endswith(".state")][0]
        return statefile

    def GetTreeFileFN(self) -> str:
        if self.GetRealPhylogenyDir() is None:
            raise Exception("No Real_Phylogeny directory.")
        treefile = [os.path.join(self.GetRealPhylogenyDir(), file) for file in os.listdir(self.GetRealPhylogenyDir()) if file.endswith(".treefile")][0]
        return treefile

    def GetIQTreeFileFN(self) -> str:
        if self.GetRealPhylogenyDir() is None:
            raise Exception("No Real_Phylogeny directory.")
        iqtreefile = [os.path.join(self.GetRealPhylogenyDir(), file) for file in os.listdir(self.GetRealPhylogenyDir()) if file.endswith(".iqtree")][0]
        return iqtreefile

    def GetMutCountMatricesFN(self) -> str:
        if self.GetMutMatrixDir() is None:
            raise Exception("No Substitution_Matrices directory.")
        mut_matrix_file = os.path.join(self.GetMutMatrixDir(), f"{self.gene_of_interest}.substituion_matrix.tsv")
        return mut_matrix_file

    def GetAccumMutCountMatricesFN(self) -> str:
        if self.GetMutMatrixDir() is None:
            raise Exception("No Substitution_Matrices directory.")
        mut_matrix_file = os.path.join(self.GetMutMatrixDir(), f"{self.gene_of_interest}.accum_substituion_matrix.txt")
        return mut_matrix_file

    def GetCombinedDNASeqsFN(self) -> str:
        if self.GetInferredSeqsDir() is None:
            raise Exception("No Inferred_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferredSeqsDir(), f"{self.gene_of_interest}.combined_dna_sequences.aln")
        return combined_seqs_file

    def GetNodeDNASeqsFN(self) -> str:
        if self.GetInferredSeqsDir() is None:
            raise Exception("No Inferred_Sequences directory.")
        node_seqs_file = os.path.join(self.GetInferredSeqsDir(), f"{self.gene_of_interest}.inferred_ancestral_sequences.dna.aln")
        return node_seqs_file

    def GetBinaryCombinedProtSeqsFN(self) -> str:
        if self.GetInferredSeqsDir() is None:
            raise Exception("No Inferred_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferredSeqsDir(), f"{self.gene_of_interest}.recur_combined_protein_sequences_bin.aln")
        return combined_seqs_file

    def GetCombinedProtSeqsFN(self) -> str:
        if self.GetInferredSeqsDir() is None:
            raise Exception("No Inferred_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferredSeqsDir(), f"{self.gene_of_interest}.recur_combined_protein_sequences.aln")
        return combined_seqs_file

    def GetNodeProtSeqsFN(self) -> str:
        if self.GetInferredSeqsDir() is None:
            raise Exception("No Inferred_Sequences directory.")
        node_seqs_file = os.path.join(self.GetInferredSeqsDir(), f"{self.gene_of_interest}.inferred_ancestral_sequences.prot.aln")
        return node_seqs_file

    def GetRecurrenceListFN(self, recDir: Optional[str] = None) -> str:
        if recDir is None:
            raise Exception("No Results directory.")

        recurrence_list_file = os.path.join(recDir, f"{self.gene_of_interest}.recur.tsv")
        return recurrence_list_file
    
    def GetRecurrenceListRealPhylogenyFN(self, recDir: Optional[str] = None) -> str:
        if recDir is None:
            raise Exception("No Results directory.")

        recurrence_list_file = os.path.join(recDir, f"{self.gene_of_interest}.recur.real_phylogeny.tsv")
        return recurrence_list_file

    def GetRecurrenceCountPhylogenyFN(self, recDir: Optional[str] = None) -> str:
        if recDir is None:
            raise Exception("No Results directory.")

        recurrence_count_file = os.path.join(recDir, f"{self.gene_of_interest}.recur.real_phylogeny.count.tsv")
        return recurrence_count_file

    def GetMCSimulationTreeFN(self) -> str:
        if self.mcs_dir is None:
            raise Exception("No Monte_Carlo_Simuation directory.")
        mcs_treefile = [os.path.join(self.mcs_dir, file) for file in os.listdir(self.mcs_dir) if file.endswith("alisim.full.treefile")]

        if len(mcs_treefile) == 0:
            raise Exception("No Monte Carlo Simulated treefile.")

        return mcs_treefile[0]

    def GetMutMatrixDir(self) -> str:
        d = os.path.join(self.rd1, "Substitution_Matrices") + os.sep
        # self.rd1 + "Substitution_Matrices" + os.sep
        if not os.path.exists(d):
            os.mkdir(d)
        return d

    def GetInferredSeqsDir(self) -> str:
        d = os.path.join(self.rd1, "Inferred_Sequences") + os.sep
        # d = self.rd1 + "Inferred_Sequences" + os.sep
        if not os.path.exists(d):
            os.mkdir(d)
        return d

    # def GetMCSimulationDir(self) -> str:
    #     d = self.rd1 + "Monte_Carlo_Similation" + os.sep
    #     if not os.path.exists(d):
    #         os.mkdir(d)
    #     return d

    def GetBinaryPhylogenyDir(self) -> str:
        d = os.path.join(self.rd1, "Indel_Estimation") + os.sep
        # d = self.rd1 + "Indel_Estimation" + os.sep
        if not os.path.exists(d):
            os.mkdir(d)
        return d

    def GetRealPhylogenyDir(self) -> str:
        d = os.path.join(self.rd1, "Real_Phylogeny") + os.sep
        # d = self.rd1 + "Real_Phylogeny" + os.sep
        if not os.path.exists(d):
            os.mkdir(d)
        return d

    def GetResultsDirectory(self) -> str:
        if self.rd1 is None: raise Exception("No rd1")
        return self.rd1


class FileReader(object):

    @staticmethod
    def ReadAlignment(fn: str) -> Tuple[Dict[str, str], int, bool]:
        msa: Dict[str, str] = dict()
        accession = ""
        seq = ""
        seq_len_set: Set[int] = set()
        dash_exist = False
        with open(fn, 'r') as infile:
            for line_ in infile:
                line = line_.replace("\n", "").strip()
                if line.startswith(">") or ">" in line:
                    if accession:
                        msa[accession] = seq
                        if len(seq_len_set) != 0:
                            if seq_len_set.pop() != len(seq):

                                print(f"The sequence length in {fn} is inconsistent.")
                                sys.stderr.flush()
                                print(traceback.format_exc())
                                sys.exit(1)

                        seq_len_set.add(len(seq))

                    # accession = line[1:].replace('_', ' ')
                    accession = line[1:]
                    if ">" in line:
                        accession = accession.replace('>', '')
                    seq = ""
                else:
                    if "-" in line:
                        dash_exist = True

                    if len(set(line) - set(util.legal_chars)):
                        print(f"RECUR failed due to invalid line found in {fn}.")
                        print(line)
                        print("Please remove any special characters from the alignment file, if present, and try again.")
                        sys.stderr.flush()
                        print(traceback.format_exc())
                        sys.exit(1)

                    seq += line
            if accession:
                msa[accession] = seq
                if "-" in seq:
                    dash_exist = True
                if seq_len_set.pop() != len(seq):
                    print(f"The sequence length in {fn} is inconsistent.")
                    sys.stderr.flush()
                    print(traceback.format_exc())
                    sys.exit(1)

        return msa, len(seq), dash_exist

    @staticmethod
    def ReadStateFile(fn: str) -> Dict[str, str]:

        node_seq_dict = {}
        with open(fn, 'r') as reader:
            for line in reader:
                if "#" in line or "Site" in line:
                    continue
                line_list = line.replace("\n", "").strip().split("\t", 3)
                node, state = line_list[0], line_list[2]
                if node not in node_seq_dict:
                    node_seq_dict[node] = state
                elif node in node_seq_dict:
                    node_seq_dict[node] += state

        return node_seq_dict

    @staticmethod
    def ReadIQTreeFile(iqtree_path: str) -> str:
        with open(iqtree_path) as reader:
            for line in reader:
                if "Best-fit model" in line or "Model of substitution" in line:
                    best_evolution_model = line.replace("\n", "").strip().split(": ")[-1]

        return best_evolution_model

    @staticmethod
    def ConverToBinary(fn: str) -> Dict[str, str]:
        msa_binary = {}
        accession = ""
        seq = ""
        with open(fn, 'r') as infile:
            for line_ in infile:
                line = line_.replace("\n", "").strip()
                if line.startswith(">") or ">" in line:
                    if accession:
                        msa_binary[accession] = seq
                    # accession = line[1:].replace('_', ' ')
                    accession = line[1:]
                    if ">" in line:
                        accession = accession.replace('>', '')
                    seq = ""
                else:

                    binary_line = [
                        1 if l in util.residues else 0 if l in util.reserved_chars else 2 for l in line
                    ]

                    if 2 in binary_line:
                        print(f"RECUR failed due to invalid line found in {fn}.")
                        print(line)
                        print("Please remove any special characters from the alignment file, if present, and try again.")
                        sys.stderr.flush()
                        print(traceback.format_exc())
                        sys.exit(1)
                    seq += ''.join(map(str, binary_line))

            if accession:
                msa_binary[accession] = seq

        return msa_binary
    
    @staticmethod    
    def ReadMCSRecurrenceCount(base_dir: str) -> List[Dict[Tuple[int, int, int], int]]:
        mcs_count_files = []
        with os.scandir(base_dir) as entries:
            for entry in entries:
                if entry.is_dir():
                    if "Monte_Carlo_Similation" in entry.name:
                        for file in os.listdir(entry):
                            file_path = os.path.join(entry, file)
                            mcs_count_files.append(file_path)

        mcs_results = []
        for file in mcs_count_files:
            mcs_count_dict = {}
            with open(file) as reader:
                for i, line in enumerate(reader):
                    if i == 0:
                        continue

                    res_loc, parent, child, count = map(int, line.strip().split("\t"))
                    mcs_count_dict[(res_loc, parent, child)] = count
            mcs_results.append(mcs_count_dict)

        return mcs_results
    
    @staticmethod    
    def ReadRecurrenceCount(recurrence_count_file: str) -> Dict[Tuple[int, int, int], int]:

        recurrence_count_dict = {}
        with open(recurrence_count_file) as reader:
            for i, line in enumerate(reader):
                if i == 0:
                    continue
                res_loc, parent, child, count = map(int, line.strip().split("\t"))
                recurrence_count_dict[(res_loc, parent, child)] = count

        return recurrence_count_dict
    
    @staticmethod    
    def ReadRecurrenceList(recurrence_list_fn: str) -> List[List[Union[str, int]]]:

        recurrence_list = []
        with open(recurrence_list_fn) as reader:
            for i, line in enumerate(reader):
                if i == 0:
                    continue
                res_loc, parent, child, recurrence, reversion = line.strip().split("\t")
                recurrence_list.append([
                    int(res_loc) - 1,
                    parent,
                    child,
                    int(recurrence),
                    int(reversion)
                ])
        
        return recurrence_list
    

class FileWriter(object):

    @staticmethod
    def WriteMutMatrix(res_loc_count_dict: Dict[Tuple[int, int, int], int],
                        residue_dict_flip: Dict[int, str],
                        protein_len: int,
                        mut_matrix_path: str,
                        accum_matrix_path: str) -> None:

        res_loc_info_dict = util.get_sorted_res_loc_info(res_loc_count_dict, protein_len)
        colname = ["Site", "Parent>Child:SubCount", "RowIndex", "ColIndex"]
        accum_mutation_matrix = np.zeros((20, 20), dtype=int)

        with open(mut_matrix_path, "w") as writer:
            writer.write("\t".join(colname) + "\n")

            for pos, item in sorted(res_loc_info_dict.items()):
                if len(item) == 0:
                    line = "\t".join((str(pos+1), "-", "-", "-"))
                else:
                    parent_ids, child_ids, counts = map(np.array, zip(*item))
                    np.add.at(accum_mutation_matrix, (parent_ids - 1, child_ids - 1), counts)

                    parent_child = []
                    rows = []
                    cols = []
                    for parent_id, child_id, data in item:
                        parent = residue_dict_flip[parent_id]
                        child = residue_dict_flip[child_id]
                        parent_child.append(">".join((parent, child)) + ":" + str(data))
                        rows.append(parent_id)
                        cols.append(child_id)

                    parent_child_str = ",".join(parent_child)
                    row_str = ",".join(map(str, np.array(rows)))
                    col_str = ",".join(map(str, np.array(cols)))
                    line = "\t".join((str(pos+1), parent_child_str, row_str, col_str))

                line += "\n"
                writer.write(line)

        with open(accum_matrix_path, "w") as writer:
            col_index = "    ".join(util.residues)
            col_index = "          ".join((" ", col_index)) + "\n"
            col_val0 = "    ".join(map(str, range(1, 10)))
            col_val1 = "   ".join(map(str, range(10, 21)))
            col_val = col_val0 + "    " + col_val1
            col_val = "          ".join((" ", col_val)) + "\n"
            writer.write(col_val)
            writer.write(col_index)

            for i in range(accum_mutation_matrix.shape[0]):
                row_str = "    ".join(map(str, accum_mutation_matrix[i]))
                row_str = "    ".join((util.residues[i], row_str))
                if i + 1 < 10:
                    row_str = "     ".join((str(i+1), row_str))
                else:
                    row_str = "    ".join((str(i+1), row_str))
                row_str += "\n"
                writer.write(row_str)

    @staticmethod
    def WriteSeqsToAln(seqs_dict: Dict[str, str], outFilename: str) -> None:
        with open(outFilename, 'w') as writer:
            for seq in seqs_dict:
                writer.write(">%s\n" % seq)
                writer.write(seqs_dict[seq] + "\n")
                writer.write("\n")

    @staticmethod
    def WriteMCSRecurrenceCountToFile(
        mcs_results: List[Dict[Tuple[int, int, int], int]],
        mcs_dir: str,
    ) -> None:
        
        colname = ["Site", "ParentIndex", "ChildIndex", "Count"]
        for i, mcs_count_dict in enumerate(mcs_results):
            outFilename = os.path.join(mcs_dir, f"mcs_recurrence_count_{i}.tsv")
            with open(outFilename, "w") as writer:
                writer.write("\t".join(colname) + "\n")
                for key, val in mcs_count_dict.items():
                    rec_loc, parent, child = key 
                    writer.write("\t".join((str(rec_loc), str(parent), str(child), str(val))) + "\n")

    @staticmethod
    def WriteRecurrenceCountToFile(
        rec_loc_count_dict: Dict[Tuple[int, int, int], int],
        outFilename: str,
    ) -> None:
        
        colname = ["Site", "ParentIndex", "ChildIndex", "Count"]
        with open(outFilename, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for key, val in rec_loc_count_dict.items():
                rec_loc, parent, child = key 
                writer.write("\t".join((str(rec_loc), str(parent), str(child), str(val))) + "\n")


    @staticmethod
    def WriteRecurrenceListRealPhylogeny(
            recurrence_list: List[List[Union[str, int]]],
            outFilename: str,
        ) -> None:
        colname = ['Site', 'Parent', 'Child', 'Recurrence', "Reversion"]
        with open(outFilename, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for rec_list in recurrence_list:
                if isinstance(rec_list[0], int):
                    rec_list[0] += 1
                writer.write("\t".join(map(str, rec_list)) + "\n")

    @staticmethod
    def WriteRecurrenceList(recurrence_list: List[List[Union[str, int, float]]],
                            outFilename: str,
                            options: process_args.Options,
                            ) -> None:

        colname = ['Site', 'Parent', 'Child', 'Recurrence',
                   "Reversion", "P-Value", "AllSiteSubs",  "SiteComposition"]

        baseFileName = outFilename.rsplit(".", 1)[0] 
        if options.name == "":
            outFilename = util.CreateNewFileName(baseFileName,
                                            options.sequence_type,
                                            extended_filename=options.extended_filename,
                                            keepprev=options.keepprev)
        else:
            outFilename = util.CreateNewFileName(baseFileName + "." + options.name,
                                                options.sequence_type,
                                                qDate=False,
                                                extended_filename=options.extended_filename,
                                                keepprev=options.keepprev)
 
        with open(outFilename, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for rec_list in recurrence_list:
                if isinstance(rec_list[0], int):
                    rec_list[0] += 1
                writer.write("\t".join(map(str, rec_list)) + "\n")
