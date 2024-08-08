#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import scipy as spy 
import numpy as np
from typing import Dict, List, Optional, Union
from recur.utils import util, process_args
from recur import __version__
import dendropy


class FileHandler(object):  

    def __init__(self):
        # self.base_matrix_format = "MUT%07d"
        self.wd_base = []               # Base: blast, species & sequence IDs, species fasta files - should not request this and then write here
        self.wd_current = ""          # Location to write out any new files
        self.rd1 = "" 
        self.species_ids_corrected = None
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
        self.wd_current = self.rd1
        if not os.path.exists(self.wd_current):
            os.makedirs(self.wd_current, exist_ok=True)
        self.wd_base = [self.wd_current]  


    def CheckFileCorrespondance(self, gene: str, state_file: str, tree_file:str) -> None:
        gene_statfile = os.path.basename(state_file).split(".", 1)[0]
        gene_treefile = os.path.basename(tree_file).split(".", 1)[0]
        if gene != gene_statfile and gene !=  gene_treefile:
            print('ERROR: Gene files have not been correctly processed')
            util.Fail()
    
    @staticmethod
    def ProcessesNewAln(alnDir: str, alnPath: Optional[str] = None) -> Dict[str, str]:
        originalALNFilenames = {}
        
        if alnPath is None:
            alnExtensions = {"aln", "fasta"} # need to add more potential file extensions
            if not os.path.exists(alnDir):
                print("\nDirectory does not exist: %s" % alnDir)
                util.Fail()
            files_in_directory = sorted([f for f in os.listdir(alnDir) if os.path.isfile(os.path.join(alnDir, f))])
            excludedFiles = []
            for f in files_in_directory:
                if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in alnExtensions and not f.startswith("._"):
                    if len(f.split(".", 1)) == 2:
                        gene_of_interest = f.split(".", 1)[0]
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
            if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower():
                if len(f.split(".", 1)) == 2:
                    gene_of_interest = f.split(".", 1)[0]
                    originalALNFilenames[gene_of_interest] = alnPath
        
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
        mut_matrix_file = os.path.join(self.GetMutMatrixDir(), f"{self.gene_of_interest}.accum_substituion_matrix.tsv")
        return mut_matrix_file 

    def GetCombinedDNASeqsFN(self) -> str:
        if self.GetInferedSeqsDir() is None: 
            raise Exception("No Infered_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.combined_dna_sequences.aln")
        return combined_seqs_file 

    
    def GetNodeDNASeqsFN(self) -> str:
        if self.GetInferedSeqsDir() is None: 
            raise Exception("No Infered_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.infered_ancestral_sequences.dna.aln")
        return combined_seqs_file 

    def GetCombinedProtSeqsFN(self) -> str:
        if self.GetInferedSeqsDir() is None: 
            raise Exception("No Infered_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.combined_protein_sequences.aln")
        return combined_seqs_file 

    
    def GetNodeProtSeqsFN(self) -> str:
        if self.GetInferedSeqsDir() is None: 
            raise Exception("No Infered_Sequences directory.")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.infered_ancestral_sequences.prot.aln")
        return combined_seqs_file 

    def GetRecurrenceListFN(self, recDir: Optional[str] = None) -> str:
        if recDir is None: 
            raise Exception("No Results directory.")
        
        recurrence_list_file = os.path.join(recDir, f"{self.gene_of_interest}.recur.tsv")
        return recurrence_list_file


    def GetMCSimulationTreeFN(self) -> str:
        if self.GetMCSimulationDir() is None: 
            raise Exception("No Monte_Carlo_Simuation directory.")
        mcs_treefile = [os.path.join(self.GetMCSimulationDir(), file) for file in os.listdir(self.GetMCSimulationDir()) if file.endswith("alisim.full.treefile")]
        
        if len(mcs_treefile) == 0:
            raise Exception("No Monte Carlo Simulated treefile.")
        
        return mcs_treefile[0]

    def GetMutMatrixDir(self) -> str:
        d = self.rd1 + "Substitution_Matrices" + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetInferedSeqsDir(self) -> str:
        d = self.rd1 + "Infered_Sequences" + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationDir(self) -> str:
        d = self.rd1 + "Monte_Carlo_Similation" + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationFADir(self) -> str:
        d = self.GetMCSimulationDir() + "FA_Files" + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetRealPhylogenyDir(self) -> str:
        d = self.rd1 + "Real_Phylogeny" + os.sep
        if not os.path.exists(d): os.mkdir(d)
        return d
        
    def GetResultsDirectory(self) -> str:
        if self.rd1 is None: raise Exception("No rd1")
        return self.rd1
        

class FileReader(object):

    @staticmethod    
    def ReadAlignment(fn: str):
        msa: Dict[str, str] = dict()
        accession = ""
        seq = ""
        with open(fn, 'r') as infile:
            for line_ in infile:
                line = line_.replace("\n", "").strip()
                if line.startswith(">") or ">" in line:
                    if accession:
                        msa[accession] = seq
                    accession = line[1:].replace('_', ' ')
                    if ">" in line:
                        accession = accession.replace('>', '')
                    seq = ""
                else:
                    seq += line
            if accession:
                msa[accession] = seq
        return msa, len(seq)
    
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
    def ReadMutMatrixFile(mut_matrix_path: str) -> Dict[int, np.typing.NDArray[np.int64]]:

        matrix_dict ={}
        mut_matrix = np.zeros((20, 20), dtype=np.int64)
        with open(mut_matrix_path, 'r') as reader:
            for i, line in enumerate(reader):
                if "Site" in line:
                    continue
                line_list = line.replace("\n", "").strip().split()
                data = line_list[2]
                if data == "-":
                    matrix_dict[i] = mut_matrix

                else:
                    data_arr = np.array([*map(np.int64, data.split(","))])
                    row = np.array([*map(np.int64, line_list[3].split(","))])
                    col = np.array([*map(np.int64, line_list[4].split(","))])#
                    coo = spy.sparse.coo_matrix((data_arr, (row, col)), shape=(20, 20))
                    matrix = coo.toarray()
                    matrix_dict[i] = matrix

        return matrix_dict

    @staticmethod
    def ReadIQTreeFile(iqtree_path: str) -> str:
        with open(iqtree_path) as reader:
            for line in reader:
                if "Best-fit model" in line or "Model of substitution" in line:
                    best_evolution_model = line.replace("\n", "").strip().split(": ")[-1]
                
        return best_evolution_model


class FileWriter(object):
    
    @staticmethod
    def WriteMutMatrix(mut_matrices: Dict[int, np.typing.NDArray[np.int64]], 
                        residue_dict_flip: Dict[int, str], 
                        mut_matrix_path: str,
                        accum_matrix_path: str) -> None: 

        # fn_path = os.path.join(mut_matrix_dir, "mutation_matrices.tsv")
        mut_matricex_tuple = sorted([(pos, mut) for pos, mut in mut_matrices.items()])
        colname = ["Site", "Parent:Child:MutCount", "RowIndex", "ColIndex"]
        accum_mutation_matrix = np.zeros((20, 20), dtype=np.int64)
        with open(mut_matrix_path, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for pos, mut in mut_matricex_tuple:
                accum_mutation_matrix += mut
                coo = spy.sparse.coo_matrix(mut)
                 
                if len(coo.data) == 0:
                    line = "\t".join((str(pos+1), "-", "-", "-"))
                    
                else:
                    parent_child = []
                    for row, col, data in zip(coo.row, coo.col, coo.data):
                        parent = residue_dict_flip[row]
                        child = residue_dict_flip[col]
                        parent_child.append(">".join((parent, child)) + ":" + str(data))
                    parent_child_str = ",".join(parent_child)
                    row_str = ",".join(map(str, coo.row))
                    col_str = ",".join(map(str, coo.col))
                    line = "\t".join((str(pos+1), parent_child_str, row_str, col_str))
                line += "\n"
                writer.write(line)

        # accum_mut_fnpath = os.path.join(mut_matrix_dir, "accum_mutation_matrix.tsv")
        with open(accum_matrix_path, "w") as writer:
            col_index = "\t".join(util.residues)
            col_index = "\t".join((" ", col_index)) + "\n"
            writer.write(col_index)
            for i in range(accum_mutation_matrix.shape[0]):
                row_str = "\t".join(map(str, accum_mutation_matrix[i]))
                row_str = "\t".join((util.residues[i], row_str))
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
    def WriteRecurrenceList(recurrence_list: List[List[Union[str, int, float]]],
                            outFilename: str) -> None:

        colname = ['Site', 'Parent', 'Child', 'Recurrence', 
                   "Reversion", "P-Value", "AllSiteSubs",  "SiteComposition"]
        with open(outFilename, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for rec_list in recurrence_list:
                if isinstance(rec_list[0], int):
                    rec_list[0] += 1
                writer.write("\t".join(map(str, rec_list)) + "\n")



        