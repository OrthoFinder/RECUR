---
layout: page
title: Usage
permalink: /usage/
---

## How to Use RECUR

- [Simple Usage](#simple-usage)
- [Advance Usage](#advanced-usage)
  - [Options Overview](#options-overview)
  - [Using a constraint tree](#using-a-constraint-tree)
  - [Providing a model of evolution](#providing-a-model-of-evolution)
  - [Running RECUR on a directory](#running-recur-on-a-directory)
- [RECUR Results](#recur-results)
  - [Recurrence List](#recurrence-list)
  - [Output Structure and Extended Results](#output-structure-and-exteneded-results)
- [A Note on Reproducibility](#a-note-on-reproducibility)
- [References](#references)

### Simple Usage

The minimal requirements of RECUR is a MSA (protein or codon) in FASTA format with the sequence type specified and a defined outgroup species or clade. e.g.,

`recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON>`

* `--outgroups`: informs RECUR how to correctly root the tree.

   You can either provide a .txt file with each outgroup species listed on a new line, or if you have a small number of outgroup species you can write the species on the command line. e.g., `--outgroups "SpeciesA,SpeciesB,SpeciesC"`.

* `-st`: signals the sequence type in the MSA.

   For a protein MSA `-st AA` should be provided. For a codon MSA, the different NCBI genetic codes can be specified (defined [here](http://www.iqtree.org/doc/Substitution-Models#codon-models)). `-st CODON` will use the standard genetic code.


<!-- ### Options Overview -->

### Advanced Usage

#### Options Overview 
In this section, we will dive deep into the options you have to run RECUR. The commands shown in this section will assume that you have RECUR installed on your machine.

```bash
  -f <dir/file>                Protein or codon alignment in FASTA format [Required]
  -st <str>                     <AA|CODON> [Required][Default: AA]
  --outgroups <dir/file/str>   List of outgroup sequences [Required]
  --num-alignments <int>       Number of simulated alignments for p-value estimation [Default: 1000]
  -te <dir/file>               Complete constraint tree [Default: Estimated from alignment]
  -m <str>                     Model of sequence evolution [Default: estimated from alignment]
  -nt <int>                    Number of threads provided to IQ-TREE
  -t <int>                     Number of threads used for RECUR internal processing
  --seed <int>                 Random starting see number [Default: 8]
  -o <txt>                     Results directory [Default: same directory as MSA files]
  -uc <int>                    Update cycle used in progress bar [Default: no progress bar]
  -bs <int>                    Batch size used in subsitution analysis of the Monte Carlo Simulated sequences [Default: no batch processing]
  -iv <str>                    IQ-TREE version [Default: iqtree2]
  -blfix                       Fix branch lengths of tree. [Default: False]
```

Please note that the default values for `-t`, `-nt` are processor dependent. If you are following the installation step mentioned in the previous section, you can run one of the following commands to find out the actual default setting for your machine.

```bash
recur
python3 recur.py
```

#### Using a constraint tree

To specify the topology of the phylogeny used by RECUR the user can provide a constraint tree using the -te flag. The argument is a file containing a tree in Newick format. E.g.,

```bash
recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON> -te <treefile>
```

#### Providing a model of evolution

A model of sequence evolution (as long as it is supported by IQ-TREE) can be provided using the `-m` flag. E.g.,

```bash
recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON> -te <treefile> -m <model_of_evolution>
```

#### Running RECUR on a directory

To free you from typing multiple commands in a terminal to run on multiple genes, RECUR provides an option to run on a folder. E.g.,

```bash
recur [options] -f <directory> --outgroups <directory> -st <AA|CODON> -te <directory>
```
For example, if you have three genes files, each have a different set of outgroups and tree files. You can place those outgroups files and tree files with the corresponding gene names inside the same MSA data folder. E.g.,
<p align="center">
  <img src="../assets/images/RECUR_input_structure_1.PNG" alt="RECUR input structure 1" width="250"/>
</p>

If your genes share the same outgroups and tree, you only need to create a single outgroups file and tree file in your data folder. Those two files will be share by all the genes for your RECUR analysis. E.g.,

<p align="center">
  <img src="../assets/images/RECUR_input_structure_2.PNG" alt="RECUR input structure 2" width="250"/>
</p>

> **Important Information when running RECUR on a directory**:
>  * `<alignment_file>`: can only have `.aln`, `fasta`, `faa`, or `fa` as the file extensions. (When running RECUR directly on an alignment file, extension requirement does not apply.)
>  * `<outgroups_file>`: needs to have `.outgroup` in the file name.
>  * `<treefile>`: needs to have `.tree` in the file name.

- Running RECUR with parallel processing

To speed up the analysis, RECUR has introduced both multiprocessing and multithreading in the package. By default, the number of threads will be automatically calculated based the configuration of your machine. Alternatively, they can be user specified using the following three options.

- `-t`: Number of threads used for the RECUR algorithms
- `-nt`: Number of threads for IQ-TREE to run in parallel

Note, the number of threads cannot exceed the computer allowance.

> * `-t` <= maximum_num_logical_threads
> * `-nt` <= maximum_num_logical_threads

For instance, if your computer has 8 logical threads, you can have the following setup,

```bash
recur -f example_alignments.aln -st AA --outgroups example_alignments.outgroups.txt -t 8 -nt 8
```

### RECUR Results
#### Recurrence List

The main output file from RECUR is a list of all recurrent substitutions inferred to have occurred across the phylogeny, which is found in the `*.recur.tsv` file.

<p align="center">
  <img src="../assets/images/RECUR_recurrence_list.PNG" alt="RECUR recurrence list" width="500"/>
</p>

The `*.recur.tsv` file contains a list of recurrent substitutions, i.e., an amino acid substitution that occurred more than once at a specific site, identified in the alignment. Here is the breakdown of each columns inside this output file:
* `Site`: indicates the position in the protein MSA.
* `Parent`: indicates the ancestral amino acid residue.
* `Child`: indicates the descent amino acid residue.
* `Recurrence`: indicates the number of times that amino acid substitution occurred at that site across the phylogeny.
* `Reversion`: indicates the number of times the reverse amino acid substitution occurred at that site across the phylogeny.
* `P-value`: signifies the probability that the substitution has not occurred by chance (for details of its calculation please read Robbins et al., 2025).

  Note: A p-value equal to 0.001 (with the default 1000 simulations) indicates the recurrence of the substitution was not observed in any of the simulated alignments and is unlikely to have occurred by chance.
* `AllSiteSubs`: provides a list of all other substitutions that occurred at that site.
* `SiteComposition`: indicates the count of each residue at that site in the protein MSA.


#### Output Structure and Extended Results

RECUR outputs additional results files from each step in the analysis which can be found in the `.recur` folder. Provided no output folder was specified using `-o`, the results are structured as seen below. 
<!-- If you do not specify the output folder using `-o`, the results are located in the same directory as the MSA files. -->

<p align="center">
  <img src="../assets/images/RECUR_output_structure.PNG" alt="RECUR output structure" width="500"/>
</p>

 <!-- Additional results of intermediate steps can be found in the `.recur` folder. -->

Inside the `*.recur` folder, you will find two `*.txt` files, i.e., `Citation.txt` and `Log.txt`, and four subdirectories, i.e., `Real_Phylogeny`, `Infered_Sequences`, `Monte_Carlo_Simulation` and `Substitution_Matirces`.

- `Real_Phylogeny`: contains the tree constructed from the alignment (`*.treefile`), the associated .iqtree file (`*.iqtree`) with information on model testing and the raw ancestral state reconstruction file (`*.state`) output from IQ-TREE.
- `Inferred_Sequences`: contains the inferred ancestral sequences with node labels corresponding to the tree in the Real_Phylogeny directory.
- `Substitution_Matirces`: contains a `*.substituion_matrix.tsv` file and a `*.accum_substituion_matrix.txt` file. The `*.substituion_matrix.tsv` summarises the substitutions that occurred at each site in the protein across the phylogeny, while the `*.accum_substituion_matrix.txt` presents the accumulated counts for each amino acid substitution type across all residues.

For example, by running RECUR on the `example_alignments.aln` file, the `example_alignments.aln.substituion_matrix.tsv` would have the following content:
<p align="center">
  <img src="../assets/images/RECUR_substitution_matrices.PNG" alt="RECUR substitution matrices" width="300"/>
</p>

`Site` represents the position in the protein alignment. Please note that the sites of an alignment **start from 1** in RECUR. The `Parent>Child:MutCount` column contains a list of amino acid substitutions that occurred at that site. For example, at site 2 a lysine (K) mutated to a serine (S) twice. A dash signifies no substitution was inferred at that site.

As for the `RowIndex` and `ColIndex` columns, they index the substitution in a mutation matrix.  The arrangement of the mutation matrix is the same as that in the accumulated subsititon matrix. Here, the accumulated subsititon matrix is obtained by summing up all the mutation matrix at each site.

<p align="center">
  <img src="../assets/images/RECUR_accum_substitution_matrices.PNG" alt="RECUR accum substitution matrices" width="700"/>
</p>

- `Monte_Carlo_Simulation`: contains the simulated alignments based on the best model of evolution,  phylogeny and root sequence of the subtree excluding the outgroup. The number of alignment files in this directory will match the number of Monte Carlo Simulation used in your analysis. It can be adjusted by setting `--num-alignments` to a different value.

### A Note on Reproducibility

According to Shen et al. (2020), irreproducibility in maximum likelihood phylogenetic inference is a significant issue [^1]. Regardless of the method related parameters that would affect the reproducibility of the final outputs, different random starting seed number, number of threads and processor type can also introduce uncertainties to the results. Such effect can be observed by setting different `-nt`, or setting different seed number using `--seed` in RECUR for each run. When running RECUR on different machine with different processor type even with the fixed `-nt` and `--seed`, the results can still be different. Nevertheless, the results should stay the same for each run when running RECUR on the same machine with the fixed input parameters.

## References

[^1]: *Shen, XX., Li, Y., Hittinger, C.T. et al.* **An investigation of irreproducibility in maximum likelihood phylogenetic inference.** Nat Commun 11, 6096 (2020). [![DOI:10.1038/s41467-020-20005-6](https://raw.githubusercontent.com/OrthoFinder/RECUR/main/docs/images/doi-badge.svg)](https://doi.org/10.1038/s41467-020-20005-6)