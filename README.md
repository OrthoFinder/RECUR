# RECUR ![Tested](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/OrthoFinder/RECUR/main/badge-status.json)

Finding Recurrent Substitutions from Multiple Sequence Alignments

## Introduction
![RECUR method workflow](./docs/images/RECUR_workflow_figure.png)

<div align="center">
  Figure 1: The RECUR workflow
</div>

The required input is either a protein or codon multiple sequence alignment (in FASTA format) and a defined outgroup species or clade. The output of RECUR is a list of recurrent amino acid substitutions, that have occurred in the inferred phylogeny (file suffix: `.recur.tsv`). Outputs of intermediate steps, i.e. model selection, tree inference, ancestral state reconstruction and site substitution matrices, can be found in the .recur output directory.

## Table of Contents

- [Getting started with RECUR](#getting-started-with-recur)
  - [Installing RECUR on Linux](#installing-recur-on-linux)
  - [Installing RECUR on Windows and MacOS](#installing-recur-on-windows-and-macos)
  - [Running RECUR in Conda](#running-recur-in-conda)
  - [Running RECUR in a Docker Container](#running-recur-in-a-docker-container)
- [How to Use RECUR](#how-to-use-recur)
  - [Options Overview](#options-overview)
  - [Simple Usage](#simple-usage)
  - [Advance Usage](#advanced-usage)
  - [Output Structure](#output-structure)
  - [Interpretation of the Recurrence List](#interpretation-of-the-recurrence-list)
  - [A Note on Reproducibility](#a-note-on-reproducibility)
- [Citations](#citations)
  - [Credits and Acknowledgements](#credits-and-acknowledgements)
- [References](#references)
- [Contributing](#contributing)

## Getting Started with RECUR
The recurrence analysis implemented by RECUR utilises IQ-TREE2 phylogenomic software package to infer a phylogeny, the ancestral node sequences and to generate simulated alignments.

### Installing RECUR on Linux

  If you are working on a Linux machine or WSL2, running RECUR is straightforward. A separate installation of IQ-TREE2 is not necessary.

  Before installing any relevant dependencies or packages, it is recommended that you create and activate a new virtual environment. You can do so by running:
  ```bash
  python3 -m venv .venv && . .venv/bin/activate
  ```
  *To deactivate the virtual environment* run
  ```bash
  deactivate
  ```
  1. Clone the repository
  ```bash
  git clone https://github.com/OrthoFinder/RECUR.git
  ```
  2. Navigate into the cloned repository
  ```bash
  cd RECUR
  ```
  Once you are inside the package's root directory, based on your requirements, you can choose either to install RECUR on your machine or just the required public python dependencies.
  * Option 1:  Install RECUR

    Running the following command will help you install RECUR in your current virtual environment.
    ```bash
    pip install .
    ```
    To test your installation, please run
    ```bash
    recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
    ```
  * Option 2: Install the requirements.txt without installing the package

    If you do not wish to install the RECUR package, you can simply run the following command to install the required dependencies.
    ```bash
    pip install -r requirements.txt
    ```
    Then run 
    ```bash
    python3 recur.py -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
    ```
  to test the installation and the environment. 
  
  > NOTE: if python3 does not work, please try python.

### Installing RECUR on Windows and MacOS 

  If you are on a Windows machine or a MAC, you need to install IQ-TREE2 before you can run RECUR. You can download the latest version of IQ-TREE2 [here](http://www.iqtree.org/#download).

  Once you have IQ-TREE2 installed, you can run the following commands based on your OS to create a virtual environment in which to run RECUR.

  - **Windows**

  If you are on windows, please open a command prompt, and run the following command to create and activate the virtual environment.
  ```bash
  python -m venv win_venv
  win_venv\Scripts\activate.bat
  ```

  *To deactivate the virtual environment* run
  ```bash
  win_venv\Scripts\deactivate.bat
  ```
  - **macOS**

  If you are on macOS, please run the following commands in a terminal.

  ```bash
  python3 -m venv mac_venv
  source mac_venv/bin/activate
  ```

  *To deactivate the virtual enviroment* run
  ```bash
  deactivate
  ```
  Now that you have successfully created and activated a new virtual environment, you can follow the three steps mentioned in [Installing RECUR on Linux](#installing-recur-on-linux) to install RECUR or just the python dependencies based on your needs. Then run one of the following commands (depending on whether you have installed RECUR or not) to test if RECUR runs on your machine.

  * With RECUR installed 
    ```bash
    recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt -iv system
    ```
  * Without RECUR installed 
    ```bash
    python3 recur.py -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt -iv system
    ```
  where `-iv` stands for IQ-TREE2 version. By default, RECUR will use the Linux binary version from the package. 

### Running RECUR in Conda

  Working in the conda environment may be easiest when you do not have access to a Linux machine. You can create and activate a new environment by running:

  ```bash
  conda create -n recur_env python=3.12
  conda activate recur_env
  ```
  Then install IQ-TREE2 package by running one of the following:
  ```bash
  conda install bioconda::iqtree
  conda install bioconda/label/cf201901::iqtree
  ```
  Now that you have a suitable environment with IQ-TREE2 installed, you can follow the previous steps to install RECUR directly or the relevant dependencies.
  * With RECUR installed 
    ```bash
    recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt -iv conda
    ```
  * Without RECUR installed 
    ```bash
    python3 recur.py -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt -iv conda
    ```

### Running RECUR in a Docker Container

  Apart from the above methods, you can also run RECUR inside a container. You can find the RECUR image on [Docker Hub](https://hub.docker.com/repository/docker/orthofinder/recur/general).

  - **Personal computer**
  
    Before you can run the RECUR container, you need to have Docker Desktop installed. 
    - Windows: https://docs/images.docker.com/desktop/install/windows-install/
    - macOS: https://docs/images.docker.com/desktop/install/mac-install/
    - Linux: https://docs/images.docker.com/desktop/install/linux-install/

    Once you have the Docker Desktop installed, please launch it and run the following command in the terminal to check if it is up and running.

    ```bash
    docker version
    ```
  - **Server**

    If you need to install the Docker Engine before you can run the RECUR container.

    Please find the right docker engine to install on your server in [here](https://docs/images.docker.com/engine/install/).

  With either Docker Desktop or Docker Engine installed on your machine, you can simply run the following command to test if you can run the RECUR container
  ```
  docker container run -it --rm orthofinder/recur:v1.0.0
  ```
  To run the RECUR container on your dataset, you will need to create a folder which contains your data in your current working directory. For instance, you have a data folder called MyData which contains a protein alignment file called `my_alignment.aln` and a file called `my_alignment.outgroups.txt` that contains all the outgroups, you can run the following command to start the RECUR container and make it run your dataset.
  ```
  docker container run -it --rm -v $(pwd)/MyData:/usr/src/RECUR/MyData:Z orthofinder/recur:v1.0.0 recur -f MyData/my_alignment.aln -st AA --outgroups MyData/my_alignment.outgroups.txt   
  ```
  Please note that arguments behind `orthofinder/recur:v1.0.0` will be the same as you run RECUR directly as we mentioned previous sections. 

## How to Use RECUR

In this section, we will dive deep into the options you can have to run RECUR. The commands shown in this section will assume that you have RECUR installed on your machine. 

### Options Overview

```bash
  -f <dir/file>                Protein or codon alignment in FASTA format [Required]
  -s <str>                     <AA|CODON> [Required][Default: CODON1]
  --outgroups <dir/file/str>   List of outgroup species [Required]
  --num-alignments <int>       Number of simulated alignments for p-value estimation [Default: 1000]
  -te <dir/file>               Complete constraint tree [Default: Estimated from alignment]
  -m <str>                     Model of sequence evolution [Default: estimated from alignment]
  -nt <int>                    Number of threads provided to IQ-TREE2
  -t <int>                     Number of threads used for RECUR internal processing 
  --seed <int>                 Random starting see number [Default: 8]
  -o <txt>                     Results directory [Default: same directory as MSA files]
  -iv <str>                    IQ-TREE2 path [Default: local]
```
Please note that the default values for `-t`, `-nt` are processor dependent. If you are following the installation step mentioned in the previous section, you can run one of the following commands to find out the actual default setting for your machine.
```bash
recur 
python3 recur.py
```
### Simple Usage

The minimal requirements of RECUR is a MSA (protein or codon) in FASTA format with the sequence type specified and a defined outgroup species or clade. e.g., 

>`recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON>`

* `--outgroups`: informs RECUR how to correctly root the tree. 

   You can either provide a .txt file with each outgroup species listed on a new line, or if you have a small number of outgroup species you can write the species on the command line. e.g., `--outgroups "SpeciesA,SpeciesB,SpeciesC"`. 

* `-st`: signals the sequence type in the MSA. 

   For a protein MSA `-st AA` should be provided. For a codon MSA, the different NCBI genetic codes can be specified (defined [here](http://www.iqtree.org/doc/Substitution-Models#codon-models)). `-st CODON` will use the standard genetic code. 

### Advanced Usage

- Using a constraint tree 

To specify the topology of the phylogeny used by RECUR the user can provide a constraint tree using the -te flag. The argument is a file containing a tree in Newick format. E.g.,

```bash
recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON> -te <treefile>
```

- Providing a model of evolution

A model of sequence evolution (as long as it is supported by IQ-TREE2) can be provided using the `-m` flag. E.g.,

```bash
recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON> -te <treefile> -m <model_of_evolution>
```

- Running RECUR on a directory 

To free you from typing multiple commands in a terminal to run on multiple genes, RECUR provides an option to run on a folder. E.g., 

```bash
recur [options] -f <directory> --outgroups <directory> -st <AA|CODON> -te <directory>
```
For example, if you have three genes files, each have a different set of outgroups and tree files. You can place those outgroups files and tree files with the corresponding gene names inside the same MSA data folder. E.g.,
<p align="center">
  <img src="./docs/images/RECUR_input_structure_1.PNG" alt="RECUR input structure 1" width="250"/>
</p>

If your genes share the same outgroups and tree, you only need to create a single outgroups file and tree file in your data folder. Those two files will be share by all the genes for your RECUR analysis. E.g., 

<p align="center">
  <img src="./docs/images/RECUR_input_structure_2.PNG" alt="RECUR input structure 2" width="250"/>
</p>

> **Important Information when running RECUR on a directory**:
>  * `<alignment_file>`: can only have `.aln`, `fasta`, or `fa` as the file extensions. (When running RECUR directly on an alignment file, extension requirement does not apply.)
>  * `<outgroups_file>`: needs to have `.outgroup` in the file name.
>  * `<treefile>`: needs to have `.tree` in the file name.

- Running RECUR with parallel processing

To speed up the analysis, RECUR has introduced both multiprocessing and multithreading in the package. By default, the number of threads will be automatically calculated based the configuration of your machine. Alternatively, they can be user specified using the following three options. 

- `-t`: Number of threads used for the RECUR algorithms
- `-nt`: Number of threads for IQ-TREE2 to run in parallel

Note, the number of threads cannot exceed the computer allowance. 

> * `-t` <= maximum_num_logical_threads
> * `-nt` <= maximum_num_logical_threads

For instance, if your computer has 8 logical threads, you can have the following setup,

```bash
recur -f example_alignments.aln -st AA --outgroups example_alignments.outgroups.txt -t 8 -nt 8
```

### Output Structure

If you do not specify the output folder using `-o`, the results are located in the same directory as the MSA files.

<p align="center">
  <img src="./docs/images/RECUR_output_structure.PNG" alt="RECUR output structure" width="500"/>
</p>

The list of recurrent substitutions is found in the `*.recur.tsv` file. Additional results of intermediate steps can be found in the `.recur` folder. 

Inside the `*.recur` folder, you will find two `*.txt` files, i.e., `Citation.txt` and `Log.txt`, and four subdirectories, i.e., `Real_Phylogeny`, `Infered_Sequences`, `Monte_Carlo_Simulation` and `Substitution_Matirces`. 

- `Real_Phylogeny`: contains the tree constructed from the alignment (`*.treefile`), the associated .iqtree file (`*.iqtree`) with information on model testing and the raw ancestral state reconstruction file (`*.state`) output from IQ-TREE2.
- `Inferred_Sequences`: contains the inferred ancestral sequences with node labels corresponding to the tree in the Real_Phylogeny directory.
- `Substitution_Matirces`: contains a `*.substituion_matrix.tsv` file and a `*.accum_substituion_matrix.txt` file. The `*.substituion_matrix.tsv` summarises the substitutions that occurred at each site in the protein across the phylogeny, while the `*.accum_substituion_matrix.txt` presents the accumulated counts for each amino acid substitution type across all residues.

For example, by running RECUR on the `example_alignments.aln` file, the `example_alignments.aln.substituion_matrix.tsv` would have the following content:
<p align="center">
  <img src="./docs/images/RECUR_substitution_matrices.PNG" alt="RECUR substitution matrices" width="300"/>
</p>

`Site` represents the position in the protein alignment. The `Parent>Child:MutCount` column contains a list of amino acid substitutions that occurred at that site. For example, at site 2 a lysine (K) mutated to a serine (S) twice. A dash signifies no substitution was inferred at that site.

As for the `RowIndex` and `ColIndex` columns, they index the substitution in a mutation matrix.  The arrangement of the mutation matrix is the same as that in the accumulated subsititon matrix. Here, the accumulated subsititon matrix is obtained by summing up all the mutation matrix at each site.

<p align="center">
  <img src="./docs/images/RECUR_accum_substitution_matrices.PNG" alt="RECUR accum substitution matrices" width="700"/>
</p>

- `Monte_Carlo_Simulation`: contains the simulated alignments based on the best model of evolution,  phylogeny and root sequence of the subtree excluding the outgroup. The number of alignment files in this directory will match the number of Monte Carlo Simulation used in your analysis. It can be adjusted by setting `--num-alignments` to a different value.

### Interpretation of the Recurrence List

<p align="center">
  <img src="./docs/images/RECUR_recurrence_list.PNG" alt="RECUR recurrence list" width="500"/>
</p>

The `*.recur.tsv` file contains a list of recurrent substitutions, i.e., an amino acid substitution that occurred more than once at a specific site, identified in the alignment. Here is the breakdown of each columns inside this output file:
* `Site`: indicates the position in the protein MSA.
* `Parent`: indicates the ancestral amino acid residue.
* `Child`: indicates the descent amino acid residue.
* `Recurrence`: indicates the number of times that amino acid substitution occurred at that site across the phylogeny. 
* `Reversion`: indicates the number of times the reverse amino acid substitution occurred at that site across the phylogeny. 
* `P-value`: signifies the probability that the substitution has not occurred by chance (for details of its calculation please read Robbins et al., 2024).
* `AllSiteSubs`: provides a list of all other substitutions that occurred at that site. 
* `SiteComposition`: indicates the count of each residue at that site in the protein MSA. 

### A Note on Reproducibility

According to Shen et al. (2020), irreproducibility in maximum likelihood phylogenetic inference is a significant issue [^1]. Regardless of the method related parameters that would affect the reproducibility of the final outputs, different random starting seed number, number of threads and processor type can also introduce uncertainties to the results. Such effect can be observed by setting different `-nt`, or setting different seed number using `--seed` in RECUR for each run. When running RECUR on different machine with different processor type even with the fixed `-nt` and `--seed`, the results can still be different. Nevertheless, the results should stay the same for each run when running RECUR on the same machine with the fixed input parameters.  

## Citations

### Credits and Acknowledgements

## Contributing

If you find a bug :bug:, please open a [bug report](https://github.com/).
If you have an idea for an improvement or new feature, please open a [feature request]().

## References
[^1]: *Shen, XX., Li, Y., Hittinger, C.T. et al.* **An investigation of irreproducibility in maximum likelihood phylogenetic inference.** Nat Commun 11, 6096 (2020). [![DOI:10.1038/s41467-020-20005-6](https://raw.githubusercontent.com/OrthoFinder/RECUR/main/docs/images/doi-badge.svg)](https://doi.org/10.1038/s41467-020-20005-6)

