# RECUR
Finding Recurrent Substitutions in Multiple Sequence Alignments

## Introduction
![RECUR method workflow](./docs/RECUR_workflow_figure.png)
Figure 1: The RECUR workflow. 

The required input is either a protein or codon multiple sequence alignment (in FASTA format) and a defined outgroup species or clade. The output of RECUR is a list of recurrent amino acid substitutions that have been inferred to have occurred within the phylogeny (file suffix .recfinder.tsv). Outputs of intermediate steps, i.e. model selection, tree inference, ancestral state reconstruction and site substitution matrices, can be found in the .recfinder output directory. 

## Table of Contents

- [Getting started with RECUR](#getting-started-with-recur)
  - [Installing RECUR on Linux](#installing-recur-on-linux)
  - [Installing RECUR on Windows and MacOS](#installing-recur-on-windows-and-macos)
  - [Running RECUR in conda](#running-recur-in-conda)
  - [Running RECUR in a Docker Conatainer](#running-recur-in-a-docker-conatainer)
- [Usage](#usage)
- [Contributing](#contributing)

## Getting started with RECUR
The recurrence analysis implemented by RECUR utilises two main feaures from IQ-TREE2, i.e., [Ancestral sequence reconstruction](http://www.iqtree.org/doc/Command-Reference#ancestral-sequence-reconstruction) to infer the extinct node sequences and [Sequence simulators](http://www.iqtree.org/doc/AliSim) to build the simulated phylogeny. 

### Installing RECUR on Linux

  If you are working on a Linux machine or WSL2, runing RECUR is straigtforward. Installation of IQ-TREE2 is not necessary.

  Before installing any relevant dependencies or packages, it is recommanded that you create and activate a new virtual environment. You can do so by running:
  ```
  python3 -m venv .venv && . .venv/bin/activate
  ```
  *To deactivate the virtual enviroment* run
  ```
  deactivate
  ```
  1. Clone the repository
  ```
  git clone https://github.com/OrthoFinder/RECUR.git
  ```
  2. Navigate into the cloned repository
  ```
  cd RECUR
  ```
  3. Install the package [Optional]
  ```
  pip install .
  ```
  To test your installation, please run

  ```
  recur -f ExampleData -st AA --outgroups ExampleData
  ```
    
  If you do not wish to install the RECUR package, you can simply run the following command to install the required dependencies.

  ```
  pip install -r requirements.txt
  ```
  Then run 
  ```
  python3 recur.py -f ExampleData -st AA --outgroups ExampleData -te ExampleData
  ```
  to test the installation and the environment. 
  
  > NOTE: if python3 does not work, please try python.

### Installing RECUR on Windows and MacOS 

  If you are on a Windows machine or a MAC, you need to install the IQ-TREE2 first before you can run RECUR. You can download the latest version of IQ-TREE2 in [here](http://www.iqtree.org/#download).

  Once you have IQ-TREE2 installed, you can run the following commands based on your OS to create a virtual environment for RECUR to run.

  - **Windows**

  If you are on windows, please open a command prompt, and run the following commnad to create and activate the virtual environment.
  ```
  python -m venv win_venv
  win_venv\Scripts\activate.bat
  ```

  *To deactivate the virtual enviroment* run
  ```
  win_venv\Scripts\deactivate.bat
  ```
  - **macOS**

  If you are on macOS, please run the following commands in a terminal.

  ```
  python3 -m venv mac_venv
  source mac_venv/bin/activate
  ```

  *To deactivate the virtual enviroment* run
  ```
  deactivate
  ```
  Now that you have successfully created and activated a new virtual environment, you can follow the three steps mentioned in [Installing RECUR on Linux](#installing-recur-on-linux) to install RECUR or just the python dependencies based on your needs. Then run one of the following commands (depending one whether you have installed RECUR or not) to test if RECUR runs on your machine.

  * With RECUR installed 
    ```
    recur -f ExampleData -st AA --outgroups ExampleData -iv system
    ```
  * Without RECUR installed 
      ```
    python3 recur.py -f ExampleData -st AA --outgroups ExampleData -iv system
    ```
  where `-iv` stands for IQ-TREE2 version. By default, RECUR will use the Linux binary version from the package. 

### Running RECUR in conda

  Working in the conda environment can be the easiest when you do not have access to a Linux machine. You can create and activate a new environment by running:

    ```
    conda create -n recur_env python=3.12
    conda activate recur_env
    ```
  Then install IQ-TREE2 package by running one of the following:
  ```
    conda install bioconda::iqtree
    conda install bioconda/label/cf201901::iqtree
  ```
  Now that you have a suitable environment with IQ-TREE2 installed, you can follow the previous steps to install RECUR directly or the relevant dependencies.
  * With RECUR installed 
    ```
    recur -f ExampleData -st AA --outgroups ExampleData -te ExampleData -iv conda
    ```
  * Without RECUR installed 
      ```
    python3 recur.py -f ExampleData -st AA --outgroups ExampleData  -te ExampleData -iv conda
    ```

### Running RECUR in a Docker Conatainer

  Apart from the above methods, you can also run RECUR inside a container. 

  - Personal computer
  
    Before you can run the RECUR container, you need to have Docker Desktop installed. 
    - Windows: https://docs.docker.com/desktop/install/windows-install/
    - macOS: https://docs.docker.com/desktop/install/mac-install/
    - Linux: https://docs.docker.com/desktop/install/linux-install/

    Once you have the Docker Desktop installed, please luanch it and run the following command in the terminal to check if it is up and running.

    ```
    docker version
    ```
  - **Server**

    If you need to install the Docker Engine before you can run the RECUR container.

    Please find the right docker engine to install on your server in [here](https://docs.docker.com/engine/install/).

  With either Docker Decktop or Docker Engine installed on your machine, you can simply run the following command to test if you can run the RECUR container
  ```
  docker container run -it --rm orthofinder/recur:v1.0.0
  ```
  To run RECUR container on your dataset, you will need to create a folder which contains your data in your current working directory. For instance, you have a data folder called MyData, you can run the following command to start the RECUR container and make it run your dataset.
  ```
  docker container run -it --rm -v $(pwd)/MyData:/usr/src/RECUR/MyData orthofinder/recur:v1.0.0 recur -f MyData -st AA --outgroups MyData   
  ```
  Please note that arguments behind `orthofinder/recur:v1.0.0` will be the same as you run RECUR directly as we mentioned previous sections.

## Usage
(Instructions on how to use the project)

In this section, we will dive deep into the options you can have to run RECUR. The command shown in this section will be based on the assumption that you have RECUR installed on your machine. 

- **Simple usage**

  RECUR can run on either a protien or a CODON alignment. 
  ```
  recur [options] -f <dir/file> --outgroups <outgroup_species/dir/file> -st <AA|CODON>
  ```


- **Advanced usage**


## Contributing

If you find a bug :bug:, please open a [bug report](https://github.com/).
If you have an idea for an improvement or new feature, please open a [feature request]().


