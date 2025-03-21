---
layout: page
title: Installation
permalink: /installation/
---

## Getting Started with RECUR
The recurrence analysis implemented by RECUR utilises IQ-TREE2 phylogenomic software package to infer a phylogeny, the ancestral node sequences and to generate simulated alignments.

### Installing RECUR on Linux

  If you are working on a Linux machine or WSL2, running RECUR is straightforward. A separate installation of IQ-TREE2 is not necessary.

  Before installing any relevant dependencies or packages, it is recommended that you create and activate a new virtual environment. You can do so by running:
  ```bash
  python3 -m venv .venv && . .venv/bin/activate
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

  Once your RECUR job has done, you can run
  ```bash
  deactivate
  ```
  to deactivate the virtual environment. The next time when you use RECUR, you can simply activate the virtual environment by running
  ```bash
  . .venv/bin/activate
  ```
  There is no need to reinstall RECUR again, as it has been installed in your virtual environment if you followed the instruction.

  > Note: If python3 doesn't work, try using python instead, or check your `/usr/bin` directory to determine which version of Python is installed on your system. Please note that the current version of RECUR requires Python version 3.9 or higher, but no greater than 3.13.

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
  Now that you have successfully created and activated a new virtual environment, you can follow the three steps mentioned in [Installing RECUR on Linux](#installing-recur-on-linux) to install RECUR or just the python dependencies based on your needs. Then run one of the following commands (depending on whether you have installed RECUR or not) to test if RECUR runs on your machine.

  * With RECUR installed
    ```bash
    recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
    ```
  * Without RECUR installed
    ```bash
    python3 recur.py -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
    ```

  *To deactivate the virtual enviroment* run
  ```bash
  deactivate
  ```
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
    
    <br>
    Once you have the Docker Desktop installed, please launch it and run the following command in the terminal to check if it is up and running.

    ```bash
    docker version
    ```
  - **Server**

    If you need to install the Docker Engine before you can run the RECUR container, please find the right docker engine to install on your server in [here](https://docs/images.docker.com/engine/install/).

  With either Docker Desktop or Docker Engine installed on your machine, you can simply run the following command to test if you can run the RECUR container
  ```
  docker container run -it --rm orthofinder/recur:v1.0.0
  ```
  To run the RECUR container on your dataset, you will need to create a folder which contains your data in your current working directory. For instance, you have a data folder called MyData which contains a protein alignment file called `my_alignment.aln` and a file called `my_alignment.outgroups.txt` that contains all the outgroups, you can run the following command to start the RECUR container and make it run your dataset.
  ```
  docker container run -it --rm -v $(pwd)/MyData:/usr/src/RECUR/MyData:Z orthofinder/recur:v1.0.0 recur -f MyData/my_alignment.aln -st AA --outgroups MyData/my_alignment.outgroups.txt
  ```
  Please note that arguments behind `orthofinder/recur:v1.0.0` will be the same as you run RECUR directly as we mentioned previous sections.
