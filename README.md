# RECUR
Finding Recurrent Substitutions in Multiple Sequence Alignments

## Introduction


## Table of Contents

- [Features](#features)
- [Usage](#usage)
- [Contributing](#contributing)

## Features (Lizzie)
* Finding recurrent substitutions in multiple sequence alignments  

## Installation
The recurrence analysis implemented by RECUR utilises two main feaures from IQ-TREE2, i.e., [Ancestral sequence reconstruction](http://www.iqtree.org/doc/Command-Reference#ancestral-sequence-reconstruction) to infer the extinct node sequences and [Sequence simulators](http://www.iqtree.org/doc/AliSim) to build the simulated phylogeny. 

- Linux

    If you are working on a Linux machine or WSL2, runing RECUR is straigtforward. Installation of IQ-TREE2 is not necessary.

    Before installing any relevant dependencies or packages, it is recommanded that you create and activate a new virtual environment. You can do so by running:
    ```
    python3 -m venv .venv && . .venv/bin/activate
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
    python3 recur.py -f ExampleData -st AA --outgroups ExampleData
    ```
    to test the installation and the environment. 
    
    > NOTE: if python3 does not work, please try python.

- Windows and macOS
  If you are on a Windows machine or MAC, you need to install the IQ-TREE2 first before you can run RECUR. You can download the latest version of IQ-TREE2 in [here](http://www.iqtree.org/#download).

  Once you have IQ-TREE2 installed, you can follow the previous steps to install RECUR or just the dependencies. Then run oen of the following commands (depending one whether you have installed RECUR or not) to test if RECUR runs on your machine
  * With RECUR installed 
    ```
    recur -f ExampleData -st AA --outgroups ExampleData -iv system
    ```
  * Without RECUR installed 
      ```
    python3 recur.py -f ExampleData -st AA --outgroups ExampleData -iv system
    ```
  where `-iv` stands for IQ-TREE2 version. By default, RECUR will use the Linux binary version from the package. 

- Conda 
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
    recur -f ExampleData -st AA --outgroups ExampleData -iv conda
    ```
  * Without RECUR installed 
      ```
    python3 recur.py -f ExampleData -st AA --outgroups ExampleData -iv conda
    ```

## Usage
(Instructions on how to use the project)


## Contributing

If you find a bug :bug:, please open a [bug report](https://github.com/).
If you have an idea for an improvement or new feature, please open a [feature request]().

## Method Workflow
![RECUR method workflow](./docs/method_workflow.tif)