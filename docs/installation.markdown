---
layout: page
title: Installation
permalink: /installation/
---

## Getting Started with RECUR

  - [Run RECUR in conda](#run-recur-in-conda)
    - [Install RECUR via bioconda](#install-recur-via-bioconda)
    - [Run RECUR in conda via Scoure Code](#run-recur-in-conda-via-scoure-code)
  - [Run RECUR on Linux](#run-recur-on-linux-macos)
    - [Install RECUR via make/gmake](#install-recur-via-make/gmake)
    - [Conventional Installation](#conventional-installation)
  - [Run RECUR in a Docker Container](#run-recur-in-a-docker-container)

  The recurrence analysis implemented by RECUR utilises IQ-TREE (v2.4.0) phylogenomic software package to infer a phylogeny, the ancestral node sequences and to generate simulated alignments. Although IQ-TREE comes with three versions, i.e., Linux, MacOS and Windows, the prefered OS for RECUR is Linux. RECUR can be run with or without installation. In either case, it is recommended that it is run inside a virtual environment. If you wish to run RECUR on a Windows or a MacOS machine, you can either install [IQ-TREE](http://www.iqtree.org/#download) before launching RECUR, or download the source package of RECUR (bundled with a matching IQ-TREE2 build) from the [RECUR release page](https://github.com/OrthoFinder/RECUR/releases/)  or [Download page](https://orthofinder.github.io/RECUR/download/). 

  >**⚠️ Note: RECUR requires IQ-TREE version 2.0 or higher to function correctly. Ensure that your environment meets this requirement before proceeding.**


### Run RECUR in conda

  The simplest way to use RECUR is to install it via conda. RECUR is available on the [bioconda channel](https://anaconda.org/bioconda/recur).

  I you have miniconda installed, you can simply run the following commands to create a `recur_env` namespace, install RECUR in that namespace, then print out the version and test it on the `ExampleData` dataset.

  ```bash
  conda create -n recur_env python=3.12
  conda activate recur_env
  conda install bioconda::recur 
  recur --version
  recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
  ```
  Once you've finished your analysis and want to exit the current environment, simply run:
  ```bash
  conda deactivate
  ```
  If you'd like to completely remove the environment where RECUR was installed, you can do so with:
  ```bash
  conda deactivate
  conda remove -n recur_env --all
  ```
  The first command deactivates the environment, and the second one deletes it entirely. If you want to use RECUR again in the future, you’ll need to recreate the environment and reinstall the package.

### Run RECUR on Linux 
#### Conventional Installation

  **Download the source code directly**
  Pre-compiled RECUR v1.0.0 archives are available for all major operating systems and CPU architectures. For the full list of files and step-by-step installation instructions, visit the [RECUR v1.0.0](https://orthofinder.github.io/RECUR/release/31/03/2025/RECUR-v1.0.0.html)

  **Install directly from GitHub**

  ```bash
  python3 -m venv recur_env 
  . recur_env/bin/activate
  pip install git+https://github.com/OrthoFinder/RECUR.git
  recur --version
  recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
  ```

  To remove RECUR, you can simply deactivate the virtual environment and remove it. 
  ```bash
  deactivate
  rm -rf recur_env
  ```

  **Clone the repository and install locally**

  ```bash
  git clone https://github.com/OrthoFinder/RECUR.git
  cd RECUR 
  python3 -m venv recur_env 
  . recur_env/bin/activate
  pip install .
  recur --version
  recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
  ```
 **Clone the repository and without installation**

  Again, if you do not wish to install RECUR please run 
  ```bash
  git clone https://github.com/OrthoFinder/RECUR.git
  cd RECUR 
  python3 -m venv recur_env 
  . recur_env/bin/activate
  pip install -r requirements.txt 
  python3 recur.py --version # running RECUR without installing it
  python3 recur.py -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
  ```
  
  >Please note that the GitHub release of RECUR includes IQ-TREE v2.4.0 compiled for 64-bit Linux. If your system uses a different OS or CPU architecture, that binary won’t run and RECUR will fail to start. In that case, either install the appropriate IQ-TREE build for your platform from the official [IQ-TREE download page](http://www.iqtree.org/#download) before launching RECUR, or download the source package of RECUR (bundled with a matching IQ-TREE build) from the [RECUR release page](https://github.com/OrthoFinder/RECUR/releases/) or [Download page](https://orthofinder.github.io/RECUR/download/).

  To deactivate the virtual environment, please run:
  ```bash
  deactivate
  ```

  If you followed the above instructions and cloned the repo locally, you can run the following command to remove RECUR 
  ```bash
  deactivate
  cd ..
  rm -rf RECUR
  ```
  This will delete both the environment and all files related to RECUR. If you plan to use RECUR again in the future, you’ll need to recreate the environment and reinstall the package.

  > Note: If python3 doesn't work, try using python instead, or check your `/usr/bin` directory to determine which version of Python is installed on your system. Please note that the current version of RECUR requires Python version 3.9 or higher, but no greater than 3.13.

#### Installation via make/gmake

  A Makefile is provided to simplify the installation process. Either `make` or `gmake` can be used to run it. Run the following commands to grab a copy of the RECUR scoure code from GitHub, and install RECUR inside a Python virtual enviroment. Again, if `git` is unavailable, please use other methods such as `wget` or `curl` to download the source code.

  ```bash
  git clone https://github.com/OrthoFinder/RECUR.git
  cd RECUR 
  make install USE_CONDA=false
  ```
  When you run `make install USE_CONDA=false`, a prompt message will pop up on your terminal and ask your permission to use or create a `~/local/bin` under your home directory. If `y`, IQ-TREE and RECUR binaries will be saved in that directory on completion of the installation process. If `n`, you will need to provide a name of the directory you wish to create or use relative to your home directory.

  Please note that even if RECUR has a IQ-TREE binary shipped with it, running `make install USE_CONDA=false` will reinstall IQ-TREE to `~/local/bin` if you do not have IQ-TREE installed globally. If you already have IQ-TREE installed globally, the isntallation of IQ-TREE will be skipped. The same as RECUR. By default, `make install USE_CONDA=false` will install the latest version of IQ-TREE and RECUR. If you wish to override your older version of those two softwares, you can run `make install USE_CONDA=false FORCE=true`. 

  You can also utilise the Makefile to only update the IQ-TREE binary. Running `make install_iqtree2 USE_CONDA=false FORCE=true` will only force to install/update the IQ-TREE binary inside the `~/local/bin`. If you have `sudo` access, you can also run `make install_iqtree2 USE_CONDA=false FORCE=true SYSTEM_WIDE=true` to install IQ-TREE inside `/usr/local/bin`. The default IQ-TREE version is `2.4.0`.

  Once the IQ-TREE binary has installed, `make install USE_CONDA=false` will create a virtual enviroment named `recur_env` inside the `RECUR` scoure code directory, then install RECUR inside that virtual environment. Having installed RECUR, it will copy the recur binary file from the virtual environment to a `~/local/bin`.  

  After running the above commands, please start a new terminal session or close your terminal and restart a new session. Then you can run RECUR globaly without needing to activate the virtual enviroment where you installed it everytime. 

  Verify that the installation is successful by running 

  ```bash
  recur --version
  recur -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
  ```
  To remove the installed RECUR package, it's environment and Python dependencies, please run 

  ```bash
  make clean_recur
  ```

  However, this command does not remove IQ-TREE in your path. You need to run a separate command to remove the binary, namely
  ```bash
  make clean_iqtree2
  ``` 

  > You can also run `make conda_install` to install RECUR inside the `recur_env` namespace inside conda (**NOT AVAILABEL at the moment**). 
  > To remove the `recur_env` namespace, you can simply run `make clean_conda_env`

### Run RECUR in a Docker Container

  Apart from the above methods, you can also run RECUR inside a container. You can find the RECUR image on [Docker Hub](https://hub.docker.com/r/orthofinder/recur).

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

    If you need to install the Docker Engine before you can run the RECUR container, please find the right docker engine to install on your server in [here](https://docs/images.docker.com/engine/install/).

  With either Docker Desktop or Docker Engine installed on your machine, you can simply run the following command to test if you can run the RECUR container after you have logged in. 
  ```
  docker container run -it --rm orthofinder/recur:v1.0.0
  ```
  To run the RECUR container on your dataset, you will need to create a folder which contains your data in your current working directory. For instance, you have a data folder called MyData which contains a protein alignment file called `my_alignment.aln` and a file called `my_alignment.outgroups.txt` that contains all the outgroups. Optionally you can also provide a tree file, e.g., `MyData/my_alignments.tree.txt`. With all the required inputs, you can run the following command to start the RECUR container and make it run your dataset.

  ```bash
  docker run -it --rm \
      -v $(pwd)/MyData:/usr/src/recur/MyData \
      -e LOCAL_UID=$(id -u) -e LOCAL_GID=$(id -g) \
      orthofinder/recur:v1.0.0 \
      -f MyData/my_alignments.aln \
      -st AA \
      --outgroups MyData/my_alignments.outgroups.txt \
      -te MyData/my_alignments.tree.txt
  ```
  Please note that arguments behind `orthofinder/recur:v1.0.0` will be the same as you run RECUR directly as we mentioned previous sections.

  Here the `--rm` flag tells Docker to automatically remove the container once it exits, ensuring no leftover container remains after the run. If you wish you remove the RECUR image, please run 

  ```bash
  docker rmi orthofinder/recur:v1.0.0
  ```