.ONESHELL:

.DEFAULT_GOAL := install

ifndef MAKE_VERSION
$(error This Makefile requires GNU Make. Please use gmake instead of make.)
endif

QUIET ?= false
QUIET := $(shell echo $(QUIET) | tr '[:upper:]' '[:lower:]')
FORCE ?= false
FORCE := $(shell echo $(FORCE) | tr '[:upper:]' '[:lower:]')
CONDA_PYTHON_VERSION ?= 3.12
PYTHON_VERSION ?= python3
RECUR_ENV_DEFAULT := recur_env
ENV_NAME ?= $(RECUR_ENV_DEFAULT)

SYSTEM_WIDE ?= false
SYSTEM_WIDE := $(shell echo $(SYSTEM_WIDE) | tr '[:upper:]' '[:lower:]')
HOME_DIR := $(if $(HOME),$(HOME),$(shell echo ~))

USE_CONDA ?= true
USE_CONDA := $(shell echo $(USE_CONDA) | tr '[:upper:]' '[:lower:]')

ifeq ($(USE_CONDA),false)
	PROMPT_USER_INSTALL_DIR = \
		@if [ -d "$(HOME_DIR)/local/bin" ]; then \
			read -p "Directory $(HOME_DIR)/local/bin exists. Do you want to use it? (y/n) " choice; \
			case $$choice in \
				[yY]*) echo "$(HOME_DIR)/local/bin";; \
				[nN]*) \
					read -p "Enter a new directory name (relative to $(HOME_DIR)): " new_dir; \
					USER_INSTALL_DIR=$(HOME_DIR)/$$new_dir; \
					echo $$USER_INSTALL_DIR; \
					mkdir -p $$USER_INSTALL_DIR || { echo "Error creating directory $$USER_INSTALL_DIR. Exiting."; exit 1; }; \
					;; \
				*) echo "Invalid choice. Exiting."; exit 1; \
			esac; \
		else \
			echo "$(HOME_DIR)/local/bin"; \
		fi
endif

USER_INSTALL_DIR := $(shell $(PROMPT_USER_INSTALL_DIR))
SYSTEM_INSTALL_DIR := /usr/local/bin
RECUR_DIR := $(USER_INSTALL_DIR)

ifeq ($(SYSTEM_WIDE),true)
    BINARY_INSTALL_DIR := $(SYSTEM_INSTALL_DIR)
    SUDO_PREFIX := $(if $(shell [ "$(USER)" = "root" ] || echo 1),sudo)
else
    BINARY_INSTALL_DIR := $(USER_INSTALL_DIR)
    SUDO_PREFIX :=
endif

IQTREE_DEFAULT_VERSION := 3.1.1
IQTREE_VERSION ?= $(IQTREE_DEFAULT_VERSION)
IQTREE3_VERSION := 3.1.1

# URLs for IQ-TREE3 urlS
IQTREE3_LINUX_INTEL := https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-Linux-intel.tar.gz
IQTREE3_LINUX_ARM := https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-Linux-arm.tar.gz

IQTREE3_MACOS_UNIVERSAL := https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-macOS.zip
IQTREE3_MACOS_INTEL := https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-macOS-intel.zip
IQTREE3_MACOS_ARM := https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-macOS-arm.zip

IQTREE3_BINARY := $(BINARY_INSTALL_DIR)/iqtree3


PYTHON := ./$(ENV_NAME)/bin/python3
PIP := ./$(ENV_NAME)/bin/pip
VENV_BIN := ./$(ENV_NAME)/bin


check_conda:
	@if ! command -v conda > /dev/null; then \
		echo "Error: Conda is not installed. Please install Conda first."; \
		exit 1; \
	fi
	@echo "Conda is installed. Activating Conda..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh || { echo "Error: Failed to initialize Conda. Exiting."; exit 1; }; \
	echo "Conda activated successfully."

create_conda_env: check_conda
	@echo "Checking if Conda environment $(ENV_NAME) exists..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh && \
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]" && [ "$(FORCE)" = "0" ]; then \
		echo "Conda environment $(ENV_NAME) already exists. Skipping creation."; \
	else \
		if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then \
			echo "Forcing recreation of Conda environment $(ENV_NAME)..."; \
			conda env remove -n $(ENV_NAME) -y || { echo "Error: Failed to remove existing Conda environment. Exiting."; exit 1; }; \
		fi; \
		echo "Creating Conda environment: $(ENV_NAME) with Python $(CONDA_PYTHON_VERSION)..."; \
		conda create -n $(ENV_NAME) python=$(CONDA_PYTHON_VERSION) -y || { echo "Error: Failed to create Conda environment. Exiting."; exit 1; }; \
		echo "Conda environment $(ENV_NAME) with Python $(CONDA_PYTHON_VERSION) created successfully."; \
	fi


conda_install_iqtree3: create_conda_env
	@echo "Checking global paths for IQ-TREE3..."; \
	iqtree3_exists=$$(command -v iqtree3 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree3_exists" = "0" ]; then \
		echo "Installing IQ-TREE3 version $(IQTREE3_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::iqtree=$(IQTREE3_VERSION) -y || { echo "Error: Failed to install IQ-TREE3. Exiting."; exit 1; }; \
		echo "IQ-TREE3 version $(IQTREE3_VERSION) installed successfully."
	elif [ "$$iqtree3_exists" = "1" ]; then \
		echo "IQ-TREE3 already exist globally. Skipping installation."; \
	fi


conda_install_recur: create_conda_env
	@echo "Checking global paths for RECUR..."; \
	recur_exists=$$(command -v recur > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$recur_exists" = "0" ]; then \
		echo "Installing  RECUR version $(RECUR_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::recur=$(RECUR_VERSION) -y || { echo "Error: Failed to install RECUR. Exiting."; exit 1; }; \
		echo " RECUR version $(RECUR_VERSION) installed successfully."
	elif [ "$$recur_exists" = "1" ]; then \
		echo " RECUR already exist globally. Skipping installation."; \
	fi

conda_install: conda_install_recur conda_install_iqtree3
	@echo "You have now installed RECUR $(RECUR_VERSION) and its dependencies in $(ENV_NAME)!"

clean_conda_env:
	@echo "Checking if Conda environment $(ENV_NAME) exists..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh && \
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then \
		echo "Removing $(ENV_NAME) from Conda..."; \
		conda remove -n $(ENV_NAME) --all -y || { echo "Error: Failed to remove $(ENV_NAME) from Conda. Exiting."; exit 1; }; \
		echo "You have successfully removed $(ENV_NAME) from Conda!"; \
	else \
		echo "Conda environment $(ENV_NAME) does not exist. Skipping removal."; \
	fi

# make_usr_bin:
# 	@if [ ! -d "$(USER_INSTALL_DIR)" ]; then \
# 		echo "Directory $(USER_INSTALL_DIR) does not exist. Creating it..."; \
# 		mkdir -p $(USER_INSTALL_DIR); \
# 	fi; \
# 	echo "Checking if $(USER_INSTALL_DIR) is already in the PATH..."; \
# 	if ! grep -qx 'export PATH="$(USER_INSTALL_DIR):$$PATH"' ~/.bashrc; then \
# 		echo "Adding $(USER_INSTALL_DIR) to the PATH permanently."; \
# 		echo 'export PATH="$(USER_INSTALL_DIR):$$PATH"' >> ~/.bashrc || { echo "Error: Failed to update PATH in ~/.bashrc. Exiting."; exit 1; }; \
# 		echo "PATH update added to ~/.bashrc. Please restart your shell or run 'source ~/.bashrc' to apply changes."; \
# 	else \
# 		echo "$(USER_INSTALL_DIR) is already in the PATH. Skipping addition to ~/.bashrc."; \
# 	fi


make_usr_bin:
	@if [ ! -d "$(USER_INSTALL_DIR)" ]; then \
		echo "Directory $(USER_INSTALL_DIR) does not exist. Creating it..."; \
		mkdir -p $(USER_INSTALL_DIR); \
	fi; \
	SHELL_NAME=$$(basename "$${SHELL}"); \
	if [ "$${SHELL_NAME}" = "zsh" ]; then \
		SHELL_RC="$${HOME}/.zshrc"; \
	elif [ "$${SHELL_NAME}" = "bash" ]; then \
		if [ "$$(uname)" = "Darwin" ]; then \
			SHELL_RC="$${HOME}/.bash_profile"; \
		else \
			SHELL_RC="$${HOME}/.bashrc"; \
		fi; \
	else \
		SHELL_RC="$${HOME}/.profile"; \
	fi; \
	echo "Using shell configuration file: $${SHELL_RC}"; \
	echo "Checking if $(USER_INSTALL_DIR) is already in the PATH..."; \
	if ! grep -qx 'export PATH="$(USER_INSTALL_DIR):$$PATH"' $${SHELL_RC}; then \
		echo "Adding $(USER_INSTALL_DIR) to the PATH permanently."; \
		echo 'export PATH="$(USER_INSTALL_DIR):$$PATH"' >> $${SHELL_RC} || { echo "Error: Failed to update PATH in $${SHELL_RC}. Exiting."; exit 1; }; \
		echo "PATH update added to $${SHELL_RC}. Please restart your shell or run 'source $${SHELL_RC}' to apply changes."; \
	else \
		echo "$(USER_INSTALL_DIR) is already in the PATH. Skipping addition to $${SHELL_RC}."; \
	fi


install_iqtree3: make_usr_bin
	@echo "Checking global paths for IQ-TREE3..."; \
	iqtree3_exists=$$(command -v iqtree3 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree3_exists" = "0" ]; then \
		echo "Detecting system architecture..."; \
		OS=$$(uname -s); ARCH=$$(uname -m); \
		if [ "$$OS" = "Linux" ]; then \
			if [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE3_LINUX_INTEL); \
			elif [ "$$ARCH" = "aarch64" ]; then \
				IQTREE_URL=$(IQTREE3_LINUX_ARM); \
				echo "Downloading IQ-TREE3 IQTREE3_LINUX_ARM version..."; \
			fi; \
		elif [ "$$OS" = "Darwin" ]; then \
			if [ "$$ARCH" = "arm64" ]; then \
				IQTREE_URL=$(IQTREE3_MACOS_ARM); \
			elif [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE3_MACOS_INTEL); \
			else \
				IQTREE_URL=$(IQTREE3_MACOS_UNIVERSAL); \
			fi; \
		else \
			echo "Error: Unsupported operating system: $$OS"; exit 1; \
		fi; \
		echo "Downloading IQ-TREE3 from $$IQTREE_URL..."; \
		temp_dir=$$(mktemp -d); \
		download_path=$$temp_dir/iqtree3-src; \
		if [ "$(QUIET)" = "true" ]; then \
			wget -O $$download_path $$IQTREE_URL > /dev/null 2>&1 || { echo "Error: Failed to download IQ-TREE3."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget -O $$download_path $$IQTREE_URL || { echo "Error: Failed to download IQ-TREE3."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Extracting IQ-TREE3..."; \
		if echo "$$IQTREE_URL" | grep -q '.tar.gz'; then \
			if [ "$(QUIET)" = "true" ]; then \
				tar -xzf $$download_path -C $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE3 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				tar -xzf $$download_path -C $$temp_dir || { echo "Error: Failed to extract IQ-TREE3 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		elif echo "$$IQTREE_URL" | grep -q '.zip'; then \
			if [ "$(QUIET)" = "true" ]; then \
				unzip -o $$download_path -d $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE3 zip file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				unzip -o $$download_path -d $$temp_dir || { echo "Error: Failed to extract IQ-TREE3 zip file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		else \
			echo "Error: Unknown file format for IQ-TREE3."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Locating extracted IQ-TREE3 binary..."; \
		iqtree3_binary=$$(find $$temp_dir -type f -name "iqtree*" -executable | head -1); \
		if [ -z "$$iqtree3_binary" ]; then \
			echo "Error: IQ-TREE3 binary not found after extraction."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Moving IQ-TREE3 binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv $$iqtree3_binary $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move IQ-TREE3 binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "IQ-TREE3 installation completed successfully."; \
	else \
		iqtree3_path=$$(command -v iqtree3); \
		echo "IQ-TREE3 already exists at: $$iqtree3_path. Skipping installation."; \
	fi


clean_iqtree3:
	@echo "Cleaning user-specific IQ-TREE3 installation..."; \
	$(SUDO_PREFIX) rm -f "$(IQTREE3_BINARY)" && \
	echo "User-specific IQ-TREE3 successfully removed." || \
	{ echo "Error: Failed to remove user-specific IQ-TREE3 binary from $(IQTREE3_BINARY). Exiting."; exit 1; }; \


venv:
	@echo "Checking for existing virtual environment $(ENV_NAME)..."
	@if [ -d "$(ENV_NAME)" ] && [ "$(FORCE)" = "0" ]; then \
		echo "Virtual environment $(ENV_NAME) already exists."; \
		echo "Activating $(ENV_NAME)..."; \
		. $(ENV_NAME)/bin/activate; \
		echo "Virtual environment $(ENV_NAME) activated successfully."; \
	else \
		echo "Creating virtual environment $(ENV_NAME) using $(PYTHON_VERSION)..."; \
		if ! command -v $(PYTHON_VERSION) > /dev/null; then \
			echo "Error: $(PYTHON_VERSION) is not installed or not in PATH. Exiting."; \
			exit 1; \
		fi; \
		rm -rf $(ENV_NAME); \
		$(PYTHON_VERSION) -m venv $(ENV_NAME) || { echo "Error: Failed to create virtual environment. Exiting."; exit 1; }; \
		chmod +x $(ENV_NAME)/bin/activate; \
		. $(ENV_NAME)/bin/activate; \
		echo "Virtual environment $(ENV_NAME) created and activated successfully."; \
	fi

install_dependencies: venv 
	$(PIP) install -r requirements.txt

install: install_iqtree3 venv make_usr_bin 
	@echo "Checking global paths for RECUR..."; \
	recur_exists=$$(command -v recur > /dev/null 2>&1 && echo 1 || echo 0); \
	if [ "$(FORCE)" = "true" ] || [ "$$recur_exists" = "0" ]; then \	
		echo "Installing RECUR..."; \
		if [ "$(QUIET)" = "true" ]; then \
			$(PIP) install -e . > /dev/null 2>&1 || { echo "Error: Failed to install RECUR. Exiting."; exit 1; }; \
		else \
			$(PIP) install -e . || { echo "Error: Failed to install RECUR. Exiting."; exit 1; }; \
		fi; \
		mkdir -p $(USER_INSTALL_DIR) || { echo "Error: Failed to create directory $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
		echo "Copying RECUR to $(USER_INSTALL_DIR)..."; \
		cp $(VENV_BIN)/recur $(USER_INSTALL_DIR)/ && \
		echo "RECUR successfully copied to $(USER_INSTALL_DIR)." || \
		{ echo "Error: Failed to copy RECUR to $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
	elif [ "$$recur_exists" = "1" ]; then \
		recur_path=$$(command -v recur); \
		echo "RECUR already exists at: $$recur_path. Skipping installation."; \
	fi


run: install_iqtree3 install
	@echo "Running RECUR..."
	@if [ -f "$(VENV_BIN)/recur" ]; then \
		$(VENV_BIN)/recur -f ExampleData -st AA --outgroups ExampleData
	elif [ -f "$(RECUR_DIR)/recur" ]; then \
		$(RECUR_DIR)/recur -f ExampleData -st AA --outgroups ExampleData
	else \
		echo "RECUR not found. Installing..."; \
		$(MAKE) install; \
	fi

clean_recur:
	@echo "Remove RECUR and it's environment and dependencies..."
	rm -rf $(ENV_NAME) 
	rm -rf **/__pycache__
	rm -rf ./src/recur/__pycache__
	rm -rf ./src/recur/utils/__pycache__
	rm -rf ./src/recur/run/__pycache__
	rm -rf ./build ./dist **/recur*.egg-info

	if [ -f "$(RECUR_DIR)/recur" ]; then \
		rm -f "$(RECUR_DIR)/recur"
		echo "$(RECUR_DIR)/recur including it's environment have been removed..."; \
	else \
		echo "RECUR not found..."; \
	fi


.PHONY: make_usr_bin clean clean_iqtree3 purge clean_conda_venv
