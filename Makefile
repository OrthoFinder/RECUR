.ONESHELL:

.DEFAULT_GOAL := install

ifndef MAKE_VERSION
$(error This Makefile requires GNU Make. Please use gmake instead of make.)
endif

QUIET ?= false
FORCE ?= false
CONDA_PYTHON_VERSION ?= 3.10
PYTHON_VERSION ?= python3
RECUR_ENV_DEFAULT := recur_venv
ENV_NAME ?= $(RECUR_ENV_DEFAULT)

SYSTEM_WIDE ?= 0
HOME_DIR := $(if $(HOME),$(HOME),$(shell echo ~))

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


USER_INSTALL_DIR := $(shell $(PROMPT_USER_INSTALL_DIR))
SYSTEM_INSTALL_DIR := /usr/local/bin
RECUR_DIR := $(USER_INSTALL_DIR)

ifeq ($(SYSTEM_WIDE),1)
    BINARY_INSTALL_DIR := $(SYSTEM_INSTALL_DIR)
    SUDO_PREFIX := $(if $(shell [ "$(USER)" = "root" ] || echo 1),sudo)
else
    BINARY_INSTALL_DIR := $(USER_INSTALL_DIR)
    SUDO_PREFIX :=
endif

IQTREE_DEFAULT_VERSION := 2.3.6
IQTREE_VERSION ?= $(IQTREE_DEFAULT_VERSION)

# URLs for IQ-TREE urlS
IQTREE_LINUX_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz
IQTREE_LINUX_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-arm.tar.gz

IQTREE_MACOS_UNIVERSAL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS.zip
IQTREE_MACOS_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-intel.zip
IQTREE_MACOS_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-arm.zip

IQTREE_BINARY := $(BINARY_INSTALL_DIR)/iqtree2

PYTHON := ./$(ENV_NAME)/bin/python3
PIP := ./$(ENV_NAME)/bin/pip
VENV_BIN := ./$(ENV_NAME)/bin/


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

conda_install_iqtree2: create_conda_env
	@echo "Checking global paths for IQ-TREE2..."; \
	iqtree2_exists=$$(command -v iqtree2 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree2_exists" = "0" ]; then \
		echo "Installing IQ-TREE2 version $(IQTREE_VERSION) in $(ENV_NAME)..."; \
		. $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda activate $(ENV_NAME) && \
		conda install bioconda::iqtree=$(IQTREE_VERSION) -y || { echo "Error: Failed to install IQ-TREE2. Exiting."; exit 1; }; \
		echo "IQ-TREE2 version $(IQTREE_VERSION) installed successfully."
	elif [ "$$iqtree2_exists" = "1" ]; then \
		echo "IQ-TREE2 already exist globally. Skipping installation."; \
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

conda_install: conda_install_recur conda_install_iqtree2
	@echo "You have now installed RECUR $(RECUR_VERSION) and its dependencies in $(ENV_NAME)!"

clean_conda_venv:
	@echo "Checking if Conda environment $(ENV_NAME) exists..."; \
	. $(shell conda info --base)/etc/profile.d/conda.sh && \
	if conda env list | grep -q "^$(ENV_NAME)[[:space:]]"; then \
		echo "Removing $(ENV_NAME) from Conda..."; \
		conda remove -n $(ENV_NAME) --all -y || { echo "Error: Failed to remove $(ENV_NAME) from Conda. Exiting."; exit 1; }; \
		echo "You have successfully removed $(ENV_NAME) from Conda!"; \
	else \
		echo "Conda environment $(ENV_NAME) does not exist. Skipping removal."; \
	fi

make_usr_bin:
	@if [ ! -d "$(USER_INSTALL_DIR)" ]; then \
		echo "Directory $(USER_INSTALL_DIR) does not exist. Creating it..."; \
		mkdir -p $(USER_INSTALL_DIR); \
	fi; \
	echo "Checking if $(USER_INSTALL_DIR) is already in the PATH..."; \
	if ! grep -qx 'export PATH="$(USER_INSTALL_DIR):$$PATH"' ~/.bashrc; then \
		echo "Adding $(USER_INSTALL_DIR) to the PATH permanently."; \
		echo 'export PATH="$(USER_INSTALL_DIR):$$PATH"' >> ~/.bashrc || { echo "Error: Failed to update PATH in ~/.bashrc. Exiting."; exit 1; }; \
		echo "PATH update added to ~/.bashrc. Please restart your shell or run 'source ~/.bashrc' to apply changes."; \
	else \
		echo "$(USER_INSTALL_DIR) is already in the PATH. Skipping addition to ~/.bashrc."; \
	fi

install_iqtree2: make_usr_bin
	@echo "Checking global paths for IQ-TREE2..."; \
	iqtree2_exists=$$(command -v iqtree2 > /dev/null && echo 1 || echo 0); \

	if [ "$(FORCE)" = "true" ] || [ "$$iqtree2_exists" = "0" ]; then \
		echo "Detecting system architecture..."; \
		OS=$$(uname -s); ARCH=$$(uname -m); \
		if [ "$$OS" = "Linux" ]; then \
			if [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE_LINUX_INTEL); \
			elif [ "$$ARCH" = "aarch64" ]; then \
				IQTREE_URL=$(IQTREE_LINUX_ARM); \
			else \
				echo "Error: Unsupported Linux architecture: $$ARCH"; exit 1; \
			fi; \
		elif [ "$$OS" = "Darwin" ]; then \
			if [ "$$ARCH" = "arm64" ]; then \
				IQTREE_URL=$(IQTREE_MACOS_ARM); \
			elif [ "$$ARCH" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE_MACOS_INTEL); \
			else \
				IQTREE_URL=$(IQTREE_MACOS_UNIVERSAL); \
			fi; \
		else \
			echo "Error: Unsupported operating system: $$OS"; exit 1; \
		fi; \
		echo "Downloading IQ-TREE2 from $$IQTREE_URL..."; \
		temp_dir=$$(mktemp -d); \
		download_path=$$temp_dir/iqtree2-src; \
		if [ "$(QUIET)" = "true" ]; then \
			wget -O $$download_path $$IQTREE_URL > /dev/null 2>&1 || { echo "Error: Failed to download IQ-TREE2."; rm -rf $$temp_dir; exit 1; }; \
		else \
			wget -O $$download_path $$IQTREE_URL || { echo "Error: Failed to download IQ-TREE2."; rm -rf $$temp_dir; exit 1; }; \
		fi; \
		echo "Extracting IQ-TREE..."; \
		if echo "$$IQTREE_URL" | grep -q '.tar.gz'; then \
			if [ "$(QUIET)" = "true" ]; then \
				tar -xzf $$download_path -C $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE2 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				tar -xzf $$download_path -C $$temp_dir || { echo "Error: Failed to extract IQ-TREE2 tar.gz file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		elif echo "$$IQTREE_URL" | grep -q '.zip'; then \
			if [ "$(QUIET)" = "true" ]; then \
				unzip -o $$download_path -d $$temp_dir > /dev/null 2>&1 || { echo "Error: Failed to extract IQ-TREE2 zip file."; rm -rf $$temp_dir; exit 1; }; \
			else \
				unzip -o $$download_path -d $$temp_dir || { echo "Error: Failed to extract IQ-TREE2 zip file."; rm -rf $$temp_dir; exit 1; }; \
			fi; \
		else \
			echo "Error: Unknown file format for IQ-TREE2."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Locating extracted IQ-TREE2 binary..."; \
		iqtree2_binary=$$(find $$temp_dir -type f -name "iqtree*" -executable | head -1); \
		if [ -z "$$iqtree2_binary" ]; then \
			echo "Error: IQ-TREE2 binary not found after extraction."; rm -rf $$temp_dir; exit 1; \
		fi; \
		echo "Moving IQ-TREE2 binary to $(BINARY_INSTALL_DIR)..."; \
		$(SUDO_PREFIX) mv $$iqtree2_binary $(BINARY_INSTALL_DIR) || { echo "Error: Failed to move IQ-TREE2 binary."; rm -rf $$temp_dir; exit 1; }; \
		rm -rf $$temp_dir; \
		echo "IQ-TREE2 installation completed successfully."; \
	else \
		iqtree2_path=$$(command -v iqtree2); \
		echo "IQ-TREE2 already exists at: $$iqtree2_path. Skipping installation."; \
	fi


clean_iqtree2:
	@echo "Cleaning user-specific IQ-TREE2 installation..."; \
	rm -f "$(IQTREE_BINARY)" && \
	echo "User-specific IQ-TREE2 successfully removed." || \
	{ echo "Error: Failed to remove user-specific IQ-TREE2 binary from $(IQTREE_BINARY). Exiting."; exit 1; }; \

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

install: install_iqtree2 venv make_usr_bin 
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
		cp $(VENV_BIN)recur $(USER_INSTALL_DIR)/ && \
		echo "RECUR successfully copied to $(USER_INSTALL_DIR)." || \
		{ echo "Error: Failed to copy RECUR to $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
	elif [ "$$recur_exists" = "1" ]; then \
		recur_path=$$(command -v recur); \
		echo "RECUR already exists at: $$recur_path. Skipping installation."; \
	fi


run: install_iqtree2 install
	@echo "Running RECUR..."
	@if [ -f "$(VENV_BIN)recur" ]; then \
		$(VENV_BIN)recur -f ExampleData -st AA --outgroups ExampleData
	elif [ -f "$(RECUR_DIR)recur" ]; then \
		$(RECUR_DIR)recur -f ExampleData -st AA --outgroups ExampleData
	else \
		echo "RECUR not found. Installing..."; \
		$(MAKE) install; \
	fi

clean:
	@echo "Cleaning build environment..."
	rm -rf **/__pycache__
	rm -rf ./build ./dist ./recur.egg-info

purge:
	@echo "Remove RECUR and it's environment and dependencies..."
	rm -rf $(ENV_NAME) 
	rm -rf **/__pycache__
	rm -rf ./build ./dist ./recur.egg-info

	if [ -f "$(RECUR_DIR)recur" ]; then \
		rm -f "$(RECUR_DIR)recur"
		echo "$(RECUR_DIR)recur including it's environment have been removed..."; \
	else \
		echo "RECUR not found..."; \
	fi


.PHONY: make_usr_bin clean clean_iqtree2 purge clean_conda_venv
