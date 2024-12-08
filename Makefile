.ONESHELL:

.DEFAULT_GOAL := install

ifndef MAKE_VERSION
$(error This Makefile requires GNU Make. Please use gmake instead of make.)
endif

SYSTEM_WIDE ?= 0
HOME_DIR := $(if $(HOME),$(HOME),$(shell echo ~))

PROMPT_USER_INSTALL_DIR = \
	@if [ -d "$(HOME_DIR)/bin" ]; then \
		read -p "Directory $(HOME_DIR)/bin exists. Do you want to use it? (y/n) " choice; \
		case $$choice in \
			[yY]*) echo "$(HOME_DIR)/bin";; \
			[nN]*) \
				read -p "Enter a new directory name (relative to $(HOME_DIR)): " new_dir; \
				USER_INSTALL_DIR=$(HOME_DIR)/$$new_dir; \
				echo $$USER_INSTALL_DIR; \
				mkdir -p $$USER_INSTALL_DIR || { echo "Error creating directory $$USER_INSTALL_DIR. Exiting."; exit 1; }; \
				;; \
			*) echo "Invalid choice. Exiting."; exit 1; \
		esac; \
	else \
		echo "$(HOME_DIR)/bin"; \
	fi

USER_INSTALL_DIR := $(shell $(PROMPT_USER_INSTALL_DIR))
SYSTEM_INSTALL_DIR := /usr/local/bin
RECUR_DIR := $(USER_INSTALL_DIR)

ifeq ($(SYSTEM_WIDE),1)
    IQTREE_INSTALL_DIR := $(SYSTEM_INSTALL_DIR)
    SUDO_PREFIX := $(if $(shell [ "$(USER)" = "root" ] || echo 1),sudo)
else
    IQTREE_INSTALL_DIR := $(USER_INSTALL_DIR)
    SUDO_PREFIX :=
endif

IQTREE_BINARY := $(IQTREE_INSTALL_DIR)/iqtree2

# URLs for IQ-TREE2 binaries
IQTREE_URL_LINUX_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-arm.tar.gz
IQTREE_URL_LINUX_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz

IQTREE_URL_MAC_UNIVERSAL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS.zip
IQTREE_URL_MAC_INTEL := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-intel.zip
IQTREE_URL_MAC_ARM := https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-macOS-arm.zip

PYTHON := ./.venv/bin/python3
PIP := ./.venv/bin/pip
VENV_BIN := ./.venv/bin/

make_usr_bin:
	@if [ ! -d "$(USER_INSTALL_DIR)" ]; then \
		echo "Directory $(USER_INSTALL_DIR) does not exist. Creating it..."; \
		mkdir -p $(USER_INSTALL_DIR); \
	fi

install_iqtree: make_usr_bin
	@if [ "$(FORCE)" = "1" ] || [ ! -f "$(IQTREE_BINARY)" ]; then \
		echo "Installing IQ-TREE2 for Linux/MacOS..."; \
		if [ "$(shell uname)" = "Darwin" ]; then \
			if [ "$(shell uname -m)" = "arm64" ]; then \
				IQTREE_URL=$(IQTREE_URL_MAC_ARM); \
			elif [ "$(shell uname -m)" = "x86_64" ]; then \
				IQTREE_URL=$(IQTREE_URL_MAC_INTEL); \
			else \
				IQTREE_URL=$(IQTREE_URL_MAC_UNIVERSAL); \
			fi; \
		else \
			if [ "$(shell uname -m)" = "aarch64" ]; then \
				IQTREE_URL=$(IQTREE_URL_LINUX_ARM); \
			else \
				IQTREE_URL=$(IQTREE_URL_LINUX_INTEL); \
			fi; \
		fi; \
		curl -L $$IQTREE_URL -o iqtree2.tar.gz || { echo "Error: Failed to download IQ-TREE2. Exiting."; exit 1; }; \
		mkdir -p "$(IQTREE_INSTALL_DIR)" || { echo "Error: Failed to create directory $(IQTREE_INSTALL_DIR). Exiting."; exit 1; }; \
		if [ ! -w "$(IQTREE_INSTALL_DIR)" ]; then \
			echo "Error: Directory $(IQTREE_INSTALL_DIR) is not writable. Check permissions or run with appropriate privileges."; \
			exit 1; \
		fi; \
		$(SUDO_PREFIX) tar --strip-components=2 -xzf iqtree2.tar.gz -C "$(IQTREE_INSTALL_DIR)" || { echo "Error: Failed to extract IQ-TREE2. Exiting."; exit 1; }; \
		rm iqtree2.tar.gz || { echo "Error: Failed to remove temporary file iqtree2.tar.gz. Exiting."; exit 1; }; \
		echo "IQ-TREE2 installed to $(IQTREE_INSTALL_DIR)."; \
		if [ "$(SYSTEM_WIDE)" = "1" ]; then \
			echo "System-wide installation complete. Binary available in $(IQTREE_INSTALL_DIR)."; \
		else \
			echo "Adding $(IQTREE_INSTALL_DIR) to the PATH permanently."; \
			if ! grep -q 'export PATH="$(IQTREE_INSTALL_DIR):$$PATH"' ~/.bashrc; then \
				echo 'export PATH="$(IQTREE_INSTALL_DIR):$$PATH"' >> ~/.bashrc || { echo "Error: Failed to update PATH in ~/.bashrc. Exiting."; exit 1; }; \
			fi; \
			echo "PATH update added to ~/.bashrc."; \
		fi; \
	else \
		echo "IQ-TREE2 already exists in $(IQTREE_INSTALL_DIR). Skipping installation."; \
	fi


clean_iqtree2:
	@if [ "$(SYSTEM_WIDE)" = "1" ]; then \
		echo "Cleaning system-wide IQ-TREE2 installation..."; \
		$(SUDO_PREFIX) rm -f "$(IQTREE_BINARY)" && \
		echo "System-wide IQ-TREE2 successfully removed." || \
		{ echo "Error: Failed to remove system-wide IQ-TREE2 binary from $(IQTREE_BINARY). Exiting."; exit 1; }; \
	else \
		echo "Cleaning user-specific IQ-TREE2 installation..."; \
		rm -f "$(IQTREE_BINARY)" && \
		echo "User-specific IQ-TREE2 successfully removed." || \
		{ echo "Error: Failed to remove user-specific IQ-TREE2 binary from $(IQTREE_BINARY). Exiting."; exit 1; }; \
	fi

.venv/bin/activate: requirements.txt
	@echo "Creating virtual environment for RECUR..."
	python3 -m venv .venv
	chmod +x .venv/bin/activate
	. .venv/bin/activate
	@echo "Virtual environment created and dependencies installed."

venv: .venv/bin/activate
	@echo "Virtual environment activated."

install: venv make_usr_bin
	@echo "Installing RECUR..."
	$(PIP) install -e . || { echo "Error: Failed to install RECUR. Exiting."; exit 1; }
	@if [ "$(SYSTEM_WIDE)" = "1" ]; then \
		if [ -d "/usr/local/bin" ]; then \
			echo "Copying RECUR to $(SYSTEM_INSTALL_DIR)..."; \
			$(SUDO_PREFIX) cp $(VENV_BIN)recur $(SYSTEM_INSTALL_DIR) && \
			echo "RECUR successfully copied to $(SYSTEM_INSTALL_DIR)." || \
			{ echo "Error: Failed to copy RECUR to $(SYSTEM_INSTALL_DIR). Exiting."; exit 1; }; \
		else \
			echo "System-wide directory $(SYSTEM_INSTALL_DIR) not found."; \
			exit 1; \
		fi; \
	else \
		mkdir -p $(USER_INSTALL_DIR) || { echo "Error: Failed to create directory $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
		echo "Copying RECUR to $(USER_INSTALL_DIR)..."; \
		cp $(VENV_BIN)recur $(USER_INSTALL_DIR)/ && \
		echo "RECUR successfully copied to $(USER_INSTALL_DIR)." || \
		{ echo "Error: Failed to copy RECUR to $(USER_INSTALL_DIR). Exiting."; exit 1; }; \
	fi

install_dependencies: venv 
	$(PIP) install -r requirements.txt

run: install install_iqtree
	@echo "Running RECUR..."
	@if [ -f "$(VENV_BIN)recur" ]; then \
		$(VENV_BIN)recur -f ExampleData -st AA --outgroups ExampleData
	elif [ -f "$(RECUR_DIR)recur" ]; then \
		$(RECUR_DIR)recur -f ExampleData -st AA --outgroups ExampleData
	elif [ -f "$(SYSTEM_INSTALL_DIR)recur" ]; then \
		$(SUDO_PREFIX) $(SYSTEM_INSTALL_DIR)recur -f ExampleData -st AA --outgroups ExampleData
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
	rm -rf .venv
	rm -rf **/__pycache__
	rm -rf ./build ./dist ./recur.egg-info
	@if [ -f "$(VENV_BIN)recur" ]; then \
	    rm -f "$(VENV_BIN)recur"
		echo "$(VENV_BIN)recur including it's environment and dependencies have been removed..."; \
	if [ -f "$(RECUR_DIR)recur" ]; then \
		rm -f "$(RECUR_DIR)recur"
		echo "$(RECUR_DIR)recur including it's environment and dependencies have been removed..."; \
	elif [ -f "$(SYSTEM_INSTALL_DIR)recur" ]; then \
		$(SUDO_PREFIX) rm -f $(SYSTEM_INSTALL_DIR)recur
		echo "$(SYSTEM_INSTALL_DIR)recur including it's environment and dependencies have been removed..."; \
	else \
		echo "RECUR not found..."; \
	fi

.PHONY: make_usr_bin clean clean_iqtree2 purge
