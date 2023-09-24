#!/bin/bash

# Activate the base conda environment.
source activate base

# Remove the existing conda environment.
conda env remove -n coredump

# Create the new conda environment from the environment.yaml file.
conda env update --file environment.yaml

# Install the kernel for Jupyter Notebook.
conda run -n coredump python -m ipykernel install --user --name coredump --display-name "coredump-kernel"

# Build the Python package.
conda run -n coredump python -m build
