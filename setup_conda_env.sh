#!/bin/bash

# SpecFWAT Conda Detection and Setup Script
# This script detects conda installation and sets up the SpecFWAT environment

set -e

echo "========================================================"
echo "SpecFWAT Conda Detection and Environment Setup"
echo "========================================================"

# Function to print colored output
print_status() {
    local status=$1
    local message=$2
    if [ "$status" = "success" ]; then
        echo "✅ $message"
    elif [ "$status" = "error" ]; then
        echo "❌ $message"
    elif [ "$status" = "info" ]; then
        echo "ℹ️  $message"
    elif [ "$status" = "warning" ]; then
        echo "⚠️  $message"
    fi
}

# Check if conda is available
check_conda() {
    print_status "info" "Checking for conda installation..."
    
    if command -v conda &> /dev/null; then
        print_status "success" "Conda found in system PATH"
        echo "Conda version: $(conda --version)"
        echo "Conda info:"
        conda info --envs
        return 0
    else
        print_status "error" "Conda not found in system PATH"
        return 1
    fi
}

# Install miniconda if not found
install_miniconda() {
    print_status "info" "Installing Miniconda..."
    
    # Detect OS
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    else
        print_status "error" "Unsupported operating system: $OSTYPE"
        exit 1
    fi
    
    # Download and install
    wget "$MINICONDA_URL" -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda3"
    export PATH="$HOME/miniconda3/bin:$PATH"
    conda init bash
    rm miniconda.sh
    
    print_status "success" "Miniconda installed successfully"
}

# Create SpecFWAT environment
create_environment() {
    print_status "info" "Creating SpecFWAT conda environment..."
    
    if [ -f "environment.yml" ]; then
        print_status "info" "Using environment.yml file"
        conda env create -f environment.yml
    else
        print_status "info" "Creating environment with basic dependencies"
        conda create -n specfwat-env -y python=3.11 \
            cmake make gfortran openmpi hdf5 hdf5-parallel \
            openblas lapack yaml-cpp numpy scipy matplotlib
    fi
    
    print_status "success" "SpecFWAT environment created"
}

# Activate and test environment
test_environment() {
    print_status "info" "Testing SpecFWAT environment..."
    
    # Activate environment
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate specfwat-env
    
    # Test key dependencies
    echo "Testing dependencies:"
    
    # CMake
    if command -v cmake &> /dev/null; then
        print_status "success" "CMake: $(cmake --version | head -n1)"
    else
        print_status "error" "CMake not found"
    fi
    
    # Fortran compiler
    if command -v gfortran &> /dev/null; then
        print_status "success" "Fortran compiler: $(gfortran --version | head -n1)"
    else
        print_status "error" "Fortran compiler not found"
    fi
    
    # MPI
    if command -v mpirun &> /dev/null; then
        print_status "success" "MPI: $(mpirun --version | head -n1)"
    else
        print_status "error" "MPI not found"
    fi
    
    # Python packages
    python -c "import numpy; print('✅ NumPy:', numpy.__version__)" 2>/dev/null || print_status "error" "NumPy not available"
    python -c "import scipy; print('✅ SciPy:', scipy.__version__)" 2>/dev/null || print_status "error" "SciPy not available"
    python -c "import h5py; print('✅ h5py:', h5py.__version__)" 2>/dev/null || print_status "warning" "h5py not available"
}

# Main execution
main() {
    echo "Starting conda detection and setup process..."
    
    # Check if conda is available
    if ! check_conda; then
        print_status "warning" "Conda not found. Would you like to install Miniconda? (y/n)"
        read -r response
        if [[ "$response" =~ ^[Yy]$ ]]; then
            install_miniconda
            # Reload shell environment
            source "$HOME/.bashrc" || true
        else
            print_status "error" "Conda is required for SpecFWAT environment setup"
            exit 1
        fi
    fi
    
    # Check if environment already exists
    ENV_NAME="specfwat-env"
    if conda env list | grep -q "$ENV_NAME"; then
        print_status "warning" "Environment '$ENV_NAME' already exists. Would you like to recreate it? (y/n)"
        read -r response
        if [[ "$response" =~ ^[Yy]$ ]]; then
            conda env remove -n "$ENV_NAME" -y
            create_environment
        else
            print_status "info" "Using existing environment"
        fi
    else
        create_environment
    fi
    
    # Test the environment
    test_environment
    
    echo "========================================================"
    print_status "success" "SpecFWAT conda environment setup complete!"
    echo ""
    echo "To activate the environment, run:"
    echo "    conda activate $ENV_NAME"
    echo ""
    echo "To build SpecFWAT, run:"
    echo "    conda activate $ENV_NAME"
    echo "    mkdir build && cd build"
    echo "    cmake .. -DCMAKE_BUILD_TYPE=Release"
    echo "    make -j\$(nproc)"
    echo "========================================================"
}

# Run main function
main "$@"