#!/bin/bash

# Simple test script to validate conda detection functionality
# This can be run locally to test the conda setup

echo "========================================"
echo "Testing SpecFWAT Conda Detection"
echo "========================================"

# Test 1: Check if conda is available
echo "Test 1: Checking conda availability..."
if command -v conda &> /dev/null; then
    echo "✅ PASS: Conda found in PATH"
    conda --version
else
    echo "❌ FAIL: Conda not found in PATH"
    echo "Please install conda or run the setup script"
fi

# Test 2: Check if environment.yml is valid
echo ""
echo "Test 2: Validating environment.yml..."
if [ -f "environment.yml" ]; then
    if python -c "import yaml; yaml.safe_load(open('environment.yml'))" 2>/dev/null; then
        echo "✅ PASS: environment.yml is valid YAML"
    else
        echo "❌ FAIL: environment.yml has syntax errors"
    fi
else
    echo "❌ FAIL: environment.yml not found"
fi

# Test 3: Check if setup script is executable
echo ""
echo "Test 3: Checking setup script..."
if [ -f "setup_conda_env.sh" ]; then
    if [ -x "setup_conda_env.sh" ]; then
        echo "✅ PASS: setup_conda_env.sh is executable"
    else
        echo "⚠️  WARNING: setup_conda_env.sh is not executable (run: chmod +x setup_conda_env.sh)"
    fi
    
    # Check syntax
    if bash -n setup_conda_env.sh; then
        echo "✅ PASS: setup_conda_env.sh has valid syntax"
    else
        echo "❌ FAIL: setup_conda_env.sh has syntax errors"
    fi
else
    echo "❌ FAIL: setup_conda_env.sh not found"
fi

# Test 4: Check if GitHub Actions workflows are valid
echo ""
echo "Test 4: Checking GitHub Actions workflows..."
workflow_files=(".github/workflows/conda-detection.yml" ".github/workflows/conda-multiplatform.yml")

for workflow in "${workflow_files[@]}"; do
    if [ -f "$workflow" ]; then
        if python -c "import yaml; yaml.safe_load(open('$workflow'))" 2>/dev/null; then
            echo "✅ PASS: $workflow is valid YAML"
        else
            echo "❌ FAIL: $workflow has syntax errors"
        fi
    else
        echo "❌ FAIL: $workflow not found"
    fi
done

# Test 5: Check if SpecFWAT environment exists (if conda is available)
echo ""
echo "Test 5: Checking for existing SpecFWAT environment..."
if command -v conda &> /dev/null; then
    if conda env list | grep -q "specfwat"; then
        echo "✅ PASS: SpecFWAT conda environment found"
        echo "Available SpecFWAT environments:"
        conda env list | grep specfwat
    else
        echo "ℹ️  INFO: No SpecFWAT conda environment found (this is expected on first run)"
    fi
else
    echo "ℹ️  INFO: Skipping environment check (conda not available)"
fi

echo ""
echo "========================================"
echo "Test Summary Complete"
echo "========================================"
echo ""
echo "To set up the SpecFWAT conda environment:"
echo "  ./setup_conda_env.sh"
echo ""
echo "To create environment manually:"
echo "  conda env create -f environment.yml"
echo "  conda activate specfwat-env"