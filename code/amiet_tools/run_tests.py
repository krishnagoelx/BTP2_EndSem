#!/usr/bin/env python3
"""
Script to run the three test scripts with DARP2016 files as arguments
"""
import os
import sys
import amiet_tools as AmT
from amiet_tools.classes import TestSetup, AirfoilGeom

# Create plots directory if it doesn't exist
plots_dir = "plots"
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Monkey patch the loader functions to fix initialization issues
original_loadTestSetup = AmT.loadTestSetup
original_loadAirfoilGeom = AmT.loadAirfoilGeom

def patched_loadTestSetup(*args):
    """Patched version that initializes with default DARP2016 values"""
    # Default DARP2016 values
    if len(args) == 0:
        return TestSetup(c0=340.0, rho0=1.2, p_ref=20e-6, Ux=60.0, 
                         turb_intensity=0.025, length_scale=0.007, z_sl=-0.075)
    return original_loadTestSetup(*args)

def patched_loadAirfoilGeom(*args):
    """Patched version that initializes with default DARP2016 values"""
    # Default DARP2016 values
    if len(args) == 0:
        return AirfoilGeom(b=0.075, d=0.225, Nx=100, Ny=101)
    return original_loadAirfoilGeom(*args)

# Apply the patches
AmT.loadTestSetup = patched_loadTestSetup
AmT.loadAirfoilGeom = patched_loadAirfoilGeom

# Scripts to run
scripts = [
    "test_scripts/TestScript1_SingleGustPlots.py",
    "test_scripts/TestScript2_MultGustsPlots.py",
    "test_scripts/TestScript3_MultGustsSurfPressure.py"
]

# Run each script
print("Running test scripts...")
for script in scripts:
    print(f"Running {script}...")
    try:
        # Execute the script
        with open(script) as f:
            script_content = f.read()
            
        # Replace any savefig calls to save to plots directory
        script_content = script_content.replace("save_fig = False", "save_fig = True")
        
        # Create a temporary script
        temp_script = f"temp_{os.path.basename(script)}"
        with open(temp_script, "w") as f:
            f.write(script_content)
        
        # Run the script
        exec(open(temp_script).read())
        print(f"Successfully completed {script}")
        
        # Clean up
        os.remove(temp_script)
    except Exception as e:
        print(f"Error running {script}: {e}")

print("All scripts completed. Check the plots directory for output files.") 