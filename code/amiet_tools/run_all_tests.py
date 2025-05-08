#!/usr/bin/env python3
"""
Script to run all test scripts and save results to whyresults folder
"""
import os
import subprocess
import shutil

# Create whyresults directory if it doesn't exist
results_dir = "whyresults"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Get current working directory
current_dir = os.getcwd()

# Path to original scripts
script_paths = [
    "test_scripts/TestScript1_SingleGustPlots.py",
    "test_scripts/TestScript2_MultGustsPlots.py",
    "test_scripts/TestScript3_MultGustsSurfPressure.py"
]

# Create modified versions of the scripts
for script_path in script_paths:
    script_name = os.path.basename(script_path)
    
    # Read the original script
    with open(script_path, 'r') as f:
        content = f.read()
    
    # Disable LaTeX rendering
    content = content.replace("plt.rc('text', usetex=True)", "plt.rc('text', usetex=False)")
    
    # Enable saving figures and redirect to whyresults directory
    content = content.replace("save_fig = False", "save_fig = True")
    
    # Change folder to whyresults
    if "fig_folder = '../plots/'" in content:
        content = content.replace("fig_folder = '../plots/'", f"fig_folder = '{current_dir}/{results_dir}/'")
    else:
        # Add fig_folder variable if not present
        content = content.replace("# flag for saving figures\nsave_fig = True", 
                                 f"# flag for saving figures\nsave_fig = True\n# folder to save figures\nfig_folder = '{current_dir}/{results_dir}/'")
    
    # Change any savefig calls to use the fig_folder
    content = content.replace(".savefig('", ".savefig(fig_folder + '")
    
    # Change output format from eps to png
    content = content.replace(".eps')", ".png')")
    
    # Fix the path to the DARP2016 files to use absolute paths
    content = content.replace("../DARP2016_TestSetup.txt", f"{current_dir}/DARP2016_TestSetup.json")
    content = content.replace("../DARP2016_AirfoilGeom.txt", f"{current_dir}/DARP2016_AirfoilGeom.json")
    
    # Fix np.int deprecation in TestScript3
    content = content.replace("np.int(", "int(")
    
    # Create a temporary modified script in the same directory as the original
    temp_script_path = os.path.join(os.path.dirname(script_path), "temp_" + script_name)
    
    # Write the modified script
    with open(temp_script_path, 'w') as f:
        f.write(content)
    
    print(f"Created modified script: {temp_script_path}")

# Run each modified script
for script_path in script_paths:
    script_name = os.path.basename(script_path)
    temp_script_name = "temp_" + script_name
    temp_script_path = os.path.join(os.path.dirname(script_path), temp_script_name)
    
    print(f"\nRunning {script_name}...")
    try:
        # Run the script
        result = subprocess.run(['python', temp_script_path], 
                               check=True)
        
        print(f"Successfully completed {script_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
    
    # Clean up the temporary script
    os.remove(temp_script_path)

print("\nAll scripts completed. Check the whyresults directory for output files.") 