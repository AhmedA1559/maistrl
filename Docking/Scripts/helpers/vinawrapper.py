import subprocess
import os
os.chdir(os.getcwd() + "/Vina-GPU+/")

def run_vina_gpu(receptor, ligand_directory, center, size, opencl_binary_path=None):
    # Build the command with the provided options
    command = [
        './QuickVina2-GPU-2-1',
        '--receptor', receptor,
        '--ligand_directory', ligand_directory,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(size),
        '--size_y', str(size),
        '--size_z', str(size),
        '--thread', '8000',
        '--opencl_binary_path','/home/victor/DOCKING/Vina-GPU+'
    ]


    # Run the Vina-GPU-2.0 program
    with open('result.txt', 'w') as output_file:
        subprocess.run(command, stdout=output_file, stderr=subprocess.STDOUT)

# Set the specified options
receptor = r'../../../d2dr.pdbqt'
ligand_directory = r'../../../Ligands'
center = (9.74, 5.37, -10.25)

size = 16

# Specify the path to the OpenCL binary (if known)
opencl_binary_path = None

# Run Vina-GPU with the specified options
run_vina_gpu(receptor, ligand_directory, center, size, opencl_binary_path)

