import subprocess
import shutil
import os
# Define the input file and chunk size
input_file = "smiles.txt"
chunk_size = 2000

# Function to process a chunk of lines
def process_chunk(chunk):
    with open("temp_input.txt", "w") as f:
        f.writelines(chunk)

    subprocess.run(["python3", "helpers/ligprep.py", "temp_input.txt"],stderr=subprocess.DEVNULL)
    subprocess.run(["python3", "helpers/vinawrapper.py"])
    subprocess.run(["python3", "helpers/serialize_outputs.py"])
    shutil.rmtree("../Ligands")
    shutil.rmtree("../Ligands_out")
    os.mkdir("../Ligands")
    os.mkdir("../Ligands_out")
# Read all lines into a list
with open(input_file, "r") as f:
    lines = f.readlines()

# Process chunks of lines
for i in range(0, len(lines), chunk_size):
    chunk = lines[i:i+chunk_size]
    process_chunk(chunk)
