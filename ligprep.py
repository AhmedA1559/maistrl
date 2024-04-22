import os
import sys
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
import shutil
import concurrent.futures
import subprocess

def protonate_and_embed_with_obabel(smiles):
    try:
        # Run obabel command to protonate
        obabel_command = f'obabel -:"{smiles}" -ismi -ocan -p7'
        obabel_output = subprocess.check_output(obabel_command, shell=True, text=True).strip()

        # Collect the results into a list
        results = []
        for state_idx, protonated_smiles in enumerate(obabel_output.split('\n')):
            mol = Chem.MolFromSmiles(protonated_smiles)
            AllChem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
            results.append((mol, state_idx))

        return results

    except Exception as e:
        #print(f"Error processing SMILES: {smiles}")
        #print(f"Error message: {e}")
        #print(f"SMILES of failed molecule: {smiles}")
        return []
def save_to_pdbqt(mol_setup, directory, filename):
    # Create the output directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)

    # Combine directory and filename to get the full path
    full_path = os.path.join(directory, filename)


    # Generate PDBQT string using Meeko PDBQTWriterLegacy
    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setup)

    # Check if the PDBQT string is not empty#






    if pdbqt_string.strip():
        # Write the PDBQT string to a file
        with open(full_path, 'w') as file:
            file.write(pdbqt_string)
        
        #print(f'File {full_path} written successfully.')
    else:
        #print(f'Error: Empty PDBQT string for {filename}. Skipping.')
        pass      

def process_molecule(mol):
    try:
        # Prepare the molecule using Meeko
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        # Process each Meeko MoleculeSetup
        for setup in mol_setups:
            # Save the PDBQT string to a file with a unique filename for each protonation state
            pdbqt_filename = f'molecule_{idx}_{state_idx}.pdbqt'
            save_to_pdbqt(setup, default_output_directory, pdbqt_filename)

    except Exception as e:
        pass
        #print(f"Error preparing molecule: {e}")
        #print(f"Skipping molecule.")
        
if __name__ == "__main__":
    # Check if the input file is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python ligprep.py input_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]

    with open(input_file, "r") as file:
        # Read SMILES strings from the file
        smiles_list = [line.strip() for line in file]

    # Specify the default output directory as a subdirectory named "PDBQT_Output"
    default_output_directory = "Ligands"
    if os.path.exists(default_output_directory):
        shutil.rmtree(default_output_directory)
    if os.path.exists("Ligands_out"):
        shutil.rmtree("Ligands_out")

    # Process each SMILES string using concurrent.futures
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for idx, results in enumerate(executor.map(protonate_and_embed_with_obabel, smiles_list)):
            for mol, state_idx in results:
                if mol is not None:
                    AllChem.SanitizeMol(mol)
                    process_molecule(mol)
                
                
                


