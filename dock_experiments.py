import os
import shutil
import subprocess
import json
import ast

def extract_info_from_pdbqt(file_path):
    smiles_scores = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                score = float(parts[3])
                smiles = None
                for sub_line in file:
                    if sub_line.startswith("REMARK SMILES"):
                        smiles = sub_line.split(maxsplit=1)[1].strip()
                        break
                if smiles:
                    smiles_scores[smiles] = score
    return smiles_scores


def process_folder(folder_path):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".txt"):
            input_file_path = os.path.join(folder_path, file_name)
            with open(input_file_path, 'r') as input_file:
                contents = input_file.read()
                smiles_set = ast.literal_eval(contents)
            with open("temp_input.txt", "w") as temp_file:
                for smile in smiles_set:
                    temp_file.write(smile + "\n")
            
            subprocess.run(["python3", "ligprep.py", "temp_input.txt"])
            subprocess.run(["python3", "vinawrapper.py"])
            
            smiles_scores = {}
            for pdbqt_file in os.listdir("Ligands_out"):
                if pdbqt_file.endswith(".pdbqt"):
                    smiles_scores.update(extract_info_from_pdbqt(os.path.join("Ligands_out", pdbqt_file)))
            
            mean_score = sum(smiles_scores.values()) / len(smiles_scores)
            median_score = sorted(smiles_scores.values())[len(smiles_scores) // 2]
            std_dev = (sum((x - mean_score) ** 2 for x in smiles_scores.values()) / len(smiles_scores)) ** 0.5
            max_score = min(smiles_scores.values())
            
            output_folder = os.path.join("final_experiment_molecules", "results")
            os.makedirs(output_folder, exist_ok=True)
            
            output_file_name = f"{file_name.replace('.txt', '')}_mean_{mean_score:.2f}_median_{median_score:.2f}_std_{std_dev:.2f}_max_{max_score:.2f}.json"
            output_file_path = os.path.join(output_folder, output_file_name)
            with open(output_file_path, "w") as output_file:
                json.dump(smiles_scores, output_file, indent=4)

            # Clear contents of directories "Ligands" and "Ligands_out"
            shutil.rmtree("Ligands")
            shutil.rmtree("Ligands_out")
            os.mkdir("Ligands")
            os.mkdir("Ligands_out")

# Get folders starting with "dpo_" or "sft_" inside the "final_experiment_results" folder
experiment_folder = "final_experiment_molecules"
folders_to_process = [os.path.join(experiment_folder, folder) for folder in os.listdir(experiment_folder) if folder.startswith("dpo_") or folder.startswith("sft_")]
for folder in folders_to_process:
    process_folder(folder)
