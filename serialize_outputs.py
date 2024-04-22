import os
import json
import datetime
from rdkit import Chem
def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def extract_score_and_smiles(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    score = None
    smiles = None
    for line in lines:
        if line.startswith('REMARK VINA RESULT:'):
            score = float(line.split()[3])
        elif line.startswith('REMARK SMILES'):
            smiles = ' '.join(line.split()[2:])
            break
    
    if smiles and score is not None:
        mol = neutralize_atoms(Chem.MolFromSmiles(smiles, sanitize=True))
        if mol:
            sanitized_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            return sanitized_smiles, score
        else:
            return None, None

def process_directory(directory):
    results = {}
    for filename in os.listdir(directory):
        if filename.endswith(".pdbqt"):
            file_path = os.path.join(directory, filename)
            smiles, score = extract_score_and_smiles(file_path)
            if smiles and score is not None:
                results[smiles] = score
    return results

def serialize_results(results):
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    output_filename = f"results_{timestamp}.json"
    with open(output_filename, 'w') as f:
        json.dump(results, f, indent=4)
    print(f"Results serialized to {output_filename}")

if __name__ == "__main__":
    directory = "Ligands_out"  # Replace with your directory path
    results = process_directory(directory)
    serialize_results(results)
