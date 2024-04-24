# Reproducing MAISTRL Experiments

This README outlines the steps to reproduce a set of machine learning experiments from our GitHub repository.

## Overview

The experiment consists of four main steps:
1. Data preparation for pre-training.
2. Pre-training the model.
3. Fine-tuning the pre-trained model.
4. Evaluating the fine-tuned experiments.

## Data Preparation

### Downloading Zinc Data

1. Navigate to [Zinc20 Data](https://zinc20.docking.org/tranches/home/).
2. Click on "3D".
3. React "Anodyne".
4. Set Purch to "Wait OK", ph to "Ref Mid", Charge to "0 +1".
5. Select all options on the grid from molecular weight 250-450 daltons and logp -1 to -5.
6. Download the selected files.
7. Concatenate the downloaded files into one text file.
8. Randomly shuffle the lines.
9. Extract the first 200 million lines, saving 1 million for fine-tuning. We named these files `train.txt` and `reserved.txt`.

### Convert Data to AIS Format

1. In the Pre-training folder, run `multithreaded_AIS_converter.py` on both `train.txt` and `reserved.txt`.
2. Follow the instructions provided by the script.
3. Concatenate the generated vocab files into `vocab.txt`, avoiding duplication.
4. Unzip the outputted files to get the encoded training textfile and reserved textfile and name it train.text

### Pre-training

1. Run `pretrain_maistral.ipynb` using the generated training file containing molecules in AIS form (`train.txt`) and the vocab file`vocab.txt`.

## Fine-tuning

### Docking Evaluation

1. Dock the molecules in the original unencoded `reserved.txt` to evaluate their binding to the selected molecular target (D2DR) using `dock_data_for_finetuning.py` to evaluate good vs bad molecules for the purposes of finetuning with SFT and DPO. This expects a file "smile.txt" in the docking folder. Make a copy of the original reserved.txt, rename it to smiles.txt and place it here. 
2. This will output a JSON file of molecules and their docking scores.

### Fine-tuning Process

1. Run `sft_trainer.ipynb` on the generated JSON file from Docking Evaluation along with `vocab.txt`.
2. The notebook will output text files with hyperparameters as filenames, along with metrics such as uniqueness, novelty, and internal diversity appended to the filename. The contents will be a set of 250 sampled outputs. 

### DPO Fine-tuning

1. Construct pairs of molecules for DPO:
   - To construct random pairs, run `random_pairing.py` with `vocab.txt` and the JSON file used in `sft_trainer.ipynb`.
   - To construct pairs based on molecular similarity, run `create_similar_pairs.ipynb` using the same JSON file. 
2. Each method will output a JSON file containing pairs that, along with `vocab.txt`, can be used with `dpo_trainer.ipynb`.


### Results Replication

1. Choose all dockers with above 98% uniqueness, and place them in folders based on the experiment (we ran DPO-similar-pairs, DPO-similar-half-pairs, DPO-random-pairs, DPO-random-halfdata-pairs, SFT, and Baseline (non finetuned model). Make sure each folder name starts with dpo_ or sft_. Place baseline in any of the folders. 
2. Run `prep.py` in the finetuning folder to append validity to the names in each folder. 
3. On all hyperparameter combinations with above 90% validity, keep them in their folder starting with dpo_ or sft_ and move all the folders to a folder containing `dock_experiments.py` from `Docking/Scripts` and run it to get the final experiment results.

Here are the hyperparameters we found to work the best:
| Method             | # Epochs| Beta    | Alpha   |
|--------------------|---------|---------|---------|
| Baseline           | n/a     | n/a     | n/a     |
| DPO                | 3       | 0.5     | 0.00005 |
| DPO half           | 2       | 0.3     | 0.0001  |
| DPO random         | 1       | 0.5     | 0.00005 |
| DPO random half    | 1       | 0.5     | 0.00005 |
| SFT                | 3       | n/a     | 0.0005  |
