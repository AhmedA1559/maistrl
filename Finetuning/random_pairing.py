import os
import json
import random
import torch
from transformers import AutoModel, AutoTokenizer
from sklearn.preprocessing import normalize

model = AutoModel.from_pretrained("ibm/MoLFormer-XL-both-10pct", deterministic_eval=True, trust_remote_code=True).to("cuda")
tokenizer = AutoTokenizer.from_pretrained("ibm/MoLFormer-XL-both-10pct", trust_remote_code=True)

# Function to read JSON files from a folder
def read_json_files(folder_path):
    data = {}
    for filename in os.listdir(folder_path):
        if filename.endswith('.json'):
            with open(os.path.join(folder_path, filename), 'r') as file:
                json_data = json.load(file)
                data.update(json_data)
    return data

# Read JSON files from the 'data' folder
folder_path = 'data'
json_data = read_json_files(folder_path)

# Sort the JSON data by score
sorted_data = sorted(json_data.items(), key=lambda x: x[1])

# Calculate the indices for the most negative 3% and most positive 15%
num_entries = len(sorted_data)
negative_threshold_index = int(0.03 * num_entries)
positive_threshold_index = int(0.85 * num_entries)  # 100% - 15% = 85%

# Extract entries with scores in the most negative 3%
negative_entries = [entry for entry in sorted_data[:negative_threshold_index]]

# Extract entries with scores in the most positive 15%
positive_entries = [entry for entry in sorted_data[positive_threshold_index:]]

# Shuffle the positive_entries list, then truncate negative
random.shuffle(positive_entries)
positive_entries = positive_entries[:len(negative_entries)]

# Function to get embeddings in batches
# Process in batches of 
def get_all_embeddings(entries):
  embeddings = []
  batch_size = 512
  for i in range(0, len(entries), batch_size):
    
    batch = entries[i:i + batch_size]
    batch = [sublist[0] for sublist in batch]

    # Tokenize the current batch
    encoded_batch = tokenizer(batch, padding='max_length', return_tensors='pt', truncation=True).to("cuda")
    
    with torch.no_grad():  # Don't compute gradients to save memory and computation
        model_output = model(**encoded_batch)

        # Assuming the model has a pooler_output attribute (like BERT)
        # Convert to float16 to save memory (if necessary)
        batch_embeddings = model_output.pooler_output.to(torch.float16).squeeze().tolist()
    
    print("Finished embedding", i)
    # Extend the embeddings list with the embeddings from the current batch
    embeddings.extend(batch_embeddings)
  return embeddings



# Get embeddings for negative and positive entries
negative_embeddings = get_all_embeddings(negative_entries)
positive_embeddings = get_all_embeddings(positive_entries)

# L2 normalize all the embeddings
negative_embeddings = normalize(negative_embeddings,norm='l2')
positive_embeddings = normalize(positive_embeddings,norm='l2')

for i in range(len(negative_entries)):
  negative_entries[i] = [negative_entries[i],negative_embeddings[i]]
for i in range(len(positive_entries)):
  positive_entries[i] = [positive_entries[i],positive_embeddings[i]]


# Construct pairs and write them to a JSON file
pairs = []
for neg_item, pos_item in zip(negative_entries, positive_entries):
    neg_key, neg_embedding = neg_item
    pos_key, pos_embedding = pos_item
    similarity_score = (neg_embedding * pos_embedding).sum()  # Dot product as similarity score
    pairs.append([neg_key[0], pos_key[0], similarity_score])
similarity_scores = [pair[2] for pair in pairs]
mean_similarity_score = sum(similarity_scores) / len(similarity_scores)
std_similarity_score = torch.tensor(similarity_scores).std().item()

# Print mean and standard deviation
print(f"Mean of similarity scores: {mean_similarity_score}")
print(f"Standard deviation of similarity scores: {std_similarity_score}")
# Write pairs to a JSON file
output_file = 'pairs.json'
with open(output_file, 'w') as file:
    json.dump(pairs, file, indent=4)

print(f"Pairs written to {output_file}")
