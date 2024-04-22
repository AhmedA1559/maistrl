import multiprocessing
import zstd
import tempfile
import os
import atomInSmiles as AIS
from datetime import datetime
import gc
import resource
import rdkit
from rdkit import Chem
def neutralize_atoms(smile):
    mol = Chem.MolFromSmiles(smile,sanitize=True)
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
    return Chem.MolToSmiles(mol)

def encode_chunk(chunk):
    encoded_chunk = [AIS.encode(neutralize_atoms(smile)) for smile in chunk]
    return encoded_chunk



def split_into_chunks(part_file):
    # Calculate the number of lines in the part file
    with open(part_file, 'r') as f:
        total_lines = sum(1 for _ in f)

    # Calculate the size of each chunk
    chunk_size = total_lines // 16

    # Initialize list to store chunks
    chunks = []

    with open(part_file, 'r') as part:
        chunk = []
        lines_read = 0
        for line in part:
            chunk.append(line.strip())
            lines_read += 1
            # If chunk size reached, add chunk to chunks list and reset chunk
            if lines_read == chunk_size:
                #print(chunk)
                #quit()
                chunks.append(chunk)
                chunk = []
                lines_read = 0

        # Add the remaining lines to the last chunk
        if chunk:
            chunks.append(chunk)

    return chunks

def split_into_parts(input_file):
    # Split input file into 12 parts using NamedTemporaryFile
    with open(input_file, 'r') as f:
        total_lines = sum(1 for _ in f)
        chunk_size = total_lines // 12
        parts = []

        for _ in range(12):
            part_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
            parts.append(part_file.name)

        f.seek(0)  # Reset file pointer
        for part_file in parts:
            with open(part_file, 'w') as part:
                for _ in range(chunk_size):
                    part.write(f.readline())

    return parts

def main():
    input_file = input("Enter the input file: ")
    parts = split_into_parts(input_file)

    global_tokens = set()

    for part_file in parts:
        with open(part_file, 'r') as part:
            chunks = split_into_chunks(part_file)

            pool = multiprocessing.Pool(16)
            encoded_chunks = pool.map(encode_chunk, chunks)
            pool.close()
            pool.join()
            del chunks
            gc.collect()
            memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print("Memory usage after coalescing chunks (GB):", memory_usage / (1024 ** 3))

            encoded_part = [item for sublist in encoded_chunks for item in sublist]
            del encoded_chunks
            memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print("Memory usage after encoding part (GB):", memory_usage / (1024 ** 3))

            # Update global set of encountered tokens
            print("finished encoding...",datetime.now().time())
            for encoded_line in encoded_part:
                   for token in encoded_line.split():     
                            global_tokens.add(token)
            encoded_str = '\n'.join(encoded_part)
            del encoded_part
            gc.collect()
            memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print("Memory usage after encoding string (GB):", memory_usage / (1024 ** 3))



        with open('encoded_part.zst', 'ab') as archive:
            compressed_data = zstd.compress(encoded_str.encode())
            archive.write(compressed_data)
        print("finished writing to archive...",datetime.now().time())
        os.remove(part_file)
        print("part complete",datetime.now().time())


        del encoded_str
        gc.collect()
        memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print("Memory usage after part completion (bytes):", memory_usage)

    with open('vocab.txt', 'w') as vocab_file:
        for token in global_tokens:
            vocab_file.write(token + '\n')
    print("vocab written")
    print("execution complete")


if __name__ == "__main__":
    main()
