import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from collections import Counter
from collections import defaultdict
import re
import pandas as pd
import os
import time
from conversion_matrix import get_category_name, get_category_number




def show_error(message):
    """Displays an error message in a popup window."""
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    messagebox.showerror("Error", message)  # Show error popup

def show_info(message):
    """Displays an error message in a popup window."""
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    messagebox.showinfo("info", message)  # Show error popup


def browse_file(sentence):
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_path = filedialog.askopenfilename(
        title=sentence,
        defaultextension=".fasta",
        filetypes=[("fasta files", "*.fasta")]
    )
    return file_path

def all_strings_same_length(strings):
    return len(set(map(len, strings))) == 1

def read_fasta(file_path):
    """Reads a FASTA file and validates all possible errors before processing."""
    if not os.path.exists(file_path):
        show_error(f"Error: The file '{file_path}' does not exist.")
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")

    headers = []
    sequences = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Error: The file '{file_path}' is empty.")
        show_error(f"Error: The file '{file_path}' is empty.")

    sequence_data = ''
    for line in lines:
        line = line.strip()

        if line.startswith(">"):
            if sequence_data:
                sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
                sequences.append(sequence_data)
                sequence_data = ''

            headers.append(line[1:])  # Store header (without '>')

        elif line:
            sequence_data += line


    if sequence_data:
        sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
        sequences.append(sequence_data)

    if not headers:
        show_error(f"Error: The file '{file_path}' does not contain any headers.")
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")

    if len(headers) != len(sequences):
        show_error(f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)}).")
        raise ValueError(f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)}).")

    if not all_strings_same_length(sequences):
        show_error(f"Error: The file '{file_path}' is not a multiple sequence alignment.")
        raise ValueError(f"Error: The file '{file_path}' is not a multiple sequence alignment.")

    return headers, sequences

def consensus_sequence(msa, threshold):
    """Compute consensus sequence from a list of aligned sequences."""
    consensus = []
    for column in zip(*msa):  # Transpose the alignment matrix
        counter = Counter(column)
        del counter['-']
        seq_number = sum(counter.values())
        if len(counter.most_common(1)) != 0:
            ratio = counter.most_common(1)[0][1]/seq_number
            if ratio >= threshold and counter.most_common(1)[0][1] >= 5:
            # if ratio >= threshold :
                most_common = counter.most_common(1)[0][0]  # Get most frequent character
                consensus.append(most_common)
            else:
                consensus.append('-')
        else:
            consensus.append('-')


    return "".join(consensus)

def count_triplets(string, reference):
    triplet_counts = defaultdict(int)
    s = "".join([reference[i] for i in range(len(string)) if string[i] != '-'])
    for i in range(len(s) - 2):  # Ensuring we have enough characters for a triplet
        triplet = s[i:i + 3]  # Extract 3 consecutive characters
        triplet_counts[triplet] += 1

    return triplet_counts

def find_mutations(seq1, seq2):
    mutations = list(filter(lambda i: seq1[i] != seq2[i] and 1 <= i < len(seq1) - 1, range(len(seq1))))
    return (mutations)

def find_mutations_with_dashes(query , reference):
    mutations = list(filter(lambda i: query[i] != reference[i] and 1 <= i < len(query) - 1 and query[i] != '-' and reference[i] != '-', range(len(query))))
    return (mutations)


def extract_valid_pairs(seq1, seq2):
    """Returns the positions where both sequences have nucleotides (no dashes)."""
    result1 = []
    result2 = []

    for a, b in zip(seq1, seq2):
        if a != '-' or b != '-':
            result1.append(a)
            result2.append(b)

    return ''.join(result1), ''.join(result2)

try:
    fasta_path = browse_file("Please choose the FASTA file that contains the multiple sequence alignment")
    base_directory = os.path.splitext(fasta_path)[0]
    viruses , msa = read_fasta(fasta_path)
    consensus = consensus_sequence(msa, 0.6)
    with open(base_directory + "_consensus_sequence.txt", "w") as file:
        file.write(consensus)


    triplet_counts = defaultdict(lambda: defaultdict(int))
    for i in range(len(viruses)):
        print("triplet_counting: {}".format(i))
        # seq = msa[i].replace("-", "")
        seq = msa[i]
        # s = "".join([consensus[i] if seq[i] != consensus[i] and consensus[i] != '-' else seq[i]  for i in range(len(seq)) ])
        s = []
        # for NA_pos in range(len(seq)):
        #     if consensus[NA_pos] != '-' and seq[NA_pos] != '-':
        #         s = s + [consensus[NA_pos]]
        #     elif consensus[NA_pos] == '-' and seq[NA_pos] != '-':
        #         s = s + [seq[NA_pos]]
        # s = "".join(s)
        s = "".join([
            consensus[i] if consensus[i] != '-' and seq[i] != '-'
            else seq[i]
            for i in range(len(seq))
            if seq[i] != '-'
        ])

        for j in range(len(s) - 2):  # Ensuring we have enough characters for a triplet
            triplet = s[j:j + 3]  # Extract 3 consecutive characters
            triplet_counts[viruses[i]][triplet] += 1

    # Convert mutation data to a Pandas DataFrame
    df = pd.DataFrame.from_dict(triplet_counts, orient="index").fillna(0)
    df.to_csv(base_directory + "_triplet_counts.csv")


    mutation_counts = defaultdict(lambda: defaultdict(int))
    mutation_category = defaultdict(set)
    for i in range(len(viruses)):
        print("mutation_counting: {}".format(i))
        seq , m_consensus = extract_valid_pairs(msa[i], consensus)
        mutation_list = find_mutations_with_dashes(seq, m_consensus)
        for mutation_point in mutation_list:
            if "-" not in seq[(mutation_point-1):(mutation_point+2)] and \
                "-" not in m_consensus[(mutation_point-1):(mutation_point+2)]:
                from_seq = m_consensus[(mutation_point - 1):(mutation_point + 2)]
                to_seq = seq[(mutation_point - 1):(mutation_point + 2)]
                reverse_mutation, category_number = get_category_number(from_seq, to_seq)
                mutation_counts[viruses[i]][category_number] += 1

                mutation_category[category_number].add(from_seq + "_" + to_seq)
                mutation_category[category_number].add(reverse_mutation)

    # Convert mutation data to a Pandas DataFrame
    df = pd.DataFrame.from_dict(mutation_counts, orient="index").fillna(0)
    df.to_csv(base_directory + "_mutation_counts.csv")

    # Convert mutation data to a Pandas DataFrame
    df = pd.DataFrame.from_dict(mutation_category, orient="index").fillna(0)
    df.to_csv(base_directory + "_mutation_category.csv")



    show_info("The results are saved in\n triplet_counts.csv and mutation_counts.csv.\n"
              " The consensus sequence is stored in consensus_sequence.txt.")



except Exception as e:
    print(str(e))












