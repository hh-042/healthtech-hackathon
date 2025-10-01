import google.generativeai as genai
import pandas as pd
import random
from Bio.SeqUtils import MeltingTemp as mt

import os

# Set API key
genai.configure(api_key="")


def generate_dna_sequences(num_sequences=5):
    """
    Generates DNA sequences using Google Gemini API with specified parameters.

    Parameters:
        - num_sequences (int): Number of sequences to generate.

    Returns:
        - List of tuples (sequence, melting_temp, P90, N10, Diff).
    """

    prompt = f"""
    Generate exactly {num_sequences} DNA template sequences.  
    Each sequence must:
    - Be **on a new line** (one per line).
    - Follow this structure:
      - Starts with 10-15 base sequence
      - Followed by 4 random A,T,C,G bases
      - Followed by 5 specific cutting bases ("GACTC")
      - Followed by a 10-15 base sequence that matches the start
    - Ensure Tm (melting temperature) is between 50-60°C.
    - Ensure Diff (N10 - P90) is positive and large.
    **Output:** Provide only the sequences, no extra text, formatted with one sequence per line.
    """

    model = genai.GenerativeModel("gemini-1.5-pro-latest")
    response = model.generate_content(prompt)

    if not response or not response.text:
        return "Failed to generate sequences."
    print("Raw API response:\n", response.text)  # Debugging

    sequences = response.text.strip().split("\n")

    filtered_sequences = []
    for seq in sequences:
        seq = seq.strip()
        if 29 <= len(seq) <= 39:
            melting_temp = mt.Tm_NN(seq)  # Calculate Tm
            if 50 <= melting_temp <= 60:
                P90 = random.uniform(5, 15)  # Simulate P90 (lower is better)
                N10 = random.uniform(P90 + 5, P90 + 20)  # Simulate N10 (higher is better)
                Diff = N10 - P90
                if Diff > 5:  # Ensures good separation
                    filtered_sequences.append((seq, melting_temp, P90, N10, Diff))

    return filtered_sequences


# Example usage
if __name__ == "__main__":
    sequences = generate_dna_sequences(num_sequences=5)


