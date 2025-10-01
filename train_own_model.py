#evaluate a Random Forest Classifier to predict the class labels (CI, CII, CIII) of
# DNA sequences based on their sequence composition and some numerical features (P90, N10, Diff)

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

#loading data in
df = pd.read_csv('/Users/angelyao/PycharmProjects/healthHackathon/healthAIHack/384_sequences.csv')
#print(df.head())

# Remove rows with invalid results (e.g., NR, CE, SD, etc.)
df_cleaned = df[~df['Results'].isin(['NR', 'CE', 'SD', 'LA', 'NA', 'NS'])]

# Remove rows with zero values in P90, N10, or Diff
df_cleaned = df_cleaned[(df_cleaned[['P90', 'N10', 'Diff']] != 0).all(axis=1)]

# Map class labels (CI, CII, CIII) to numeric values -1 is good, 2 soso, 3 bad
class_mapping = {'CI': 1, 'CII': 3, 'CIII': 2}
df_cleaned['Class'] = df_cleaned['Results'].map(class_mapping)

# Check for rows where 'Class' is NaN (due to missing or unexpected 'Results' values)
print(f"Rows with NaN in 'Class': {df_cleaned['Class'].isna().sum()}")

# Remove rows with NaN in the 'Class' column
df_cleaned = df_cleaned.dropna(subset=['Class'])


def one_hot_encode(sequence, max_length=None):
    """Dynamically one-hot encodes a DNA sequence, ensuring all sequences have the same length."""

    if max_length is None:
        max_length = len(sequence)  # Use the sequence's own length

    encoding = np.zeros((max_length, 4))  # Create a (max_length, 4) matrix
    base_map = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

    for i, base in enumerate(sequence[:max_length]):  # Truncate if needed
        if base in base_map:
            encoding[i, base_map[base]] = 1

    return encoding.flatten()  # Returns a 1D array

# **Determine the max length dynamically from your dataset**
max_seq_length = max(df_cleaned['Template_Sequence'].apply(len))  # Find longest sequence

# **Encode all sequences with the same max length**
encoded_sequences = np.array(
    [one_hot_encode(seq, max_length=max_seq_length) for seq in df_cleaned['Template_Sequence']])

# Normalize P90, N10, Diff
scaler = StandardScaler()
df_cleaned[['P90', 'N10', 'Diff']] = scaler.fit_transform(df_cleaned[['P90', 'N10', 'Diff']])


# **Print shape to verify**
print(f"Encoded shape: {encoded_sequences.shape}")  # Should be (num_samples, max_length * 4)
#print(encoded_sequences)

print("Feature means after scaling:")
print(df_cleaned[['P90', 'N10', 'Diff']].mean())

print("\nFeature standard deviations after scaling:")
print(df_cleaned[['P90', 'N10', 'Diff']].std())

# Prepare features and labels
X = np.hstack([encoded_sequences, df_cleaned[['P90', 'N10', 'Diff']].values])  # Combine encoded sequences with other features
y = df_cleaned['Class']

# Split into train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train Random Forest Classifier
rf_model = RandomForestClassifier(n_estimators=100)
rf_model.fit(X_train, y_train)

# Evaluate the model
print(f"Train Accuracy: {rf_model.score(X_train, y_train)}")
print(f"Test Accuracy: {rf_model.score(X_test, y_test)}")
