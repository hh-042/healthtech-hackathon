#PWM/SVM only

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt

# Loading data in
df = pd.read_csv('/Users/angelyao/PycharmProjects/healthHackathon/healthAIHack/384_sequences.csv')

# Remove rows with invalid results (e.g., NR, CE, SD, etc.)
df_cleaned = df[~df['Results'].isin(['NR', 'CE', 'SD', 'LA', 'NA', 'NS'])]

# Remove rows with zero values in P90, N10, or Diff
df_cleaned = df_cleaned[(df_cleaned[['P90', 'N10', 'Diff']] != 0).all(axis=1)]

# Class labels (CI = good, CII = bad)
df_cleaned['Class'] = df_cleaned['Results'].apply(lambda x: 1 if x == 'CI' else 0)

# real pwm from etseq github https://github.com/expartools/ETSeq/blob/master/Source_code_0.5.2.zip
pwm_p90 = np.array([
    [0.8473,0.9985,1.1451,1.0116,0.6419,0.5465,0.7376,0.6506,1.9924,
     1.0986,0.0000,0.0606,0.3483,0.0000,0.8473,0.9985,1.1451,1.0116,
     0.6419,0.5465,0.7376,0.6506,1.9924,1.0986],
    [0.2719, 0.0000, 0.2007, 0.6061, 0.3054, 0.2364, 0.8329, 0.4274,
     0.3830, 0.0000, 0.1178, 0.0000, 0.0000, 0.1335, 0.2719, 0.0000,
     0.2007, 0.6061, 0.3054, 0.2364, 0.8329, 0.4274, 0.3830, 0.0000],
    [0.3365, 0.0541, 0.0000, 0.1466, 0.0541, 0.1719, 0.3023, 0.7376,
     0.0465, 0.5390, 0.1178, 0.3483, 0.0606, 0.0645, 0.3365, 0.0541,
     0.0000, 0.1466, 0.0541, 0.1719, 0.3023, 0.7376, 0.0465, 0.5390],
    [0.0000, 0.1112, 0.4520, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     0.0000, 0.4700, 0.4925, 0.0606, 0.0606, 0.0000, 0.0000, 0.1112,
     0.4520, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4700]
    ])

pwm_diff = np.array([
    [0.470003629245736, 0.628608659422374, 0.191055236762709, 0.969400557188104,
     0.887303195000903, 0.594707107746693, 1.12846525181779, 1.09861228866811,
     0.810930216216329, 1.09861228866811, 0, 0, 0.550046336919272, 0, 0.470003629245736,
     0.628608659422374, 0.191055236762709, 0.969400557188104, 0.887303195000903, 0.594707107746693,
     1.12846525181779, 1.09861228866811, 0.810930216216329, 1.09861228866811],
    [0.287682072451781, 0.693147180559945, 0, 0.231801614057324, 1.04145387482816, 0.371563556432483,
     0.961411167154625, 0.587786664902119, 0.204794412646013, 0, 0.559615787935423, 0.0953101798043249,
     0.122602322092332, 0.424883193965266, 0.287682072451781, 0.693147180559945, 0, 0.231801614057324,
     1.04145387482816, 0.371563556432483, 0.961411167154625, 0.587786664902119, 0.204794412646013, 0],
    [0, 0.405465108108164, 0, 0.476924072090309, 0.481838086892738, 0.594707107746693,
     0.390866308687012, 1.01856958099457, 0.300104592450338, 0.916290731874155, 0.287682072451781,
     0.146603474191875, 0.424883193965266, 0.367724780125317, 0, 0.405465108108164, 0,
     0.476924072090309, 0.481838086892738, 0.594707107746693, 0.390866308687012, 1.01856958099457,
     0.300104592450338, 0.916290731874155],
    [0, 0, 0.362905493689368, 0, 0, 0, 0, 0, 0, 0.0339015516756814, 0.559615787935423,
     0.0953101798043249, 0, 0.262364264467491, 0, 0, 0.362905493689368, 0, 0, 0, 0, 0, 0,
     0.0339015516756814]
])

# Feature selection (P90 and Diff as input for classification)
X = df_cleaned[['P90', 'Diff']].values
y = df_cleaned['Class'].values

# Standardizing features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split data into training and testing sets (80% training, 20% testing)
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Train the SVM classifier on the training data
svm_classifier = SVC(kernel='linear')
svm_classifier.fit(X_train, y_train)

# Make predictions on the testing data
y_pred = svm_classifier.predict(X_test)

# Compute the accuracy
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy * 100:.2f}%")

def compute_features(sequence):
    # Truncate or pad the sequence to 15 bases
    sequence = sequence[:22]  # Take the first 15 bases
    sequence = sequence.ljust(22, 'A')  # Pad with 'N' if it's shorter than 15

    seq = Seq(sequence)

    # Calculate melting temperature as a feature
    melting_temp = mt.Tm_NN(seq)

    # Compute P90, N10, and Diff scores using PWM
    def compute_pwm_score(sequence, pwm):
        score = 0
        for i, nucleotide in enumerate(sequence):
            if nucleotide == 'A':
                score += pwm[0][i]
            elif nucleotide == 'T':
                score += pwm[1][i]
            elif nucleotide == 'G':
                score += pwm[2][i]
            elif nucleotide == 'C':
                score += pwm[3][i]
        return score / len(sequence)  # Normalize score by length

    p90 = compute_pwm_score(sequence, pwm_p90)
    diff = compute_pwm_score(sequence, pwm_diff)

    return p90, diff, melting_temp

def classify_sequence(sequence):
    p90, diff, melting_temp = compute_features(sequence)

    print(f"p90: {p90}, diff: {diff}") #add this line

    # Transform input
    input_scaled = scaler.transform([[p90, diff]])
    print(f"input_scaled: {input_scaled}") #add this line

    # Predict using SVM
    prediction = svm_classifier.predict(input_scaled)
    return 'CI' if prediction[0] == 1 else 'CII'

# Example: Classify new sequence
new_sequence = "TGCGAGTGGAGCGGGACTCTTGCGAGTGGA"  # Replace with any sequence
print(classify_sequence(new_sequence))  # Output: CI or CII

# Create a grid of points to plot the decision boundary
xx, yy = np.meshgrid(np.linspace(-1, 2, 100), np.linspace(-1, 2, 100))

# Predict the SVM output for each point on the grid
Z = svm_classifier.predict(np.c_[xx.ravel(), yy.ravel()])

# Reshape the prediction results to match the grid
Z = Z.reshape(xx.shape)

# Plot the decision boundary
plt.contourf(xx, yy, Z, alpha=0.75)
plt.scatter([x[0] for x in X], [x[1] for x in X], c=y, edgecolors='k', marker='o', cmap='coolwarm')
plt.xlabel('P90 Score')
plt.ylabel('Diff Score')
plt.title('SVM Decision Boundary')
plt.show()
