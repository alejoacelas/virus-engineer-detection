#!/usr/bin/env python3
"""
K-mer Frequency Baseline with Naive Bayes
Simple and fast approach using 6-mer counting
"""

import pandas as pd
import numpy as np
from collections import Counter
from sklearn.naive_bayes import MultinomialNB
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, top_k_accuracy_score

def extract_kmers(sequence, k=6):
    """Extract k-mers from DNA sequence"""
    sequence = sequence.upper().replace('N', '')  # Remove N's
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def sequence_to_kmer_counts(sequence, kmer_vocab, k=6):
    """Convert sequence to k-mer count vector"""
    kmers = extract_kmers(sequence, k)
    kmer_counts = Counter(kmers)
    return np.array([kmer_counts.get(kmer, 0) for kmer in kmer_vocab])

def train_kmer_baseline(csv_path, k=6, max_kmers=1000):
    """Train k-mer baseline model"""
    # Load data
    df = pd.read_csv(csv_path)
    
    # Extract sequences and labels
    sequences = df['sequence'].values
    labels = df['label'].values
    
    # Build k-mer vocabulary from all sequences
    all_kmers = set()
    for seq in sequences:
        all_kmers.update(extract_kmers(seq, k))
    
    # Keep most common k-mers to limit feature space
    kmer_vocab = list(all_kmers)[:max_kmers]
    print(f"Using {len(kmer_vocab)} k-mers of size {k}")
    
    # Convert sequences to k-mer count matrix
    X = np.array([sequence_to_kmer_counts(seq, kmer_vocab, k) for seq in sequences])
    
    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(labels)
    
    # Train Naive Bayes
    model = MultinomialNB(alpha=1.0)
    model.fit(X, y)
    
    # Calculate accuracy
    y_pred = model.predict(X)
    acc = accuracy_score(y, y_pred)
    top10_acc = top_k_accuracy_score(y, model.predict_proba(X), k=10)
    
    print(f"Training Accuracy: {acc:.3f}")
    print(f"Top-10 Accuracy: {top10_acc:.3f}")
    
    return model, le, kmer_vocab

if __name__ == "__main__":
    model, label_encoder, vocab = train_kmer_baseline("data/example.csv")