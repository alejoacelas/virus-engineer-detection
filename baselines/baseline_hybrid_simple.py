#!/usr/bin/env python3
"""
Hybrid Simple Baseline
Combines k-mer features with metadata for ensemble prediction
"""

import pandas as pd
import numpy as np
from collections import Counter
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, top_k_accuracy_score

def extract_kmers(sequence, k=6):
    """Extract k-mers from DNA sequence"""
    sequence = sequence.upper().replace('N', '')
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def get_kmer_features(sequences, max_kmers=500):
    """Get k-mer count features for all sequences"""
    # Build vocabulary
    all_kmers = set()
    for seq in sequences:
        all_kmers.update(extract_kmers(seq))
    
    # Take most frequent k-mers
    kmer_vocab = list(all_kmers)[:max_kmers]
    
    # Convert to count matrix
    kmer_matrix = []
    for seq in sequences:
        kmer_counts = Counter(extract_kmers(seq))
        features = [kmer_counts.get(kmer, 0) for kmer in kmer_vocab]
        kmer_matrix.append(features)
    
    return np.array(kmer_matrix), kmer_vocab

def train_hybrid_baseline(csv_path):
    """Train hybrid baseline combining k-mers and metadata"""
    # Load data
    df = pd.read_csv(csv_path)
    
    sequences = df['sequence'].values
    labels = df['label'].values
    
    # Get metadata features
    feature_cols = [col for col in df.columns 
                   if col not in ['sequence_id', 'sequence', 'label']]
    metadata_features = df[feature_cols].values
    
    # Get k-mer features
    print("Extracting k-mer features...")
    kmer_features, kmer_vocab = get_kmer_features(sequences)
    
    print(f"K-mer features: {kmer_features.shape[1]}")
    print(f"Metadata features: {metadata_features.shape[1]}")
    
    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(labels)
    
    # Split data
    indices = np.arange(len(sequences))
    train_idx, val_idx = train_test_split(indices, test_size=0.2, random_state=42, stratify=y)
    
    # Train k-mer model
    kmer_model = MultinomialNB(alpha=1.0)
    kmer_model.fit(kmer_features[train_idx], y[train_idx])
    
    # Train metadata model
    scaler = StandardScaler()
    metadata_scaled = scaler.fit_transform(metadata_features)
    
    meta_model = LogisticRegression(max_iter=1000, class_weight='balanced')
    meta_model.fit(metadata_scaled[train_idx], y[train_idx])
    
    # Get predictions from both models
    kmer_train_proba = kmer_model.predict_proba(kmer_features[train_idx])
    kmer_val_proba = kmer_model.predict_proba(kmer_features[val_idx])
    
    meta_train_proba = meta_model.predict_proba(metadata_scaled[train_idx])
    meta_val_proba = meta_model.predict_proba(metadata_scaled[val_idx])
    
    # Simple ensemble: average probabilities
    ensemble_train_proba = (kmer_train_proba + meta_train_proba) / 2
    ensemble_val_proba = (kmer_val_proba + meta_val_proba) / 2
    
    # Predictions
    ensemble_train_pred = np.argmax(ensemble_train_proba, axis=1)
    ensemble_val_pred = np.argmax(ensemble_val_proba, axis=1)
    
    # Calculate metrics
    train_acc = accuracy_score(y[train_idx], ensemble_train_pred)
    val_acc = accuracy_score(y[val_idx], ensemble_val_pred)
    
    train_top10 = top_k_accuracy_score(y[train_idx], ensemble_train_proba, k=10)
    val_top10 = top_k_accuracy_score(y[val_idx], ensemble_val_proba, k=10)
    
    print(f"\nHybrid Ensemble Results:")
    print(f"Training Accuracy: {train_acc:.3f}, Top-10: {train_top10:.3f}")
    print(f"Validation Accuracy: {val_acc:.3f}, Top-10: {val_top10:.3f}")
    
    # Individual model performance for comparison
    kmer_train_acc = accuracy_score(y[train_idx], kmer_model.predict(kmer_features[train_idx]))
    meta_train_acc = accuracy_score(y[train_idx], meta_model.predict(metadata_scaled[train_idx]))
    
    print(f"\nIndividual Models (Training):")
    print(f"K-mer only: {kmer_train_acc:.3f}")
    print(f"Metadata only: {meta_train_acc:.3f}")
    
    models = {
        'kmer_model': kmer_model,
        'meta_model': meta_model,
        'scaler': scaler,
        'label_encoder': le,
        'kmer_vocab': kmer_vocab,
        'feature_cols': feature_cols
    }
    
    return models

if __name__ == "__main__":
    models = train_hybrid_baseline("data/example.csv")