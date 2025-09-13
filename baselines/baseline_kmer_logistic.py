#!/usr/bin/env python3
"""
K-mer Frequency Baseline with Logistic Regression
Simple and fast approach using 6-mer counting
"""

import pandas as pd
import numpy as np
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, top_k_accuracy_score
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
from scipy.sparse import csr_matrix, hstack
import re

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def preprocess_sequence(sequence, include_reverse_complement=False):
    """Preprocess DNA sequence: uppercase, remove N's, optionally add reverse complement"""
    clean_seq = re.sub(r'[^ATGC]', '', sequence.upper())
    
    if include_reverse_complement:
        # Combine forward and reverse complement k-mers for normalization
        rc_seq = reverse_complement(clean_seq)
        return clean_seq + rc_seq  # Concatenate for k-mer extraction
    
    return clean_seq

def create_kmer_analyzer(k=6):
    """Create a k-mer analyzer function for CountVectorizer"""
    def kmer_analyzer(text):
        """Extract k-mers from preprocessed DNA sequence"""
        if len(text) < k:
            return []
        return [text[i:i+k] for i in range(len(text)-k+1)]
    return kmer_analyzer

def extract_kmers_vectorized(sequences, k=6, max_features=1000, use_tfidf=False, include_reverse_complement=False):
    """
    Extract k-mers using sklearn vectorizers (optimized version)
    
    Args:
        sequences: List of DNA sequences
        k: k-mer size
        max_features: Maximum number of k-mers to keep
        use_tfidf: Whether to use TF-IDF weighting instead of raw counts
        include_reverse_complement: Whether to include reverse complement normalization
    
    Returns:
        X: Sparse matrix of k-mer features
        kmer_vocab: List of k-mer vocabulary
        vectorizer: Fitted vectorizer for future use
    """
    # Preprocess all sequences
    processed_sequences = [preprocess_sequence(seq, include_reverse_complement) for seq in sequences]
    
    # Choose vectorizer type
    if use_tfidf:
        vectorizer = TfidfVectorizer(
            analyzer=create_kmer_analyzer(k),
            max_features=max_features,
            lowercase=False,
            token_pattern=None,  # We handle tokenization in our analyzer
            sublinear_tf=True,   # Apply sublinear TF scaling
            norm='l2'            # L2 normalization
        )
    else:
        vectorizer = CountVectorizer(
            analyzer=create_kmer_analyzer(k),
            max_features=max_features,
            lowercase=False,
            token_pattern=None  # We handle tokenization in our analyzer
        )
    
    # Fit and transform sequences to k-mer features
    X = vectorizer.fit_transform(processed_sequences)
    
    # Get the k-mer vocabulary
    kmer_vocab = vectorizer.get_feature_names_out().tolist()
    
    return X, kmer_vocab, vectorizer

def extract_multi_kmer_features(sequences, k_values=[4, 6, 8], max_features=5000, use_tfidf=True, include_reverse_complement=True):
    """
    Extract k-mer features for multiple k values and concatenate them
    
    Args:
        sequences: List of DNA sequences
        k_values: List of k-mer sizes to extract
        max_features: Maximum features per k-mer size
        use_tfidf: Whether to use TF-IDF weighting
        include_reverse_complement: Whether to normalize with reverse complement
    
    Returns:
        X_combined: Combined sparse matrix with all k-mer features
        feature_info: Dictionary with information about each k-mer size
    """
    feature_matrices = []
    feature_info = {}
    
    for k in k_values:
        print(f"Extracting {k}-mer features...")
        X_k, vocab_k, vectorizer_k = extract_kmers_vectorized(
            sequences, k=k, max_features=max_features, 
            use_tfidf=use_tfidf, include_reverse_complement=include_reverse_complement
        )
        
        feature_matrices.append(X_k)
        feature_info[f'{k}mer'] = {
            'vocab': vocab_k,
            'vectorizer': vectorizer_k,
            'num_features': X_k.shape[1]
        }
        print(f"  {k}-mer: {X_k.shape[1]} features")
    
    # Concatenate all k-mer feature matrices horizontally
    X_combined = hstack(feature_matrices)
    
    print(f"Combined feature matrix: {X_combined.shape}")
    return X_combined, feature_info

# Keep old functions for backward compatibility
def extract_kmers(sequence, k=6):
    """Extract k-mers from DNA sequence (legacy function)"""
    sequence = preprocess_sequence(sequence)
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def sequence_to_kmer_counts(sequence, kmer_vocab, k=6):
    """Convert sequence to k-mer count vector (legacy function)"""
    kmers = extract_kmers(sequence, k)
    kmer_counts = Counter(kmers)
    return np.array([kmer_counts.get(kmer, 0) for kmer in kmer_vocab])

def train_kmer_baseline(csv_path, k=6, max_kmers=1000, test_size=0.2, use_optimized=True):
    """Train k-mer baseline model with optional optimization"""
    # Load data
    df = pd.read_csv(csv_path)
    
    # Extract sequences and labels
    sequences = df['sequence'].values
    labels = df['label'].values.astype(int)
    
    if use_optimized:
        # Use optimized vectorized approach
        print(f"Building {k}-mer features using optimized vectorizer...")
        X, kmer_vocab, vectorizer = extract_kmers_vectorized(sequences, k=k, max_features=max_kmers)
        print(f"Using {len(kmer_vocab)} most common {k}-mers")
        print(f"Feature matrix shape: {X.shape} (sparse: {X.format})")
    else:
        # Use original approach for comparison
        print(f"Building {k}-mer vocabulary...")
        all_kmers = Counter()
        for seq in sequences:
            all_kmers.update(extract_kmers(seq, k))
        
        # Keep most common k-mers to limit feature space
        kmer_vocab = [kmer for kmer, count in all_kmers.most_common(max_kmers)]
        print(f"Using {len(kmer_vocab)} most common {k}-mers")
        
        # Convert sequences to k-mer count matrix
        print("Converting sequences to k-mer features...")
        X = np.array([sequence_to_kmer_counts(seq, kmer_vocab, k) for seq in sequences])
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, labels, test_size=test_size, random_state=42, stratify=labels
    )
    
    # Train Logistic Regression with class weighting for imbalanced data
    print("Training Logistic Regression model...")
    # Use sparse matrix support for LogisticRegression
    model = LogisticRegression(
        random_state=42, 
        max_iter=1000, 
        class_weight='balanced',
        solver='liblinear'  # Better for sparse matrices
    )
    model.fit(X_train, y_train)
    
    # Get probabilities (threshold calibration will be done externally)
    y_train_proba = model.predict_proba(X_train)
    y_test_proba = model.predict_proba(X_test)
    
    # Default predictions using 0.5 threshold for accuracy calculation
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    
    # Calculate accuracies
    train_acc = accuracy_score(y_train, y_train_pred)
    test_acc = accuracy_score(y_test, y_test_pred)
    
    # Top-k accuracy (k=min(10, num_classes))
    k_top = min(10, len(np.unique(labels)))
    # Get unique classes present in the data to avoid the binary/multiclass mismatch
    unique_classes = np.unique(labels)
    
    # For binary classification, use k=1 for top-k accuracy
    if len(unique_classes) == 2:
        k_top = 1
    
    try:
        train_topk = top_k_accuracy_score(y_train, y_train_proba, k=k_top, labels=unique_classes)
        test_topk = top_k_accuracy_score(y_test, y_test_proba, k=k_top, labels=unique_classes)
    except ValueError as e:
        # Fallback: if top-k fails, just use regular accuracy
        print(f"Warning: Top-k accuracy failed ({e}), using regular accuracy instead")
        train_topk = train_acc
        test_topk = test_acc
    
    print("\nK-mer Logistic Regression Results:")
    print(f"Training Accuracy: {train_acc:.3f}, Top-{k_top}: {train_topk:.3f}")
    print(f"Test Accuracy: {test_acc:.3f}, Top-{k_top}: {test_topk:.3f}")
    
    return {
        'model': model,
        'kmer_vocab': kmer_vocab,
        'y_train_true': y_train,
        'y_train_pred': y_train_pred,
        'y_test_true': y_test,
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'train_acc': train_acc,
        'test_acc': test_acc,
        'train_topk': train_topk,
        'test_topk': test_topk,
        'k': k
    }

def train_multi_kmer_baseline(csv_path, k_values=[4, 6, 8], max_kmers=5000, test_size=0.2, 
                             use_tfidf=True, include_reverse_complement=True, ensemble_method='average'):
    """
    Train multi-scale k-mer baseline with ensemble of different k values
    
    Args:
        csv_path: Path to training data CSV
        k_values: List of k-mer sizes to use
        max_kmers: Maximum k-mers per k value
        test_size: Test split fraction
        use_tfidf: Whether to use TF-IDF weighting
        include_reverse_complement: Whether to normalize with reverse complement
        ensemble_method: How to combine predictions ('average', 'weighted', 'separate')
    """
    # Load data
    df = pd.read_csv(csv_path)
    sequences = df['sequence'].values
    labels = df['label'].values.astype(int)
    
    if ensemble_method == 'combined':
        # Extract multi-k features and train single model
        print(f"Building combined multi-k-mer features: {k_values}")
        X, feature_info = extract_multi_kmer_features(
            sequences, k_values=k_values, max_features=max_kmers,
            use_tfidf=use_tfidf, include_reverse_complement=include_reverse_complement
        )
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, labels, test_size=test_size, random_state=42, stratify=labels
        )
        
        # Train single model on combined features
        print("Training combined multi-k-mer model...")
        model = LogisticRegression(
            random_state=42, max_iter=1000, class_weight='balanced',
            solver='liblinear'
        )
        model.fit(X_train, y_train)
        
        # Get predictions
        y_train_proba = model.predict_proba(X_train)
        y_test_proba = model.predict_proba(X_test)
        y_train_pred = model.predict(X_train)
        y_test_pred = model.predict(X_test)
        
        # Calculate metrics
        train_acc = accuracy_score(y_train, y_train_pred)
        test_acc = accuracy_score(y_test, y_test_pred)
        
        # Top-k accuracy
        k_top = min(10, len(np.unique(labels)))
        unique_classes = np.unique(labels)
        
        if len(unique_classes) == 2:
            k_top = 1
        
        try:
            train_topk = top_k_accuracy_score(y_train, y_train_proba, k=k_top, labels=unique_classes)
            test_topk = top_k_accuracy_score(y_test, y_test_proba, k=k_top, labels=unique_classes)
        except ValueError as e:
            print(f"Warning: Top-k accuracy failed ({e}), using regular accuracy instead")
            train_topk = train_acc
            test_topk = test_acc
        
        print(f"\nCombined Multi-K-mer Results:")
        print(f"K-values: {k_values}")
        print(f"Total features: {X.shape[1]}")
        for k_val in k_values:
            print(f"  {k_val}-mer: {feature_info[f'{k_val}mer']['num_features']} features")
        print(f"Training Accuracy: {train_acc:.3f}, Top-{k_top}: {train_topk:.3f}")
        print(f"Test Accuracy: {test_acc:.3f}, Top-{k_top}: {test_topk:.3f}")
        
        return {
            'model': model,
            'feature_info': feature_info,
            'y_train_true': y_train,
            'y_train_pred': y_train_pred,
            'y_test_true': y_test,
            'y_test_pred': y_test_pred,
            'y_train_proba': y_train_proba,
            'y_test_proba': y_test_proba,
            'train_acc': train_acc,
            'test_acc': test_acc,
            'train_topk': train_topk,
            'test_topk': test_topk,
            'k_values': k_values,
            'ensemble_method': ensemble_method
        }
    
    else:
        # Train separate models for each k and ensemble predictions
        print(f"Training separate models for k-values: {k_values}")
        models = {}
        results = {}
        
        # Train individual models
        for k in k_values:
            print(f"\n--- Training {k}-mer model ---")
            result = train_kmer_baseline(
                csv_path, k=k, max_kmers=max_kmers, test_size=test_size, use_optimized=True
            )
            models[k] = result['model']
            results[k] = result
        
        # Get predictions from all models
        train_probas = []
        test_probas = []
        
        for k in k_values:
            train_probas.append(results[k]['y_train_proba'])
            test_probas.append(results[k]['y_test_proba'])
        
        # Ensemble predictions
        if ensemble_method == 'average':
            # Simple average
            y_train_proba_ensemble = np.mean(train_probas, axis=0)
            y_test_proba_ensemble = np.mean(test_probas, axis=0)
        elif ensemble_method == 'weighted':
            # Weight by individual model accuracy
            weights = [results[k]['test_acc'] for k in k_values]
            weights = np.array(weights) / np.sum(weights)
            
            y_train_proba_ensemble = np.average(train_probas, axis=0, weights=weights)
            y_test_proba_ensemble = np.average(test_probas, axis=0, weights=weights)
        else:
            raise ValueError(f"Unknown ensemble method: {ensemble_method}")
        
        # Get ensemble predictions
        y_train_pred_ensemble = np.argmax(y_train_proba_ensemble, axis=1)
        y_test_pred_ensemble = np.argmax(y_test_proba_ensemble, axis=1)
        
        # Calculate ensemble metrics
        y_train = results[k_values[0]]['y_train_true']  # Same for all models
        y_test = results[k_values[0]]['y_test_true']
        
        train_acc = accuracy_score(y_train, y_train_pred_ensemble)
        test_acc = accuracy_score(y_test, y_test_pred_ensemble)
        
        # Top-k accuracy
        k_top = min(10, len(np.unique(labels)))
        unique_classes = np.unique(labels)
        
        if len(unique_classes) == 2:
            k_top = 1
        
        try:
            train_topk = top_k_accuracy_score(y_train, y_train_proba_ensemble, k=k_top, labels=unique_classes)
            test_topk = top_k_accuracy_score(y_test, y_test_proba_ensemble, k=k_top, labels=unique_classes)
        except ValueError as e:
            print(f"Warning: Top-k accuracy failed ({e}), using regular accuracy instead")
            train_topk = train_acc
            test_topk = test_acc
        
        print(f"\n=== ENSEMBLE RESULTS ===")
        print(f"Method: {ensemble_method}")
        print(f"K-values: {k_values}")
        
        # Show individual model performance
        for k in k_values:
            print(f"  {k}-mer individual: Train {results[k]['train_acc']:.3f}, Test {results[k]['test_acc']:.3f}")
        
        print(f"Ensemble: Train {train_acc:.3f}, Test {test_acc:.3f}")
        print(f"Ensemble Top-{k_top}: Train {train_topk:.3f}, Test {test_topk:.3f}")
        
        return {
            'models': models,
            'individual_results': results,
            'y_train_true': y_train,
            'y_train_pred': y_train_pred_ensemble,
            'y_test_true': y_test,
            'y_test_pred': y_test_pred_ensemble,
            'y_train_proba': y_train_proba_ensemble,
            'y_test_proba': y_test_proba_ensemble,
            'train_acc': train_acc,
            'test_acc': test_acc,
            'train_topk': train_topk,
            'test_topk': test_topk,
            'k_values': k_values,
            'ensemble_method': ensemble_method
        }

if __name__ == "__main__":
    results = train_kmer_baseline("data/dummy_train.csv")