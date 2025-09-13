#!/usr/bin/env python3
"""
K-mer Logistic Regression Baseline
Standalone implementation for DNA sequence classification
"""

import numpy as np
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from scipy.sparse import hstack
import re

def preprocess_sequence(sequence):
    """Clean DNA sequence - uppercase, remove N's"""
    return re.sub(r'[^ATGC]', '', sequence.upper())

def create_kmer_analyzer(k):
    """Create k-mer analyzer for vectorizer"""
    def analyzer(text):
        if len(text) < k:
            return []
        return [text[i:i+k] for i in range(len(text)-k+1)]
    return analyzer

def extract_kmer_features(sequences, k=6, max_features=1000, use_tfidf=False):
    """Extract k-mer features using sklearn vectorizers"""
    processed_seqs = [preprocess_sequence(seq) for seq in sequences]

    if use_tfidf:
        vectorizer = TfidfVectorizer(
            analyzer=create_kmer_analyzer(k),
            max_features=max_features,
            lowercase=False,
            token_pattern=None
        )
    else:
        vectorizer = CountVectorizer(
            analyzer=create_kmer_analyzer(k),
            max_features=max_features,
            lowercase=False,
            token_pattern=None
        )

    X = vectorizer.fit_transform(processed_seqs)
    return X, vectorizer

def extract_multi_kmer_features(sequences, k_values=[4, 6, 8], max_features=500, use_tfidf=True):
    """Extract and combine multiple k-mer sizes"""
    feature_matrices = []

    for k in k_values:
        X_k, _ = extract_kmer_features(sequences, k, max_features, use_tfidf)
        feature_matrices.append(X_k)

    return hstack(feature_matrices)

def train_kmer_baseline(X_train, X_test, y_train, y_test, k=6, max_features=1000, use_tfidf=False):
    """Train k-mer logistic regression baseline"""
    # Extract features from sequences
    train_sequences = X_train['sequence']
    test_sequences = X_test['sequence']

    X_train_features, vectorizer = extract_kmer_features(
        train_sequences, k, max_features, use_tfidf
    )
    X_test_features = vectorizer.transform([preprocess_sequence(seq) for seq in test_sequences])

    # Train model
    model = LogisticRegression(
        random_state=42,
        max_iter=1000,
        class_weight='balanced',
        solver='liblinear'
    )
    model.fit(X_train_features, y_train)

    # Predictions
    y_train_pred = model.predict(X_train_features)
    y_test_pred = model.predict(X_test_features)
    y_train_proba = model.predict_proba(X_train_features)
    y_test_proba = model.predict_proba(X_test_features)

    return {
        'model': model,
        'vectorizer': vectorizer,
        'y_train_pred': y_train_pred,
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'params': {'k': k, 'max_features': max_features, 'use_tfidf': use_tfidf}
    }

def train_multi_kmer_baseline(X_train, X_test, y_train, y_test, k_values=[4, 6, 8],
                             max_features=500, use_tfidf=True):
    """Train multi-k-mer baseline with combined features"""
    # Extract features
    train_sequences = X_train['sequence']
    test_sequences = X_test['sequence']

    X_train_features = extract_multi_kmer_features(
        train_sequences, k_values, max_features, use_tfidf
    )
    X_test_features = extract_multi_kmer_features(
        test_sequences, k_values, max_features, use_tfidf
    )

    # Train model
    model = LogisticRegression(
        random_state=42,
        max_iter=1000,
        class_weight='balanced',
        solver='liblinear'
    )
    model.fit(X_train_features, y_train)

    # Predictions
    y_train_pred = model.predict(X_train_features)
    y_test_pred = model.predict(X_test_features)
    y_train_proba = model.predict_proba(X_train_features)
    y_test_proba = model.predict_proba(X_test_features)

    return {
        'model': model,
        'y_train_pred': y_train_pred,
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'params': {'k_values': k_values, 'max_features': max_features, 'use_tfidf': use_tfidf}
    }