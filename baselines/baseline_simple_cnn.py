#!/usr/bin/env python3
"""
Simple CNN Baseline
Basic convolutional neural network for DNA sequence classification
"""

import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

def dna_to_onehot(sequence, max_len=5000):
    """Convert DNA sequence to one-hot encoding"""
    # Map nucleotides to integers
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    # Truncate or pad sequence
    sequence = sequence.upper()[:max_len]
    seq_encoded = [mapping.get(nuc, 4) for nuc in sequence]
    
    # Pad with N's if needed
    if len(seq_encoded) < max_len:
        seq_encoded.extend([4] * (max_len - len(seq_encoded)))
    
    # Convert to one-hot (exclude N channel for simplicity)
    onehot = np.zeros((max_len, 4))
    for i, nuc_idx in enumerate(seq_encoded):
        if nuc_idx < 4:  # Skip N's
            onehot[i, nuc_idx] = 1
    
    return onehot

def create_simple_cnn(input_shape, num_classes):
    """Create simple CNN architecture"""
    model = tf.keras.Sequential([
        tf.keras.layers.Conv1D(64, 15, activation='relu', input_shape=input_shape),
        tf.keras.layers.GlobalMaxPooling1D(),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dropout(0.3),
        tf.keras.layers.Dense(num_classes, activation='softmax')
    ])
    
    return model

def train_simple_cnn(csv_path, max_len=5000, epochs=10):
    """Train simple CNN baseline"""
    # Load data
    df = pd.read_csv(csv_path)
    
    # Extract sequences and labels
    sequences = df['sequence'].values
    labels = df['label'].values
    
    # Convert to one-hot encoding
    print("Converting sequences to one-hot encoding...")
    X = np.array([dna_to_onehot(seq, max_len) for seq in sequences])
    
    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(labels)
    num_classes = len(le.classes_)
    
    # Convert to categorical
    y_categorical = tf.keras.utils.to_categorical(y, num_classes)
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y_categorical, test_size=0.2, random_state=42
    )
    
    # Create and compile model
    model = create_simple_cnn((max_len, 4), num_classes)
    model.compile(
        optimizer='adam',
        loss='categorical_crossentropy',
        metrics=['accuracy', tf.keras.metrics.TopKCategoricalAccuracy(k=10, name='top10_acc')]
    )
    
    print(f"Training on {len(X_train)} sequences, {num_classes} classes")
    
    # Train model
    history = model.fit(
        X_train, y_train,
        batch_size=32,
        epochs=epochs,
        validation_data=(X_val, y_val),
        verbose=1
    )
    
    # Evaluate
    train_loss, train_acc, train_top10 = model.evaluate(X_train, y_train, verbose=0)
    val_loss, val_acc, val_top10 = model.evaluate(X_val, y_val, verbose=0)
    
    print(f"Training Accuracy: {train_acc:.3f}, Top-10: {train_top10:.3f}")
    print(f"Validation Accuracy: {val_acc:.3f}, Top-10: {val_top10:.3f}")
    
    return model, le

if __name__ == "__main__":
    model, label_encoder = train_simple_cnn("data/example.csv")