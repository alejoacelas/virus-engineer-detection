#!/usr/bin/env python3
"""
Enhanced Data Generator for Virus Engineering Detection
Creates datasets with configurable engineering fractions and categorical features
"""

import pandas as pd
import numpy as np
from pathlib import Path

def create_segments(sequences, n_samples, mean_length=187, sd_length=50):
    """Generate random segments from viral sequences (optimized)"""
    # Pre-allocate array for samples
    samples = [None] * n_samples

    # Pre-generate all random choices
    sequence_indices = np.random.randint(0, len(sequences), n_samples)
    segment_lengths = np.abs(np.random.normal(mean_length, sd_length, n_samples)).astype(int)
    segment_lengths = np.clip(segment_lengths, 10, None)  # Minimum length 10

    # Process in batches for memory efficiency
    batch_size = min(10000, n_samples)

    for batch_start in range(0, n_samples, batch_size):
        batch_end = min(batch_start + batch_size, n_samples)
        batch_indices = sequence_indices[batch_start:batch_end]
        batch_lengths = segment_lengths[batch_start:batch_end]

        for i, (seq_idx, length) in enumerate(zip(batch_indices, batch_lengths)):
            seq = sequences[seq_idx]
            length = max(1, min(len(seq), length))
            max_start = max(0, len(seq) - length)
            start = np.random.randint(0, max_start + 1) if max_start > 0 else 0
            samples[batch_start + i] = seq[start:start + length]

    return samples

def engineer_sequence(sequence, method='random_replacement'):
    """Engineer a sequence using specified method"""
    if method == 'random_replacement':
        nucleotides = np.array(['A', 'T', 'G', 'C'])
        seq_len = len(sequence)
        max_replace = int(seq_len * 0.66)
        min_replace = int(seq_len * 0.10)
        if max_replace <= min_replace:
            max_replace = min_replace + 1
        replace_len = np.random.randint(min_replace, max_replace + 1)
        start_pos = np.random.randint(0, seq_len - replace_len + 1)
        replacement = ''.join(np.random.choice(nucleotides, replace_len))
        return sequence[:start_pos] + replacement + sequence[start_pos + replace_len:]

    elif method == 'codon_optimization':
        # Simple codon optimization simulation
        codons = {'TTT': 'TTC', 'TTA': 'CTG', 'TCT': 'AGC', 'TAT': 'TAC'}
        result = sequence
        for old_codon, new_codon in codons.items():
            result = result.replace(old_codon, new_codon)
        return result

    elif method == 'motif_insertion':
        # Insert artificial restriction sites
        motifs = ['GAATTC', 'AAGCTT', 'GGATCC']
        motif = np.random.choice(motifs)
        insert_pos = np.random.randint(0, len(sequence) + 1)
        return sequence[:insert_pos] + motif + sequence[insert_pos:]

    return sequence

def create_engineering_dataset(csv_path, n_samples=1000, engineering_fraction=0.02,
                             methods=['random_replacement'], mean_length=187, sd_length=50):
    """Create dataset with engineered sequences and categorical features"""
    # Load and filter viral sequences
    df = pd.read_csv(csv_path)
    sequences = [seq for seq in df['genome_sequence']
                if len(seq) >= 500 and seq.count('N') / len(seq) < 0.1]

    if len(sequences) == 0:
        raise ValueError("No valid sequences found")

    # Generate segments
    segments = create_segments(sequences, n_samples, mean_length, sd_length)

    # Determine which to engineer
    n_engineered = int(n_samples * engineering_fraction)
    engineered_indices = set(np.random.choice(n_samples, n_engineered, replace=False))

    # Pre-allocate arrays for better performance
    final_sequences = [None] * n_samples
    labels = np.zeros(n_samples, dtype=int)
    engineering_methods = [None] * n_samples
    sequence_lengths = np.zeros(n_samples, dtype=int)

    # Process in batches for memory efficiency
    batch_size = min(10000, n_samples)

    for batch_start in range(0, n_samples, batch_size):
        batch_end = min(batch_start + batch_size, n_samples)

        for i in range(batch_start, batch_end):
            if i in engineered_indices:
                method = np.random.choice(methods)
                final_sequences[i] = engineer_sequence(segments[i], method)
                labels[i] = 1
                engineering_methods[i] = method
            else:
                final_sequences[i] = segments[i]
                engineering_methods[i] = 'natural'

            sequence_lengths[i] = len(final_sequences[i])

    # Create length bins
    length_bins = pd.cut(sequence_lengths, bins=3, labels=['short', 'medium', 'long'])

    return pd.DataFrame({
        'sequence_id': [f"SEQ_{i:05d}" for i in range(n_samples)],
        'sequence': final_sequences,
        'label': labels,
        'engineering_method': engineering_methods,
        'sequence_length': sequence_lengths,
        'length_bin': length_bins
    })

def split_dataset(df, test_size=0.2, stratify_cols=['label']):
    """Split dataset maintaining stratification"""
    from sklearn.model_selection import train_test_split

    if len(stratify_cols) == 1:
        stratify = df[stratify_cols[0]]
    else:
        stratify = df[stratify_cols].apply(lambda x: '_'.join(x.astype(str)), axis=1)

    train_df, test_df = train_test_split(
        df, test_size=test_size, random_state=42, stratify=stratify
    )

    return train_df.reset_index(drop=True), test_df.reset_index(drop=True)