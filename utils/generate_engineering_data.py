#!/usr/bin/env python3
"""
Engineering Detection Data Generator
Creates viral genome samples with some sequences containing engineered segments
"""

import pandas as pd
import numpy as np
from pathlib import Path

def generate_random_segment_samples(sequences, n_samples, mean_length=187, sd_length=50):
    """
    Generate random segments from viral sequences with variable lengths.
    Optimized version using vectorized operations.
    
    Args:
        sequences: List of viral genome sequences
        n_samples: Total number of samples to generate
        mean_length: Mean length of generated segments
        sd_length: Standard deviation of segment lengths
    
    Returns:
        List of sequence segments
    """
    # Pre-allocate array for samples
    samples = [None] * n_samples
    
    # Pre-generate all random choices
    sequence_indices = np.random.randint(0, len(sequences), n_samples)
    segment_lengths = np.random.normal(mean_length, sd_length, n_samples).astype(int)
    
    # Process sequences in batches for better memory efficiency
    batch_size = min(10000, n_samples)
    
    for batch_start in range(0, n_samples, batch_size):
        batch_end = min(batch_start + batch_size, n_samples)
        batch_indices = sequence_indices[batch_start:batch_end]
        batch_lengths = segment_lengths[batch_start:batch_end]
        
        for i, (seq_idx, segment_length) in enumerate(zip(batch_indices, batch_lengths)):
            seq = sequences[seq_idx]
            
            # Ensure segment length is valid
            segment_length = max(1, min(len(seq), segment_length))
            
            # Choose random start position
            max_start = max(0, len(seq) - segment_length)
            start_pos = np.random.randint(0, max_start + 1)
            
            # Extract segment
            samples[batch_start + i] = seq[start_pos:start_pos + segment_length]
    
    return samples

def engineer_sequence(sequence, max_replacement_fraction=2/3):
    """
    Replace part of a sequence with random nucleotides.
    Optimized version using bulk random generation.
    
    Args:
        sequence: Original sequence
        max_replacement_fraction: Maximum fraction of sequence to replace
    
    Returns:
        Engineered sequence with same length
    """
    nucleotides = np.array(['A', 'T', 'G', 'C'])
    # nucleotides = np.array(['A', 'A', 'A', 'A'])
    seq_len = len(sequence)
    
    # Choose random replacement length (up to max_replacement_fraction)
    max_replace_len = int(seq_len * max_replacement_fraction)
    if max_replace_len <= 0:
        max_replace_len = 1
    replace_len = np.random.randint(1, max_replace_len + 1)
    
    # Choose random start position for replacement
    start_pos = np.random.randint(0, seq_len - replace_len + 1)
    
    # Generate random replacement using vectorized operations
    replacement = ''.join(np.random.choice(nucleotides, replace_len))
    
    # Create engineered sequence
    engineered = sequence[:start_pos] + replacement + sequence[start_pos + replace_len:]
    
    return engineered

def load_engineering_detection_data(csv_path, n_samples=1000, mean_length=187, sd_length=50, 
                                  engineering_rate=0.02, min_sequence_length=500):
    """
    Load viral genome data and create samples for engineering detection.
    
    Args:
        csv_path: Path to viral_genome_df_prelim.csv
        n_samples: Total number of samples to generate
        mean_length: Mean length of generated segments
        sd_length: Standard deviation of segment lengths
        engineering_rate: Fraction of samples to engineer (default 2%)
        min_sequence_length: Minimum length for source sequences
    
    Returns:
        DataFrame with sequence_id, sequence, label (0=original, 1=engineered)
    """
    print(f"Loading viral genome data from {csv_path}...")
    
    # Load the viral genome data
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} viral sequences")
    
    # Filter by length and remove sequences with too many N's
    df = df[df['genome_sequence'].str.len() >= min_sequence_length].copy()
    df = df[df['genome_sequence'].str.count('N') / df['genome_sequence'].str.len() < 0.1].copy()
    
    print(f"After filtering: {len(df)} sequences available for sampling")
    
    if len(df) == 0:
        raise ValueError("No sequences available after filtering")
    
    # Generate random segments
    print(f"Generating {n_samples} random segments (mean length: {mean_length}, sd: {sd_length})...")
    sequences = df['genome_sequence'].tolist()
    samples = generate_random_segment_samples(sequences, n_samples, mean_length, sd_length)
    
    # Determine which samples to engineer
    n_engineered = int(n_samples * engineering_rate)
    engineered_indices = np.random.choice(n_samples, n_engineered, replace=False)
    
    print(f"Engineering {n_engineered} sequences ({engineering_rate*100:.1f}%)...")
    
    # Pre-allocate arrays for better performance
    labels = np.zeros(n_samples, dtype=int)
    final_sequences = [None] * n_samples
    
    # Set labels for engineered sequences
    labels[engineered_indices] = 1
    
    # Process sequences in batches for memory efficiency
    batch_size = min(10000, n_samples)
    
    for batch_start in range(0, n_samples, batch_size):
        batch_end = min(batch_start + batch_size, n_samples)
        batch_indices = np.arange(batch_start, batch_end)
        
        for i in batch_indices:
            if i in engineered_indices:
                # Engineer this sequence
                final_sequences[i] = engineer_sequence(samples[i])
            else:
                # Keep original
                final_sequences[i] = samples[i]
    
    # Create final dataset using vectorized operations
    sequence_ids = [f"SEQ_{i:05d}" for i in range(n_samples)]
    result_df = pd.DataFrame({
        'sequence_id': sequence_ids,
        'sequence': final_sequences,
        'label': labels
    })
    
    print(f"Created dataset with {n_samples} samples")
    print(f"Label distribution: {pd.Series(labels).value_counts().sort_index().to_dict()}")
    # Calculate average sequence length more efficiently
    avg_length = np.mean([len(s) for s in final_sequences])
    print(f"Average sequence length: {avg_length:.1f}")
    
    return result_df

def save_engineering_detection_data(output_path="data/engineering_detection_train.csv", skip_if_exists=False, **kwargs):
    """
    Save engineering detection data to CSV file
    
    Args:
        output_path: Path to save the CSV file
        skip_if_exists: If True, skip generation if output file already exists
        **kwargs: Arguments passed to load_engineering_detection_data
    
    Returns:
        DataFrame with the data (either newly generated or loaded from existing file)
    """
    output_path = Path(output_path)
    
    # Check if file exists and skip_if_exists is True
    if skip_if_exists and output_path.exists():
        print(f"File {output_path} already exists, skipping generation")
        # Load and return the existing data
        detection_data = pd.read_csv(output_path)
        return detection_data
    
    # Generate new data
    data = load_engineering_detection_data(**kwargs)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    data.to_csv(output_path, index=False)
    print(f"Saved engineering detection data to {output_path}")
    
    return data

if __name__ == "__main__":
    # Example usage
    detection_data = save_engineering_detection_data(
        csv_path="data/viral_genome_df_prelim.csv",
        output_path="data/engineering_detection_train.csv",
        n_samples=2000,
        mean_length=187,
        sd_length=50
    )