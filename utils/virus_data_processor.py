#!/usr/bin/env python3
"""
Virus Data Processor
Processes engineered virus JSON data and original virus FASTA files into CSV datasets
with a similar interface to the existing data_generator.py
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import re

def load_engineered_data(json_path: str) -> pd.DataFrame:
    """
    Load and process engineered virus data from JSON file

    Args:
        json_path: Path to the restriction_virus_engineered.json file

    Returns:
        DataFrame with processed engineered virus data
    """
    print(f"Loading engineered data from {json_path}...")

    with open(json_path, 'r') as f:
        data = json.load(f)

    print(f"Loaded {len(data)} engineered entries")

    # Extract relevant fields
    processed_data = []
    for entry in data:
        # Handle cases where start/end might be missing for some engineering methods
        start_pos = entry.get('start', -1)  # -1 indicates no specific region defined
        end_pos = entry.get('end', -1)

        processed_data.append({
            'sequence': entry['full_genome'],
            'virus': entry['virus'],
            'virus_key': entry['virus_key'],
            'virus_family': entry['virus_family'],
            'example_id': entry['example_id'],
            'engineering_method': entry['method'],
            'sequence_length': len(entry['full_genome']),
            'start': start_pos,
            'end': end_pos,
            'label': 1  # Engineered
        })

    return pd.DataFrame(processed_data)

def load_original_virus_data(data_dir: str) -> pd.DataFrame:
    """
    Load and process original virus data from FASTA files

    Args:
        data_dir: Path to directory containing FASTA files

    Returns:
        DataFrame with processed original virus data
    """
    print(f"Loading original virus data from {data_dir}...")

    data_path = Path(data_dir)
    fasta_files = list(data_path.glob("*.fasta"))

    print(f"Found {len(fasta_files)} FASTA files")

    processed_data = []

    # Create virus family mapping based on common patterns
    family_mapping = create_virus_family_mapping()

    for fasta_file in fasta_files:
        # Extract virus name from filename
        filename = fasta_file.stem
        virus_name = extract_virus_name_from_filename(filename)
        virus_key = create_virus_key(virus_name)

        # Read FASTA file
        sequence = read_fasta_sequence(fasta_file)

        if sequence:
            # Determine virus family from filename patterns
            virus_family = determine_virus_family(filename, family_mapping)

            processed_data.append({
                'sequence': sequence,
                'virus': virus_name,
                'virus_key': virus_key,
                'virus_family': virus_family,
                'example_id': f"{virus_key}_original",
                'engineering_method': 'natural',
                'sequence_length': len(sequence),
                'start': -1,  # No engineering region for original sequences
                'end': -1,
                'label': 0  # Natural
            })

    return pd.DataFrame(processed_data)

def create_virus_family_mapping() -> Dict[str, str]:
    """Create mapping of filename patterns to virus families"""
    return {
        'adenovirus': 'adenovirus',
        'chikungunya': 'alphavirus',
        'astrovirus': 'astrovirus',
        'coronavirus': 'coronavirus',
        'ebola': 'filovirus',
        'henipavirus': 'henipavirus',
        'herpesvirus': 'herpesvirus',
        'alphaherpesvirus': 'herpesvirus',
        'influenza': 'influenza',
        'respirovirus': 'paramyxovirus',
        'parainfluenza': 'paramyxovirus',
        'rhinovirus': 'picornavirus',
        'poxvirus': 'poxvirus',
        'rubella': 'rubivirus'
    }

def extract_virus_name_from_filename(filename: str) -> str:
    """Extract virus name from FASTA filename"""
    # Remove common suffixes and convert underscores to spaces
    name = filename.replace('_NC_', ' (NC_').replace('.fasta', '')

    # Add closing parenthesis if NC_ pattern found
    if '(NC_' in name and not name.endswith(')'):
        # Find the NC_ part and add appropriate ending
        parts = name.split('(NC_')
        if len(parts) == 2:
            nc_part = parts[1]
            # Handle different NC_ formats
            if '.' in nc_part:
                name = parts[0] + f"(NC_{nc_part})"
            else:
                name = parts[0] + f"(NC_{nc_part})"

    # Clean up spacing and formatting
    name = re.sub(r'_+', ' ', name)
    name = re.sub(r'\s+', ' ', name)

    return name.strip()

def create_virus_key(virus_name: str) -> str:
    """Create virus key from virus name"""
    # Convert to lowercase, replace spaces and special chars with underscores
    key = virus_name.lower()
    key = re.sub(r'[^\w\s]', '', key)  # Remove punctuation
    key = re.sub(r'\s+', '_', key)     # Replace spaces with underscores
    key = re.sub(r'_+', '_', key)      # Collapse multiple underscores
    return key.strip('_')

def determine_virus_family(filename: str, family_mapping: Dict[str, str]) -> str:
    """Determine virus family from filename using pattern matching"""
    filename_lower = filename.lower()

    for pattern, family in family_mapping.items():
        if pattern in filename_lower:
            return family

    # Default to unknown if no pattern matches
    return 'unknown'

def read_fasta_sequence(fasta_path: Path) -> str:
    """Read sequence from FASTA file"""
    try:
        with open(fasta_path, 'r') as f:
            lines = f.readlines()

        # Skip header line and concatenate sequence lines
        sequence_lines = [line.strip() for line in lines[1:] if line.strip()]
        sequence = ''.join(sequence_lines)

        return sequence.upper()  # Ensure uppercase nucleotides

    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return ""

def create_segments_with_positions(sequences: List[str], n_samples: int, mean_length: int = 187, sd_length: int = 50) -> Tuple[List[str], List[int], List[int]]:
    """
    Generate random segments from viral sequences and return their positions

    Args:
        sequences: List of viral sequences
        n_samples: Number of segments to generate
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length

    Returns:
        Tuple of (sequence_segments, start_positions, end_positions)
    """
    if not sequences:
        raise ValueError("No sequences provided")

    # Pre-allocate arrays for samples
    samples = [None] * n_samples
    start_positions = [None] * n_samples
    end_positions = [None] * n_samples

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
            end = start + length

            samples[batch_start + i] = seq[start:end]
            start_positions[batch_start + i] = start
            end_positions[batch_start + i] = end

    return samples, start_positions, end_positions

def create_segments(sequences: List[str], n_samples: int, mean_length: int = 187, sd_length: int = 50) -> List[str]:
    """
    Generate random segments from viral sequences (backward compatibility)
    """
    segments, _, _ = create_segments_with_positions(sequences, n_samples, mean_length, sd_length)
    return segments

def calculate_overlap(seg_start: int, seg_end: int, eng_start: int, eng_end: int) -> int:
    """
    Calculate overlap between segment and engineering region

    Args:
        seg_start: Segment start position
        seg_end: Segment end position
        eng_start: Engineering region start position
        eng_end: Engineering region end position

    Returns:
        Number of overlapping positions
    """
    if eng_start == -1 or eng_end == -1:  # Original sequence, no engineering region
        return 0

    overlap_start = max(seg_start, eng_start)
    overlap_end = min(seg_end, eng_end)

    return max(0, overlap_end - overlap_start)

def create_virus_dataset(n_samples: int = 1000,
                        engineering_fraction: float = 0.02,
                        mean_length: int = 187,
                        sd_length: int = 50,
                        engineered_data_path: str = None,
                        original_data_path: str = None) -> pd.DataFrame:
    """
    Create dataset with engineered and original virus sequences
    Similar interface to create_engineering_dataset from data_generator.py

    Args:
        n_samples: Total number of samples to generate
        engineering_fraction: Fraction of samples that should be engineered
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length
        engineered_data_path: Path to preprocessed engineered data CSV
        original_data_path: Path to preprocessed original data CSV

    Returns:
        DataFrame with virus dataset
    """
    # Load preprocessed data
    if engineered_data_path and Path(engineered_data_path).exists():
        engineered_df = pd.read_csv(engineered_data_path)
    else:
        raise ValueError("Preprocessed engineered data not found. Run preprocessing first.")

    if original_data_path and Path(original_data_path).exists():
        original_df = pd.read_csv(original_data_path)
    else:
        raise ValueError("Preprocessed original data not found. Run preprocessing first.")

    # Determine sample counts
    n_engineered = int(n_samples * engineering_fraction)
    n_original = n_samples - n_engineered

    print(f"Generating {n_samples} samples ({n_engineered} engineered, {n_original} original)")

    # Sample from engineered data
    if n_engineered > 0:
        # Sample rows from engineered data to get full metadata
        sampled_engineered_rows = engineered_df.sample(n=n_engineered, replace=True, random_state=42)

        # Extract sequences and generate segments with positions
        engineered_sequences = sampled_engineered_rows['sequence'].tolist()
        engineered_segments, seg_starts, seg_ends = create_segments_with_positions(
            engineered_sequences, n_engineered, mean_length, sd_length
        )

        # Calculate overlaps and adjust labels
        labels = []
        for i, (seg_start, seg_end, _, row) in enumerate(zip(
            seg_starts, seg_ends, engineered_segments, sampled_engineered_rows.itertuples()
        )):
            overlap = calculate_overlap(seg_start, seg_end, row.start, row.end)
            # If no specific engineering region is defined (start=-1, end=-1),
            # keep the original label (1) since the whole sequence is engineered
            if row.start == -1 and row.end == -1:
                labels.append(1)  # Whole sequence is engineered
            else:
                # Label as engineered (1) only if overlap is at least 30 positions
                labels.append(1 if overlap >= 30 else 0)

        # Create engineered sample metadata
        engineered_sample_df = pd.DataFrame({
            'sequence_id': [f"ENG_{i:05d}" for i in range(n_engineered)],
            'sequence': engineered_segments,
            'label': labels,
            'engineering_method': sampled_engineered_rows['engineering_method'].values,
            'virus_key': sampled_engineered_rows['virus_key'].values,
            'virus_family': sampled_engineered_rows['virus_family'].values,
            'segment_start': seg_starts,
            'segment_end': seg_ends,
            'eng_start': sampled_engineered_rows['start'].values,
            'eng_end': sampled_engineered_rows['end'].values,
            'overlap': [calculate_overlap(seg_starts[i], seg_ends[i],
                                        sampled_engineered_rows.iloc[i]['start'],
                                        sampled_engineered_rows.iloc[i]['end'])
                       for i in range(n_engineered)]
        })
    else:
        engineered_sample_df = pd.DataFrame()

    # Sample from original data
    if n_original > 0:
        # Sample rows from original data to get full metadata
        sampled_original_rows = original_df.sample(n=n_original, replace=True, random_state=42)

        # Extract sequences and generate segments with positions
        original_sequences = sampled_original_rows['sequence'].tolist()
        original_segments, seg_starts, seg_ends = create_segments_with_positions(
            original_sequences, n_original, mean_length, sd_length
        )

        # Create original sample metadata
        original_sample_df = pd.DataFrame({
            'sequence_id': [f"ORG_{i:05d}" for i in range(n_original)],
            'sequence': original_segments,
            'label': [0] * n_original,  # Always 0 for original sequences
            'engineering_method': ['natural'] * n_original,
            'virus_key': sampled_original_rows['virus_key'].values,
            'virus_family': sampled_original_rows['virus_family'].values,
            'segment_start': seg_starts,
            'segment_end': seg_ends,
            'eng_start': [-1] * n_original,  # No engineering region
            'eng_end': [-1] * n_original,
            'overlap': [0] * n_original  # No overlap possible
        })
    else:
        original_sample_df = pd.DataFrame()

    # Combine datasets
    final_df = pd.concat([engineered_sample_df, original_sample_df], ignore_index=True)

    # Add sequence length and length bins
    final_df['sequence_length'] = final_df['sequence'].str.len()
    final_df['length_bin'] = pd.cut(final_df['sequence_length'],
                                   bins=3,
                                   labels=['short', 'medium', 'long'])

    # Shuffle the dataset
    final_df = final_df.sample(frac=1, random_state=42).reset_index(drop=True)

    return final_df

def preprocess_and_save_data(json_path: str,
                           original_data_dir: str,
                           output_dir: str = "data"):
    """
    Preprocess both engineered and original virus data and save as CSVs

    Args:
        json_path: Path to restriction_virus_engineered.json
        original_data_dir: Path to directory with original FASTA files
        output_dir: Directory to save processed CSV files
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Process engineered data
    print("Processing engineered virus data...")
    engineered_df = load_engineered_data(json_path)
    engineered_output = output_path / "processed_engineered_virus.csv"
    engineered_df.to_csv(engineered_output, index=False)
    print(f"Saved {len(engineered_df)} engineered entries to {engineered_output}")

    # Process original data
    print("\nProcessing original virus data...")
    original_df = load_original_virus_data(original_data_dir)
    original_output = output_path / "processed_original_virus.csv"
    original_df.to_csv(original_output, index=False)
    print(f"Saved {len(original_df)} original entries to {original_output}")

    print(f"\nPreprocessing complete! CSV files saved in {output_dir}/")
    return str(engineered_output), str(original_output)

if __name__ == "__main__":
    # Example usage
    json_path = "data/restriction_virus_engineered.json"
    original_data_dir = "data/original_virus"

    # Preprocess data
    eng_path, orig_path = preprocess_and_save_data(json_path, original_data_dir)

    # Create sample dataset
    sample_df = create_virus_dataset(
        n_samples=1000,
        engineering_fraction=0.02,
        engineered_data_path=eng_path,
        original_data_path=orig_path
    )

    print(f"\nGenerated dataset with {len(sample_df)} samples")
    print(f"Engineering fraction: {sample_df['label'].mean():.3f}")
    print("\nDataset columns:", list(sample_df.columns))
    print("\nSample entries:")
    print(sample_df.head())