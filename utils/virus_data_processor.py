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

def load_split_mapping(split_dataset_path: str) -> Dict[str, str]:
    """
    Load split dataset and create accession_id -> split_id mapping

    Args:
        split_dataset_path: Path to family_split_dataset.csv or random_split_dataset.csv

    Returns:
        Dictionary mapping accession_id to split_id (train/test/val)
    """
    print(f"Loading split mapping from {split_dataset_path}...")

    split_df = pd.read_csv(split_dataset_path)
    split_mapping = dict(zip(split_df['accession_id'], split_df['split_id']))

    print(f"Loaded split mapping for {len(split_mapping)} accession IDs")
    print(f"Split distribution: {split_df['split_id'].value_counts().to_dict()}")

    return split_mapping

def extract_accession_from_virus_name(virus_name: str) -> Optional[str]:
    """
    Extract accession ID from virus name like 'Human rhinovirus A (NC_001490.1)'

    Args:
        virus_name: Virus name containing accession ID in parentheses

    Returns:
        Accession ID or None if not found
    """
    import re
    match = re.search(r'\((NC_[^)]+)\)', virus_name)
    return match.group(1) if match else None

def extract_accession_from_filename(filename: str) -> Optional[str]:
    """
    Extract accession ID from FASTA filename like 'Human_rhinovirus_A_NC_001490.1.fasta'

    Args:
        filename: FASTA filename

    Returns:
        Accession ID or None if not found
    """
    import re
    match = re.search(r'(NC_[^.]+\.[^.]+)', filename)
    return match.group(1) if match else None

def load_engineered_data(json_path: str, split_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Load and process engineered virus data from JSON file

    Args:
        json_path: Path to the restriction_virus_engineered.json file
        split_mapping: Optional dict mapping accession_id to split_id

    Returns:
        DataFrame with processed engineered virus data (includes split_id if mapping provided)
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

        # Extract accession ID and get split_id if mapping provided
        accession_id = extract_accession_from_virus_name(entry['virus'])
        split_id = split_mapping.get(accession_id, 'unknown') if split_mapping and accession_id else 'unknown'

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
            'accession_id': accession_id,
            'split_id': split_id,
            'label': 1  # Engineered
        })

    return pd.DataFrame(processed_data)

def load_original_virus_data(data_dir: str, split_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Load and process original virus data from FASTA files

    Args:
        data_dir: Path to directory containing FASTA files
        split_mapping: Optional dict mapping accession_id to split_id

    Returns:
        DataFrame with processed original virus data (includes split_id if mapping provided)
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

            # Extract accession ID and get split_id if mapping provided
            accession_id = extract_accession_from_filename(filename)
            split_id = split_mapping.get(accession_id, 'unknown') if split_mapping and accession_id else 'unknown'

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
                'accession_id': accession_id,
                'split_id': split_id,
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

def create_virus_dataset_with_splits(n_samples_per_split: Dict[str, int],
                                   engineering_fraction: float = 0.02,
                                   mean_length: int = 187,
                                   sd_length: int = 50,
                                   engineered_data_path: str = None,
                                   original_data_path: str = None) -> Dict[str, pd.DataFrame]:
    """
    Create datasets split by train/test/val with engineered and original virus sequences

    Args:
        n_samples_per_split: Dict with keys train/test/val and values as sample counts
        engineering_fraction: Fraction of samples that should be engineered
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length
        engineered_data_path: Path to preprocessed engineered data CSV
        original_data_path: Path to preprocessed original data CSV

    Returns:
        Dict with keys train/test/val and DataFrames as values
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

    split_datasets = {}

    for split_name, n_samples in n_samples_per_split.items():
        print(f"\nGenerating {n_samples} samples for {split_name} split...")

        # Filter data for this split
        split_engineered = engineered_df[engineered_df['split_id'] == split_name]
        split_original = original_df[original_df['split_id'] == split_name]

        print(f"  Available: {len(split_engineered)} engineered, {len(split_original)} original")

        # Determine sample counts
        n_engineered = int(n_samples * engineering_fraction)
        n_original = n_samples - n_engineered

        print(f"  Targeting: {n_engineered} engineered, {n_original} original")

        # Sample from engineered data for this split
        if n_engineered > 0 and len(split_engineered) > 0:
            sampled_engineered_rows = split_engineered.sample(n=min(n_engineered, len(split_engineered)),
                                                             replace=True, random_state=42)

            engineered_sequences = sampled_engineered_rows['sequence'].tolist()
            engineered_segments, seg_starts, seg_ends = create_segments_with_positions(
                engineered_sequences, len(sampled_engineered_rows), mean_length, sd_length
            )

            # Calculate overlaps and adjust labels
            labels = []
            for i, (seg_start, seg_end, _, row) in enumerate(zip(
                seg_starts, seg_ends, engineered_segments, sampled_engineered_rows.itertuples()
            )):
                overlap = calculate_overlap(seg_start, seg_end, row.start, row.end)
                if row.start == -1 and row.end == -1:
                    labels.append(1)  # Whole sequence is engineered
                else:
                    labels.append(1 if overlap >= 30 else 0)

            engineered_sample_df = pd.DataFrame({
                'sequence_id': [f"{split_name.upper()}_ENG_{i:05d}" for i in range(len(engineered_segments))],
                'sequence': engineered_segments,
                'label': labels,
                'split_id': [split_name] * len(engineered_segments),
                'engineering_method': sampled_engineered_rows['engineering_method'].values,
                'virus_key': sampled_engineered_rows['virus_key'].values,
                'virus_family': sampled_engineered_rows['virus_family'].values,
                'accession_id': sampled_engineered_rows['accession_id'].values,
                'segment_start': seg_starts,
                'segment_end': seg_ends,
                'eng_start': sampled_engineered_rows['start'].values,
                'eng_end': sampled_engineered_rows['end'].values,
                'overlap': [calculate_overlap(seg_starts[i], seg_ends[i],
                                            sampled_engineered_rows.iloc[i]['start'],
                                            sampled_engineered_rows.iloc[i]['end'])
                           for i in range(len(engineered_segments))]
            })
        else:
            engineered_sample_df = pd.DataFrame()

        # Sample from original data for this split
        if n_original > 0 and len(split_original) > 0:
            sampled_original_rows = split_original.sample(n=min(n_original, len(split_original)),
                                                         replace=True, random_state=42)

            original_sequences = sampled_original_rows['sequence'].tolist()
            original_segments, seg_starts, seg_ends = create_segments_with_positions(
                original_sequences, len(sampled_original_rows), mean_length, sd_length
            )

            original_sample_df = pd.DataFrame({
                'sequence_id': [f"{split_name.upper()}_ORG_{i:05d}" for i in range(len(original_segments))],
                'sequence': original_segments,
                'label': [0] * len(original_segments),
                'split_id': [split_name] * len(original_segments),
                'engineering_method': ['natural'] * len(original_segments),
                'virus_key': sampled_original_rows['virus_key'].values,
                'virus_family': sampled_original_rows['virus_family'].values,
                'accession_id': sampled_original_rows['accession_id'].values,
                'segment_start': seg_starts,
                'segment_end': seg_ends,
                'eng_start': [-1] * len(original_segments),
                'eng_end': [-1] * len(original_segments),
                'overlap': [0] * len(original_segments)
            })
        else:
            original_sample_df = pd.DataFrame()

        # Combine datasets for this split
        split_df = pd.concat([engineered_sample_df, original_sample_df], ignore_index=True)

        if len(split_df) > 0:
            # Add sequence length and length bins
            split_df['sequence_length'] = split_df['sequence'].str.len()
            split_df['length_bin'] = pd.cut(split_df['sequence_length'],
                                           bins=3,
                                           labels=['short', 'medium', 'long'])

            # Shuffle the dataset
            split_df = split_df.sample(frac=1, random_state=42).reset_index(drop=True)

            print(f"  Generated {len(split_df)} samples for {split_name}")
            print(f"  Engineering fraction: {split_df['label'].mean():.3f}")

        split_datasets[split_name] = split_df

    return split_datasets

def preprocess_and_save_data(json_path: str,
                           original_data_dir: str,
                           output_dir: str = "data",
                           split_dataset_path: str = None,
                           split_type: str = None):
    """
    Preprocess both engineered and original virus data and save as CSVs

    Args:
        json_path: Path to restriction_virus_engineered.json
        original_data_dir: Path to directory with original FASTA files
        output_dir: Directory to save processed CSV files
        split_dataset_path: Optional path to split dataset (family_split_dataset.csv or random_split_dataset.csv)
        split_type: Optional split type identifier ('family' or 'random') for output naming

    Returns:
        Tuple of (engineered_output_path, original_output_path)
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Load split mapping if provided
    split_mapping = None
    if split_dataset_path:
        split_mapping = load_split_mapping(split_dataset_path)

    # Process engineered data
    print("Processing engineered virus data...")
    engineered_df = load_engineered_data(json_path, split_mapping)

    # Choose output filename based on split type
    if split_type:
        engineered_output = output_path / f"processed_engineered_virus_{split_type}.csv"
    else:
        engineered_output = output_path / "processed_engineered_virus.csv"

    engineered_df.to_csv(engineered_output, index=False)
    print(f"Saved {len(engineered_df)} engineered entries to {engineered_output}")

    # Process original data
    print("\nProcessing original virus data...")
    original_df = load_original_virus_data(original_data_dir, split_mapping)

    # Choose output filename based on split type
    if split_type:
        original_output = output_path / f"processed_original_virus_{split_type}.csv"
    else:
        original_output = output_path / "processed_original_virus.csv"

    original_df.to_csv(original_output, index=False)
    print(f"Saved {len(original_df)} original entries to {original_output}")

    # Print split distribution if available
    if split_mapping:
        print(f"\nSplit distribution in engineered data:")
        eng_split_counts = engineered_df['split_id'].value_counts()
        print(eng_split_counts)

        print(f"\nSplit distribution in original data:")
        orig_split_counts = original_df['split_id'].value_counts()
        print(orig_split_counts)

    print(f"\nPreprocessing complete! CSV files saved in {output_dir}/")
    return str(engineered_output), str(original_output)

def preprocess_and_create_split_datasets(json_path: str,
                                        original_data_dir: str,
                                        split_dataset_path: str,
                                        split_type: str,
                                        output_dir: str = "data",
                                        n_samples_per_split: Dict[str, int] = None,
                                        engineering_fraction: float = 0.02):
    """
    Complete pipeline to preprocess data with split information and create train/test/val datasets

    Args:
        json_path: Path to restriction_virus_engineered.json
        original_data_dir: Path to directory with original FASTA files
        split_dataset_path: Path to split dataset (family_split_dataset.csv or random_split_dataset.csv)
        split_type: Split type identifier ('family' or 'random')
        output_dir: Directory to save processed CSV files
        n_samples_per_split: Dict with keys train/test/val and values as sample counts
        engineering_fraction: Fraction of samples that should be engineered

    Returns:
        Dict with keys train/test/val and DataFrames as values
    """
    if n_samples_per_split is None:
        n_samples_per_split = {'train': 8000, 'test': 1000, 'val': 1000}

    print(f"=== Processing virus data with {split_type} splits ===")

    # Preprocess data with split information
    eng_path, orig_path = preprocess_and_save_data(
        json_path=json_path,
        original_data_dir=original_data_dir,
        output_dir=output_dir,
        split_dataset_path=split_dataset_path,
        split_type=split_type
    )

    # Create split datasets
    print(f"\n=== Creating split datasets ===")
    split_datasets = create_virus_dataset_with_splits(
        n_samples_per_split=n_samples_per_split,
        engineering_fraction=engineering_fraction,
        engineered_data_path=eng_path,
        original_data_path=orig_path
    )

    # Save split datasets
    output_path = Path(output_dir)
    saved_paths = {}

    for split_name, df in split_datasets.items():
        if len(df) > 0:
            output_file = output_path / f"virus_dataset_{split_type}_{split_name}.csv"
            df.to_csv(output_file, index=False)
            saved_paths[split_name] = str(output_file)
            print(f"\nSaved {split_name} dataset: {len(df)} samples to {output_file}")
            print(f"  Engineering fraction: {df['label'].mean():.3f}")
        else:
            print(f"\nWarning: No samples generated for {split_name} split")

    return split_datasets, saved_paths

if __name__ == "__main__":
    # File paths
    json_path = "data/restriction_virus_engineered.json"
    original_data_dir = "data/original_virus"

    # Example 1: Original functionality (backward compatibility)
    print("=== Example 1: Original functionality ===")
    eng_path, orig_path = preprocess_and_save_data(json_path, original_data_dir)

    sample_df = create_virus_dataset(
        n_samples=1000,
        engineering_fraction=0.02,
        engineered_data_path=eng_path,
        original_data_path=orig_path
    )

    print(f"Generated dataset with {len(sample_df)} samples")
    print(f"Engineering fraction: {sample_df['label'].mean():.3f}")

    # Example 2: Family split functionality
    print(f"\n{'='*60}")
    print("=== Example 2: Family split functionality ===")
    family_split_datasets, family_paths = preprocess_and_create_split_datasets(
        json_path=json_path,
        original_data_dir=original_data_dir,
        split_dataset_path="data/family_split_dataset.csv",
        split_type="family",
        n_samples_per_split={'train': 800, 'test': 100, 'val': 100},
        engineering_fraction=0.02
    )

    # Example 3: Random split functionality
    print(f"\n{'='*60}")
    print("=== Example 3: Random split functionality ===")
    random_split_datasets, random_paths = preprocess_and_create_split_datasets(
        json_path=json_path,
        original_data_dir=original_data_dir,
        split_dataset_path="data/random_split_dataset.csv",
        split_type="random",
        n_samples_per_split={'train': 800, 'test': 100, 'val': 100},
        engineering_fraction=0.02
    )

    print(f"\n{'='*60}")
    print("=== Summary ===")
    print("Original functionality still works - no breaking changes")
    print(f"Family split datasets saved: {list(family_paths.keys())}")
    print(f"Random split datasets saved: {list(random_paths.keys())}")
    print("\nSplit-aware datasets include accession_id and split_id columns for reproducible train/test splits")