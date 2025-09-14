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

def extract_start_end_positions(entry: Dict) -> Tuple[int, int]:
    """
    Extract start and end positions where engineered sequence differs from original

    For each engineering method, returns the span in the final sequence that differs
    from what would be expected in the original sequence:

    1. region_inversion: The inverted region (start, end)
    2. region_duplication: The inserted duplicated sequence (insert_position, insert_position + duplicated_length)
    3. region_deletion: The deletion point in final sequence (deleted_start, deleted_start - effectively point)
    4. insertions: The inserted sequence span (insert_position, insert_position + inserted_length)
    5. point mutations: The mutated position (position, position)
    6. random_substitution: The substituted region (substitution_start, substitution_end)

    Returns:
        Tuple of (start, end) positions where sequence differs from original
    """
    method = entry.get('method', '')

    # region_inversion: inverted region differs from original
    if method == 'region_inversion':
        return (entry['start'], entry['end'])

    # region_duplication: inserted duplicate creates new sequence span
    elif method == 'region_duplication':
        insert_pos = entry['insert_position']
        duplicated_length = entry['source_end'] - entry['source_start']
        return (insert_pos, insert_pos + duplicated_length)

    # region_deletion: in final sequence, deletion point has no span (collapsed)
    elif method == 'region_deletion':
        # In the final sequence, the deletion creates a boundary point
        return (entry['deleted_start'], entry['deleted_start'] + 100)

    # Various insertions: inserted sequence creates new span
    elif method in ['gfp_insertion', 'tn5_insertion', 'igem_gfp_insertion', 'crispr_insertion', 'fallback_insertion']:
        insert_pos = entry['insert_position']
        inserted_length = entry.get('inserted_length', 0)
        return (insert_pos, insert_pos + inserted_length)

    # Point mutations: single position differs
    elif method in ['frameshift_mutation', 'single_nucleotide_change']:
        pos = entry['position']
        return (pos, pos)

    # random_substitution: substituted region differs
    elif method == 'random_substitution':
        return (entry['substitution_start'], entry['substitution_end'])

    else:
        # Fallback: try to detect pattern automatically
        if 'start' in entry and 'end' in entry:
            return (entry['start'], entry['end'])
        elif 'insert_position' in entry:
            insert_pos = entry['insert_position']
            inserted_length = entry.get('inserted_length', 0)
            return (insert_pos, insert_pos + inserted_length)
        elif 'position' in entry:
            pos = entry['position']
            return (pos, pos)
        elif 'deleted_start' in entry and 'deleted_end' in entry:
            return (entry['deleted_start'], entry['deleted_start'])
        elif 'substitution_start' in entry and 'substitution_end' in entry:
            return (entry['substitution_start'], entry['substitution_end'])
        else:
            print(f"Warning: Unknown coordinate pattern for method '{method}': {sorted(entry.keys())}")
            return (-1, -1)

def load_engineered_data(json_path: str, family_split_path: str = None, random_split_path: str = None) -> pd.DataFrame:
    """
    Load and process engineered virus data from JSON file

    Args:
        json_path: Path to the restriction_virus_engineered.json file
        family_split_path: Path to family_split_dataset.csv for split information
        random_split_path: Path to random_split_dataset.csv for split information

    Returns:
        DataFrame with processed engineered virus data including split columns
    """
    print(f"Loading engineered data from {json_path}...")

    with open(json_path, 'r') as f:
        data = json.load(f)

    print(f"Loaded {len(data)} engineered entries")

    # Load split information if paths provided
    family_splits = {}
    random_splits = {}
    
    if family_split_path and Path(family_split_path).exists():
        print("Loading family split information...")
        family_df = pd.read_csv(family_split_path)
        family_splits = dict(zip(family_df['accession_id'], family_df['split_id']))
        print(f"Loaded {len(family_splits)} family split mappings")
    
    if random_split_path and Path(random_split_path).exists():
        print("Loading random split information...")
        random_df = pd.read_csv(random_split_path)
        random_splits = dict(zip(random_df['accession_id'], random_df['split_id']))
        print(f"Loaded {len(random_splits)} random split mappings")

    # Extract relevant fields
    processed_data = []
    
    for entry in data:
        # Extract start/end positions where sequence differs from original
        start_pos, end_pos = extract_start_end_positions(entry)

        # Extract accession ID from virus name
        virus_name = entry['virus']
        accession_id = extract_accession_id_from_virus_name(virus_name)
        
        # Get split information
        family_split = family_splits.get(accession_id, None)
        random_split = random_splits.get(accession_id, None)

        processed_data.append({
            'sequence': entry['full_genome'],
            'virus': virus_name,
            'virus_key': entry['virus_key'],
            'virus_family': entry['virus_family'],
            'example_id': entry['example_id'],
            'engineering_method': entry['method'],
            'sequence_length': len(entry['full_genome']),
            'start': start_pos,
            'end': end_pos,
            'label': 1,  # Engineered
            'accession_id': accession_id,
            'family_split': family_split,
            'random_split': random_split
        })

    return pd.DataFrame(processed_data)

def extract_accession_id_from_virus_name(virus_name: str) -> str:
    """Extract accession ID from virus name string"""
    # Look for accession ID pattern like (NC_001490.1) or (MT160087.1)
    match = re.search(r'\(([A-Z]{2}_[0-9]+\.[0-9]+)\)', virus_name)
    if match:
        return match.group(1)
    return None

def extract_accession_id_from_filename(filename: str) -> str:
    """Extract accession ID from FASTA filename"""
    # Look for accession ID pattern like NC_001498.1 or MT160087.1
    match = re.search(r'([A-Z]{2}_[0-9]+\.[0-9]+)', filename)
    if match:
        return match.group(1)
    return None

def load_original_virus_data(data_dir: str, family_split_path: str = None, random_split_path: str = None) -> pd.DataFrame:
    """
    Load and process original virus data from FASTA files

    Args:
        data_dir: Path to directory containing FASTA files
        family_split_path: Path to family_split_dataset.csv for split information
        random_split_path: Path to random_split_dataset.csv for split information

    Returns:
        DataFrame with processed original virus data including split columns
    """
    print(f"Loading original virus data from {data_dir}...")

    data_path = Path(data_dir)
    fasta_files = list(data_path.glob("*.fasta"))

    print(f"Found {len(fasta_files)} FASTA files")

    # Load split information if paths provided
    family_splits = {}
    random_splits = {}
    
    if family_split_path and Path(family_split_path).exists():
        print("Loading family split information...")
        family_df = pd.read_csv(family_split_path)
        family_splits = dict(zip(family_df['accession_id'], family_df['split_id']))
        print(f"Loaded {len(family_splits)} family split mappings")
    
    if random_split_path and Path(random_split_path).exists():
        print("Loading random split information...")
        random_df = pd.read_csv(random_split_path)
        random_splits = dict(zip(random_df['accession_id'], random_df['split_id']))
        print(f"Loaded {len(random_splits)} random split mappings")

    processed_data = []

    # Create virus family mapping based on common patterns
    family_mapping = create_virus_family_mapping()

    for fasta_file in fasta_files:
        # Extract virus name from filename
        filename = fasta_file.stem
        virus_name = extract_virus_name_from_filename(filename)
        virus_key = create_virus_key(virus_name)

        # Extract accession ID from filename
        accession_id = extract_accession_id_from_filename(filename)
        
        # Get split information
        family_split = family_splits.get(accession_id, None)
        random_split = random_splits.get(accession_id, None)

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
                'label': 0,  # Natural
                'accession_id': accession_id,
                'family_split': family_split,
                'random_split': random_split
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

def load_ged_test_dataset(csv_path: str, n_samples: int = None, mean_length: int = 187, sd_length: int = 50) -> pd.DataFrame:
    """
    Load GED_test.csv and generate segments compatible with train_and_test.py format

    Args:
        csv_path: Path to GED_test.csv
        n_samples: Number of samples to generate (None for all)
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length

    Returns:
        DataFrame compatible with train_and_test.py test dataset format
    """
    print(f"Loading GED test data from {csv_path}...")

    # Load the CSV
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} sequences from GED test dataset")

    sequences = df['sequence'].tolist()
    labels = df['engineered'].tolist()

    # If n_samples not specified, use all data
    if n_samples is None:
        n_samples = len(sequences)
    elif n_samples > len(sequences):
        # If more samples requested than available, sample with replacement
        print(f"Requested {n_samples} samples but only {len(sequences)} available, sampling with replacement")
        sample_indices = np.random.choice(len(sequences), n_samples, replace=True)
        sequences = [sequences[i] for i in sample_indices]
        labels = [labels[i] for i in sample_indices]
    else:
        # Sample without replacement
        sample_indices = np.random.choice(len(sequences), n_samples, replace=False)
        sequences = [sequences[i] for i in sample_indices]
        labels = [labels[i] for i in sample_indices]

    # Generate segments with positions
    segments, seg_starts, seg_ends = create_segments_with_positions(
        sequences, len(sequences), mean_length, sd_length
    )

    # Create DataFrame in expected format
    result_df = pd.DataFrame({
        'sequence_id': [f"GED_{i:05d}" for i in range(len(segments))],
        'sequence': segments,
        'label': labels,
        'engineering_method': ['unknown' if label == 0 else 'engineered' for label in labels],
        'virus_key': ['ged_test'] * len(segments),
        'virus_family': ['unknown'] * len(segments),
        'segment_start': seg_starts,
        'segment_end': seg_ends,
        'eng_start': [-1] * len(segments),  # Unknown engineering positions
        'eng_end': [-1] * len(segments),
        'overlap': [0] * len(segments),
        'accession_id': ['unknown'] * len(segments),
        'family_split': [None] * len(segments),
        'random_split': [None] * len(segments),
        'sequence_length': [len(seg) for seg in segments],
        'length_bin': pd.cut([len(seg) for seg in segments], bins=3, labels=['short', 'medium', 'long'])
    })

    print(f"Generated {len(result_df)} segments with engineering fraction: {result_df['label'].mean():.3f}")
    return result_df

def load_viral_vectors_dataset(csv_path: str, n_samples: int = None, mean_length: int = 187, sd_length: int = 50) -> pd.DataFrame:
    """
    Load viral_vectors.csv and generate segments compatible with train_and_test.py format
    All sequences get label=1 (engineered)

    Args:
        csv_path: Path to viral_vectors.csv
        n_samples: Number of samples to generate (None for all)
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length

    Returns:
        DataFrame compatible with train_and_test.py test dataset format
    """
    print(f"Loading viral vectors data from {csv_path}...")

    # Load the CSV
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} sequences from viral vectors dataset")

    sequences = df['sequence'].tolist()
    names = df['name'].tolist() if 'name' in df.columns else ['unknown'] * len(sequences)

    # If n_samples not specified, use all data
    if n_samples is None:
        n_samples = len(sequences)
    elif n_samples > len(sequences):
        # If more samples requested than available, sample with replacement
        print(f"Requested {n_samples} samples but only {len(sequences)} available, sampling with replacement")
        sample_indices = np.random.choice(len(sequences), n_samples, replace=True)
        sequences = [sequences[i] for i in sample_indices]
        names = [names[i] for i in sample_indices]
    else:
        # Sample without replacement
        sample_indices = np.random.choice(len(sequences), n_samples, replace=False)
        sequences = [sequences[i] for i in sample_indices]
        names = [names[i] for i in sample_indices]

    # Generate segments with positions
    segments, seg_starts, seg_ends = create_segments_with_positions(
        sequences, len(sequences), mean_length, sd_length
    )

    # Create DataFrame in expected format - all labels are 1 (engineered)
    result_df = pd.DataFrame({
        'sequence_id': [f"VV_{i:05d}" for i in range(len(segments))],
        'sequence': segments,
        'label': [1] * len(segments),  # All viral vectors are labeled as engineered
        'engineering_method': ['viral_vector'] * len(segments),
        'virus_key': ['viral_vector'] * len(segments),
        'virus_family': ['viral_vector'] * len(segments),
        'segment_start': seg_starts,
        'segment_end': seg_ends,
        'eng_start': [-1] * len(segments),  # Unknown engineering positions
        'eng_end': [-1] * len(segments),
        'overlap': [0] * len(segments),
        'accession_id': ['unknown'] * len(segments),
        'family_split': [None] * len(segments),
        'random_split': [None] * len(segments),
        'sequence_length': [len(seg) for seg in segments],
        'length_bin': pd.cut([len(seg) for seg in segments], bins=3, labels=['short', 'medium', 'long'])
    })

    print(f"Generated {len(result_df)} segments, all labeled as engineered (label=1)")
    return result_df

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
                        original_data_path: str = None,
                        split_filter: str = None,
                        split_type: str = 'family',
                        engineering_method_filter: str = None) -> pd.DataFrame:
    """
    Create dataset with engineered and original virus sequences
    Similar interface to create_engineering_dataset from data_generator.py

    Args:
        n_samples: Total number of samples to generate
        engineering_fraction: Fraction of samples that should have label=1 in final output
        mean_length: Mean segment length
        sd_length: Standard deviation of segment length
        engineered_data_path: Path to preprocessed engineered data CSV
        original_data_path: Path to preprocessed original data CSV
        split_filter: Filter by split ('train' or 'test'), None for no filtering
        split_type: Type of split to filter by ('family' or 'random')
        engineering_method_filter: If provided, only use engineered sequences with this engineering method

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

    # Apply split filtering if specified
    if split_filter is not None:
        split_column = f"{split_type}_split"
        if split_column not in engineered_df.columns or split_column not in original_df.columns:
            raise ValueError(f"Split column '{split_column}' not found in data. Run preprocessing with split information first.")

        print(f"Filtering data by {split_type} split: {split_filter}")
        engineered_df = engineered_df[engineered_df[split_column] == split_filter].copy()
        original_df = original_df[original_df[split_column] == split_filter].copy()

        print(f"After filtering: {len(engineered_df)} engineered, {len(original_df)} original samples")

    # Filter for specific engineering method if provided
    if engineering_method_filter:
        print(f"Filtering engineered data for {engineering_method_filter} method only...")
        engineered_df = engineered_df[engineered_df['engineering_method'] == engineering_method_filter].copy()
        print(f"After engineering method filter: {len(engineered_df)} engineered samples")

    target_positive_samples = int(n_samples * engineering_fraction)
    n_original = n_samples - target_positive_samples

    print(f"Generating {n_samples} samples (target: {target_positive_samples} positive, {n_original} original)")

    # Oversample engineered sequences until we get target_positive_samples with label=1
    positive_samples = []
    batch_size = max(100, target_positive_samples * 100)  # Start with 100x oversampling

    while len(positive_samples) < target_positive_samples:
        # Sample batch from engineered data
        batch_rows = engineered_df.sample(n=batch_size, replace=True)
        batch_sequences = batch_rows['sequence'].tolist()
        batch_segments, seg_starts, seg_ends = create_segments_with_positions(
            batch_sequences, batch_size, mean_length, sd_length
        )

        # Keep only segments that result in label=1
        for seg_start, seg_end, segment, row in zip(
            seg_starts, seg_ends, batch_segments, batch_rows.itertuples()
        ):
            if len(positive_samples) >= target_positive_samples:
                break

            overlap = calculate_overlap(seg_start, seg_end, row.start, row.end)
            label = 1 if (row.start == -1 and row.end == -1) or overlap >= 30 else 0

            if label == 1:
                positive_samples.append({
                    'sequence_id': f"ENG_{len(positive_samples):05d}",
                    'sequence': segment,
                    'label': 1,
                    'engineering_method': row.engineering_method,
                    'virus_key': row.virus_key,
                    'virus_family': row.virus_family,
                    'segment_start': seg_start,
                    'segment_end': seg_end,
                    'eng_start': row.start,
                    'eng_end': row.end,
                    'overlap': overlap,
                    'accession_id': row.accession_id,
                    'family_split': row.family_split,
                    'random_split': row.random_split
                })

    engineered_sample_df = pd.DataFrame(positive_samples)

    # Sample from original data
    if n_original > 0:
        sampled_original_rows = original_df.sample(n=n_original, replace=True, random_state=42)
        original_sequences = sampled_original_rows['sequence'].tolist()
        original_segments, seg_starts, seg_ends = create_segments_with_positions(
            original_sequences, n_original, mean_length, sd_length
        )

        original_sample_df = pd.DataFrame({
            'sequence_id': [f"ORG_{i:05d}" for i in range(n_original)],
            'sequence': original_segments,
            'label': [0] * n_original,
            'engineering_method': ['natural'] * n_original,
            'virus_key': sampled_original_rows['virus_key'].values,
            'virus_family': sampled_original_rows['virus_family'].values,
            'segment_start': seg_starts,
            'segment_end': seg_ends,
            'eng_start': [-1] * n_original,
            'eng_end': [-1] * n_original,
            'overlap': [0] * n_original,
            'accession_id': sampled_original_rows['accession_id'].values,
            'family_split': sampled_original_rows['family_split'].values,
            'random_split': sampled_original_rows['random_split'].values
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
                           output_dir: str = "data",
                           family_split_path: str = None,
                           random_split_path: str = None):
    """
    Preprocess both engineered and original virus data and save as CSVs

    Args:
        json_path: Path to restriction_virus_engineered.json
        original_data_dir: Path to directory with original FASTA files
        output_dir: Directory to save processed CSV files
        family_split_path: Path to family_split_dataset.csv for split information
        random_split_path: Path to random_split_dataset.csv for split information
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Process engineered data
    print("Processing engineered virus data...")
    engineered_df = load_engineered_data(json_path, family_split_path, random_split_path)
    engineered_output = output_path / "processed_engineered_virus.csv"
    engineered_df.to_csv(engineered_output, index=False)
    print(f"Saved {len(engineered_df)} engineered entries to {engineered_output}")

    # Process original data
    print("\nProcessing original virus data...")
    original_df = load_original_virus_data(original_data_dir, family_split_path, random_split_path)
    original_output = output_path / "processed_original_virus.csv"
    original_df.to_csv(original_output, index=False)
    print(f"Saved {len(original_df)} original entries to {original_output}")

    print(f"\nPreprocessing complete! CSV files saved in {output_dir}/")
    return str(engineered_output), str(original_output)

if __name__ == "__main__":
    # Example usage
    json_path = "data/restriction_virus_engineered.json"
    original_data_dir = "data/original_virus"
    family_split_path = "data/family_split_dataset.csv"
    random_split_path = "data/random_split_dataset.csv"

    # Preprocess data
    eng_path, orig_path = preprocess_and_save_data(
        json_path, 
        original_data_dir,
        family_split_path=family_split_path,
        random_split_path=random_split_path
    )

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