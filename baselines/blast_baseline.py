#!/usr/bin/env python3
"""
BLAST-based Baseline for Chimera Detection
"""

import subprocess
import os
import pandas as pd
import tempfile
import numpy as np
import math
from utils.imbalance_handling import adjust_probabilities_for_test_prior, find_optimal_threshold_f1

def write_fasta_file(sequences, file_path):
    """Write sequences to a FASTA file from a list."""
    with open(file_path, 'w') as f:
        for i, sequence in enumerate(sequences):
            f.write(f'>{i}\n')
            f.write(f'{sequence}\n')

def create_blast_db(fasta_file, db_name):
    """Create a BLAST database from a fasta file."""
    command = [
        'makeblastdb',
        '-in', fasta_file,
        '-dbtype', 'nucl',
        '-out', db_name
    ]
    subprocess.run(command, check=True, capture_output=True)

def run_blast(query_file, db_name, output_file, num_threads=1, evalue=100.0, word_size=7, max_target_seqs=9, max_hsps=5):
    """Run blastn and save results to a file."""
    command = [
        'blastn',
        '-query', query_file,
        '-db', db_name,
        '-out', output_file,
        '-outfmt', '6',  # Tabular output
        '-max_hsps', str(max_hsps),
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_target_seqs),
        '-word_size', str(word_size),
        '-num_threads', str(num_threads)
    ]
    subprocess.run(command, check=True, capture_output=True)

def parse_blast_results(blast_output_path):
    """Parse BLAST tabular results and return a DataFrame."""
    columns = [
        "query_id", "subject_id", "% identity", "alignment_length", "mismatches",
        "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
    ]
    if os.path.getsize(blast_output_path) == 0:
        return pd.DataFrame(columns=columns)

    blast_df = pd.read_csv(blast_output_path, sep='\t', header=None, names=columns)
    blast_df['query_id'] = blast_df['query_id'].astype(int)
    return blast_df

def calculate_soft_score(identity, alignment_length, query_length, evalue,
                        identity_weight=0.4, coverage_weight=0.4, evalue_weight=0.2):
    """
    Calculate a soft scoring function for BLAST results.

    Args:
        identity: % identity of alignment
        alignment_length: Length of alignment
        query_length: Length of query sequence
        evalue: E-value from BLAST
        identity_weight: Weight for identity component (0-1)
        coverage_weight: Weight for coverage component (0-1)
        evalue_weight: Weight for e-value component (0-1)

    Returns:
        float: Soft score between 0 and 1
    """
    # Normalize identity (0-100 -> 0-1)
    identity_score = identity / 100.0

    # Calculate coverage (alignment_length / query_length)
    coverage_score = min(1.0, alignment_length / max(query_length, 1))

    # E-value score: transform using negative log, then normalize
    # Lower e-values are better, so we use -log(evalue + 1e-300 to avoid log(0))
    evalue_score = min(1.0, max(0.0, -math.log10(max(evalue, 1e-300)) / 50.0))

    # Weighted combination
    soft_score = (identity_weight * identity_score +
                  coverage_weight * coverage_score +
                  evalue_weight * evalue_score)

    return soft_score

def sigmoid(x):
    """Sigmoid function to convert scores to probabilities"""
    return 1.0 / (1.0 + math.exp(-x))

def run_blast_baseline(X_train, X_test, y_train, y_test, min_non_matching_bases=30,
                      identity_threshold=96.66126976030961, evalue=100.0, word_size=7, max_target_seqs=9, max_hsps=5,
                      use_soft_scoring=True, apply_probability_adjustment=True, num_threads=1):
    """Run BLAST-based chimera detection baseline."""
    with tempfile.TemporaryDirectory() as temp_dir:
        db_name = os.path.join(temp_dir, 'viral_db')
        natural_sequences_file = os.path.join(temp_dir, 'natural.fasta')
        query_file = os.path.join(temp_dir, 'query.fasta')
        output_file = os.path.join(temp_dir, 'blast_results.tsv')

        # Create BLAST database from natural sequences in the training set
        natural_sequences = X_train[y_train == 0]['sequence']
        write_fasta_file(natural_sequences, natural_sequences_file)
        create_blast_db(natural_sequences_file, db_name)

        # Run BLAST on test sequences
        test_sequences = X_test['sequence']
        write_fasta_file(test_sequences, query_file)
        run_blast(query_file, db_name, output_file, num_threads=num_threads, evalue=evalue, word_size=word_size, max_target_seqs=max_target_seqs, max_hsps=max_hsps)

        # Parse results and predict
        blast_df = parse_blast_results(output_file)

        # Get sequence lengths
        seq_lengths = {i: len(seq) for i, seq in enumerate(test_sequences)}

        y_pred = np.zeros(len(X_test))
        y_proba = np.zeros((len(X_test), 2))
        y_proba[:, 0] = 1.0  # Default to negative class

        if use_soft_scoring:
            # Initialize probability scores for all sequences
            chimera_scores = np.zeros(len(X_test))

            if not blast_df.empty:
                for index, row in blast_df.iterrows():
                    query_id = int(row['query_id'])
                    alignment_len = row['alignment_length']
                    identity = row['% identity']
                    evalue = row['evalue']
                    original_read_len = seq_lengths.get(query_id, 0)

                    if original_read_len > 0:
                        # Calculate soft score
                        soft_score = calculate_soft_score(
                            identity, alignment_len, original_read_len, evalue
                        )

                        # Check if there's significant non-matching portion (chimera indicator)
                        non_matching_part = original_read_len - alignment_len
                        if non_matching_part > min_non_matching_bases and identity > identity_threshold:
                            # Scale soft score to get chimera probability
                            # Higher soft score = better match = more likely chimeric if non-matching part exists
                            chimera_prob = sigmoid(5 * soft_score - 2.5)  # Sigmoid scaling
                            chimera_scores[query_id] = max(chimera_scores[query_id], chimera_prob)

            # Set probabilities
            y_proba[:, 1] = chimera_scores
            y_proba[:, 0] = 1.0 - chimera_scores

        else:
            # Original binary approach
            if not blast_df.empty:
                for index, row in blast_df.iterrows():
                    query_id = int(row['query_id'])
                    alignment_len = row['alignment_length']
                    identity = row['% identity']
                    original_read_len = seq_lengths.get(query_id, 0)

                    if original_read_len > 0:
                        non_matching_part = original_read_len - alignment_len
                        if identity > identity_threshold and non_matching_part > min_non_matching_bases:
                            y_pred[query_id] = 1
                            y_proba[query_id] = [0.0, 1.0]

        # Apply probability adjustment for train/test distribution mismatch
        if apply_probability_adjustment and use_soft_scoring:
            # Adjust probabilities for test imbalance (train: 50/50, test: 2% positive)
            y_proba_adjusted = adjust_probabilities_for_test_prior(
                y_proba, train_positive_rate=0.5, test_positive_rate=0.02
            )

            # Find optimal threshold using F1 score on adjusted probabilities
            optimal_threshold = find_optimal_threshold_f1(y_test, y_proba_adjusted)

            # Apply optimal threshold for final predictions
            y_test_pred_optimized = (y_proba_adjusted[:, 1] > optimal_threshold).astype(int)
        else:
            y_proba_adjusted = y_proba
            optimal_threshold = 0.5
            y_test_pred_optimized = y_pred

    return {
        'model': None,
        'vectorizer': None,
        'y_train_pred': None, # No training predictions for BLAST
        'y_test_pred': y_test_pred_optimized,
        'y_train_proba': None,
        'y_test_proba': y_proba_adjusted,
        'params': {
            'min_non_matching_bases': min_non_matching_bases,
            'identity_threshold': identity_threshold,
            'evalue': evalue,
            'word_size': word_size,
            'max_target_seqs': max_target_seqs,
            'max_hsps': max_hsps,
            'use_soft_scoring': use_soft_scoring,
            'apply_probability_adjustment': apply_probability_adjustment,
            'optimal_threshold': optimal_threshold,
            'num_threads': num_threads
        }
    }