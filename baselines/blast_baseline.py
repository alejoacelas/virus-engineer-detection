#!/usr/bin/env python3
"""
BLAST-based Baseline for Chimera Detection
"""

import subprocess
import os
import pandas as pd
import tempfile
import numpy as np

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

def run_blast(query_file, db_name, output_file):
    """Run blastn and save results to a file."""
    command = [
        'blastn',
        '-query', query_file,
        '-db', db_name,
        '-out', output_file,
        '-outfmt', '6',  # Tabular output
        '-max_hsps', '1',
        '-evalue', '10',
        '-max_target_seqs', '1',
        '-word_size', '7'
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

def run_blast_baseline(X_train, X_test, y_train, y_test, min_non_matching_bases=30):
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
        run_blast(query_file, db_name, output_file)

        # Parse results and predict
        blast_df = parse_blast_results(output_file)

        # Get sequence lengths
        seq_lengths = {i: len(seq) for i, seq in enumerate(test_sequences)}

        y_pred = np.zeros(len(X_test))
        y_proba = np.zeros((len(X_test), 2))
        y_proba[:, 0] = 1.0

        if not blast_df.empty:
            for index, row in blast_df.iterrows():
                query_id = int(row['query_id'])
                alignment_len = row['alignment_length']
                identity = row['% identity']
                original_read_len = seq_lengths.get(query_id, 0)

                if original_read_len > 0:
                    non_matching_part = original_read_len - alignment_len
                    if identity > 95 and non_matching_part > min_non_matching_bases:
                        y_pred[query_id] = 1
                        y_proba[query_id] = [0.0, 1.0]


    return {
        'model': None,
        'vectorizer': None,
        'y_train_pred': None, # No training predictions for BLAST
        'y_test_pred': y_pred,
        'y_train_proba': None,
        'y_test_proba': y_proba,
        'params': {'min_non_matching_bases': min_non_matching_bases}
    }
