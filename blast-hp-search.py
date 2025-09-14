#!/usr/bin/env python3
"""
BLAST Hyperparameter Search Script
Performs random search on BLAST hyperparameters for virus engineering detection
"""

import pandas as pd
import numpy as np
import json
import os
from utils.virus_data_processor import create_virus_dataset
from utils.metrics import evaluate_model_performance
from baselines.blast_baseline import run_blast_baseline

def generate_random_hyperparameters():
    """Generate random hyperparameters for BLAST."""
    return {
        'min_non_matching_bases': np.random.randint(10, 100),  # 10-99
        'identity_threshold': np.random.uniform(75, 99),  # 75-99%
        'evalue': np.random.choice([0.001, 0.01, 0.1, 1, 10, 100, 1000]),  # Log scale
        'word_size': np.random.randint(4, 15),  # 4-14
        'max_target_seqs': np.random.randint(1, 11),  # 1-10
        'max_hsps': np.random.randint(1, 6),  # 1-5
        'use_soft_scoring': np.random.choice([True, False]),  # Boolean choice
        'apply_probability_adjustment': np.random.choice([True, False])  # Boolean choice
    }

def evaluate_hyperparameters(X_train, X_test, y_train, y_test, hyperparams):
    """Evaluate a set of hyperparameters."""
    blast_results = run_blast_baseline(X_train, X_test, y_train, y_test, **hyperparams)
    evaluation = evaluate_model_performance(
        y_test.values, blast_results['y_test_pred'], None, blast_results['y_test_proba']
    )
    return evaluation['overall']['model']['f1'], evaluation

def run_hyperparameter_search(n_trials=100, data_size=1000):
    """Run random hyperparameter search for BLAST."""

    print("=" * 60)
    print("BLAST HYPERPARAMETER SEARCH")
    print("=" * 60)

    # Create single dataset for evaluation
    print(f"\nCreating dataset with {data_size} samples...")

    # Create training data (for BLAST database creation)
    train_data = create_virus_dataset(
        n_samples=int(data_size * 0.8),
        engineered_data_path="data/processed_engineered_virus.csv",
        original_data_path="data/processed_original_virus.csv",
        engineering_fraction=0.5,
        split_filter='train',
        split_type='random',
        random_substitution_only=True
    )

    # Create test data (for evaluation)
    test_data = create_virus_dataset(
        n_samples=int(data_size * 0.2),
        engineering_fraction=0.02,  # Use higher fraction for better testing
        engineered_data_path="data/processed_engineered_virus.csv",
        original_data_path="data/processed_original_virus.csv",
        split_filter='test',
        split_type='random',
        random_substitution_only=True
    )

    X_train, y_train = train_data.drop('label', axis=1), train_data['label']
    X_test, y_test = test_data.drop('label', axis=1), test_data['label']

    print(f"Train: {len(train_data)} samples (eng. frac: {y_train.mean():.3f})")
    print(f"Test: {len(test_data)} samples (eng. frac: {y_test.mean():.3f})")

    # Random search
    best_score = 0
    best_params = None
    best_evaluation = None
    results = []

    print(f"\nStarting random search with {n_trials} trials...")

    for trial in range(n_trials):
        try:
            hyperparams = generate_random_hyperparameters()
            f1_score, evaluation = evaluate_hyperparameters(X_train, X_test, y_train, y_test, hyperparams)

            results.append({
                'trial': trial,
                'f1_score': f1_score,
                'hyperparams': hyperparams,
                'metrics': evaluation['overall']['model']
            })

            if f1_score > best_score:
                best_score = f1_score
                best_params = hyperparams.copy()
                best_evaluation = evaluation
                print(f"Trial {trial+1:3d}/{n_trials}: New best F1={f1_score:.4f}")
            else:
                print(f"Trial {trial+1:3d}/{n_trials}: F1={f1_score:.4f}")

        except Exception as e:
            print(f"Trial {trial+1:3d}/{n_trials}: Error - {str(e)}")
            import traceback
            traceback.print_exc()
            continue

    return best_params, best_evaluation, results

def save_results(best_params, best_evaluation, results, experiment_dir):
    """Save hyperparameter search results to experiment directory."""
    os.makedirs(experiment_dir, exist_ok=True)

    # Save best hyperparameters
    with open(os.path.join(experiment_dir, 'best_hyperparameters.json'), 'w') as f:
        json.dump(best_params, f, indent=2, default=str)

    # Save best evaluation metrics
    with open(os.path.join(experiment_dir, 'best_evaluation.json'), 'w') as f:
        json.dump(best_evaluation, f, indent=2, default=str)

    # Save all trial results
    with open(os.path.join(experiment_dir, 'all_trials.json'), 'w') as f:
        json.dump(results, f, indent=2, default=str)

    # Create summary report
    if best_evaluation is None:
        summary = {
            'best_f1_score': 0.0,
            'best_precision': 0.0,
            'best_recall': 0.0,
            'best_accuracy': 0.0,
            'best_hyperparameters': None,
            'total_trials': len(results),
            'successful_trials': 0,
            'status': 'All trials failed'
        }
    else:
        summary = {
            'best_f1_score': best_evaluation['overall']['model']['f1'],
            'best_precision': best_evaluation['overall']['model']['precision'],
            'best_recall': best_evaluation['overall']['model']['recall'],
            'best_accuracy': best_evaluation['overall']['model']['accuracy'],
            'best_hyperparameters': best_params,
            'total_trials': len(results),
            'successful_trials': len([r for r in results if 'f1_score' in r and r['f1_score'] is not None])
        }

    with open(os.path.join(experiment_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {experiment_dir}")
    if best_evaluation is not None:
        print(f"Best F1 Score: {summary['best_f1_score']:.4f}")
        print(f"Best Hyperparameters: {best_params}")
    else:
        print("All trials failed - no successful runs found")

if __name__ == "__main__":
    # Set random seeds for reproducibility
    import torch
    import random

    SEED = 42
    torch.manual_seed(SEED)
    torch.cuda.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    np.random.seed(SEED)
    random.seed(SEED)

    # Make CUDA operations deterministic
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # Set Python hash seed
    os.environ['PYTHONHASHSEED'] = str(SEED)

    # Run hyperparameter search
    best_params, best_evaluation, results = run_hyperparameter_search(
        n_trials=50,  # Full hyperparameter search
        data_size=20000
    )

    # Save results to experiments directory
    experiment_dir = "experiments/blast_hyperparameter_search"
    save_results(best_params, best_evaluation, results, experiment_dir)