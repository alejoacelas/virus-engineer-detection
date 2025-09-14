#!/usr/bin/env python3
"""
Unified Training and Testing Script for Virus Engineering Detection
Trains all baselines with consistent data splits and comprehensive evaluation
"""

import pandas as pd
# from utils.data_generator import create_engineering_dataset
from utils.virus_data_processor import create_virus_dataset
from utils.metrics import evaluate_model_performance, print_metrics_comparison
from utils.experiment_logging import log_experiment, save_predictions, prompt_user_notes, generate_experiment_id
from baselines.kmer_baseline import train_kmer_baseline
from baselines.cnn_baseline import train_cnn_baseline
from baselines.blast_baseline import run_blast_baseline
from utils.virus_data_processor import load_ged_test_dataset, load_viral_vectors_dataset

def run_experiment(engineered_data_path="data/processed_engineered_virus.csv",
                  original_data_path="data/processed_original_virus.csv",
                  n_samples=5000,
                  engineering_fraction=0.02,
                  split_type='family',
                  random_substitution_only=False
                  ):
    """Run complete training and testing experiment"""

    print("=" * 60)
    print("VIRUS ENGINEERING DETECTION - TRAINING & TESTING")
    print("=" * 60)

    # Generate experiment ID
    experiment_id = generate_experiment_id()
    print(f"Experiment ID: {experiment_id}")

    print("\n1. Creating training dataset")
    train_data = create_virus_dataset(
        n_samples=int(n_samples * 0.8),  # 80% for training
        engineering_fraction=0.5,
        engineered_data_path=engineered_data_path,
        original_data_path=original_data_path,
        split_filter='train',
        split_type=split_type,
        random_substitution_only=random_substitution_only
    )

    print("\n2. Creating test dataset")
    test_data = create_virus_dataset(
        n_samples=int(n_samples * 0.2),  # 20% for testing
        engineering_fraction=engineering_fraction,
        engineered_data_path=engineered_data_path,
        original_data_path=original_data_path,
        split_filter='test',
        split_type=split_type,
        random_substitution_only=random_substitution_only
    )
    # test_data = load_ged_test_dataset("data/GED_test.csv", n_samples=int(n_samples * 0.2))
    # test_data = load_viral_vectors_dataset("data/GED_test.csv", n_samples=int(n_samples * 0.2))

    print("\n")
    print(f"Training engineering fraction: {train_data['label'].mean():.3f}")
    print(f"Test engineering fraction: {test_data['label'].mean():.3f}")

    X_train, y_train = train_data.drop('label', axis=1), train_data['label']
    X_test, y_test = test_data.drop('label', axis=1), test_data['label']

    print(f"Train: {len(train_data)} samples, Test: {len(test_data)} samples")

    # Data parameters for logging
    data_params = {
        'n_samples': n_samples,
        'engineering_fraction_train': 0.5,
        'engineering_fraction_test': engineering_fraction,
        'split_type': split_type,
        'train_size': len(train_data),
        'test_size': len(test_data),
        'train_methods': train_data['engineering_method'].value_counts().to_dict(),
        'test_methods': test_data['engineering_method'].value_counts().to_dict()
    }

    # Train and evaluate models
    models_results = {}

    print("\nTraining K-mer Logistic Regression...")
    kmer_results = train_kmer_baseline(X_train, X_test, y_train, y_test, k=6, max_features=1000)
    kmer_eval = evaluate_model_performance(
        y_test.values, kmer_results['y_test_pred'], test_data, kmer_results['y_test_proba']
    )
    models_results['kmer'] = {'results': kmer_results, 'evaluation': kmer_eval}

    print("K-mer Results:")
    print_metrics_comparison(kmer_eval['overall'])

    print("\nRunning BLAST-based Baseline...")
    blast_results = run_blast_baseline(X_train, X_test, y_train, y_test, num_threads=18)
    blast_eval = evaluate_model_performance(
        y_test.values, blast_results['y_test_pred'], test_data, blast_results['y_test_proba']
    )
    models_results['blast'] = {'results': blast_results, 'evaluation': blast_eval}

    print("BLAST Results:")
    print_metrics_comparison(blast_eval['overall'])

    # CNN baseline
    print("\nTraining CNN Baseline...")
    cnn_results = train_cnn_baseline(X_train, X_test, y_train, y_test, epochs=10, use_improved=False)
    cnn_eval = evaluate_model_performance(
        y_test.values, cnn_results['y_test_pred'], test_data, cnn_results['y_test_proba']
    )
    models_results['cnn'] = {'results': cnn_results, 'evaluation': cnn_eval}

    print("CNN Results:")
    print_metrics_comparison(cnn_eval['overall'])

    # Improved CNN baseline
    # print("\n5. Training Improved CNN Baseline...")
    # improved_cnn_results = train_cnn_baseline(X_train, X_test, y_train, y_test, epochs=10, use_improved=True)
    # improved_cnn_eval = evaluate_model_performance(
    #     y_test.values, improved_cnn_results['y_test_pred'], test_data, improved_cnn_results['y_test_proba']
    # )
    # models_results['improved_cnn'] = {'results': improved_cnn_results, 'evaluation': improved_cnn_eval}

    # print("Improved CNN Results:")
    # print_metrics_comparison(improved_cnn_eval['overall'])

    # Summary comparison
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    for model_name, model_data in models_results.items():
        eval_data = model_data['evaluation']['overall']
        model_metrics = eval_data['model']
        print(f"{model_name:15} - Prec: {model_metrics['precision']:.3f}, "
              f"Rec: {model_metrics['recall']:.3f}, "
              f"Acc: {model_metrics['accuracy']:.3f}")

    # User notes
    user_notes = prompt_user_notes(interactive=False)  # Set to True for interactive mode

    # Log experiments
    print("\nLogging experiments...")
    for model_name, model_data in models_results.items():
        log_experiment(
            experiment_id, model_name,
            model_data['evaluation'], data_params,
            model_data['results']['params'], user_notes
        )

        # Save predictions
        save_predictions(
            experiment_id, model_name, test_data,
            y_test.values, model_data['results']['y_test_pred'],
            model_data['results']['y_test_proba']
        )

    print(f"Experiment {experiment_id} completed and logged!")
    return experiment_id, models_results

if __name__ == "__main__":
    # Set random seeds for reproducibility
    import torch
    import numpy as np
    import random
    import os

    # Set seeds
    SEED = 42
    torch.manual_seed(SEED)
    torch.cuda.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)  # For multi-GPU
    np.random.seed(SEED)
    random.seed(SEED)

    # Make CUDA operations deterministic
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # Set Python hash seed
    os.environ['PYTHONHASHSEED'] = str(SEED)
    
    experiment_id, results = run_experiment(
        engineered_data_path="data/processed_engineered_virus.csv",
        original_data_path="data/processed_original_virus.csv",
        n_samples=20_000,
        engineering_fraction=0.02,
        split_type='random'
    )
