#!/usr/bin/env python3
"""
Unified Training and Testing Script for Virus Engineering Detection
Trains all baselines with consistent data splits and comprehensive evaluation
"""

import pandas as pd
from utils.data_generator import create_engineering_dataset, split_dataset
from utils.metrics import evaluate_model_performance, print_metrics_comparison
from utils.experiment_logging import log_experiment, save_predictions, prompt_user_notes, generate_experiment_id
from baselines.kmer_baseline import train_kmer_baseline, train_multi_kmer_baseline
from baselines.cnn_baseline import train_cnn_baseline
from baselines.blast_baseline import run_blast_baseline

def run_experiment(csv_path="data/viral_genome_df_prelim.csv", n_samples=5000,
                  engineering_fraction=0.02, methods=['random_replacement']):
    """Run complete training and testing experiment"""

    print("=" * 60)
    print("VIRUS ENGINEERING DETECTION - TRAINING & TESTING")
    print("=" * 60)

    # Generate experiment ID
    experiment_id = generate_experiment_id()
    print(f"Experiment ID: {experiment_id}")

    # Create training dataset with 0.5 engineering fraction
    print("\n1. Creating training dataset...")
    train_data = create_engineering_dataset(
        csv_path=csv_path,
        n_samples=int(n_samples * 0.8),  # 80% for training
        engineering_fraction=0.5,
        methods=methods
    )

    # Create test dataset with 0.02 engineering fraction
    print("\n2. Creating test dataset...")
    test_data = create_engineering_dataset(
        csv_path=csv_path,
        n_samples=int(n_samples * 0.2),  # 20% for testing
        engineering_fraction=0.02,
        methods=methods
    )

    print(f"Created {len(train_data)} training samples")
    print(f"Training engineering fraction: {train_data['label'].mean():.3f}")
    print(f"Created {len(test_data)} test samples")
    print(f"Test engineering fraction: {test_data['label'].mean():.3f}")
    print(f"Methods: {train_data['engineering_method'].value_counts().to_dict()}")

    X_train, y_train = train_data.drop('label', axis=1), train_data['label']
    X_test, y_test = test_data.drop('label', axis=1), test_data['label']

    print(f"Train: {len(train_data)} samples, Test: {len(test_data)} samples")

    # Data parameters for logging
    data_params = {
        'n_samples': n_samples,
        'engineering_fraction_train': 0.5,
        'engineering_fraction_test': 0.02,
        'methods': methods,
        'train_size': len(train_data),
        'test_size': len(test_data)
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
    blast_results = run_blast_baseline(X_train, X_test, y_train, y_test)
    blast_eval = evaluate_model_performance(
        y_test.values, blast_results['y_test_pred'], test_data, blast_results['y_test_proba']
    )
    models_results['blast'] = {'results': blast_results, 'evaluation': blast_eval}

    print("BLAST Results:")
    print_metrics_comparison(blast_eval['overall'])

    # # Multi K-mer baseline
    # print("\nTraining Multi-K-mer Baseline...")
    # multi_kmer_results = train_multi_kmer_baseline(
    #     X_train, X_test, y_train, y_test, k_values=[4, 6, 8], max_features=500
    # )
    # multi_kmer_eval = evaluate_model_performance(
    #     y_test.values, multi_kmer_results['y_test_pred'], test_data, multi_kmer_results['y_test_proba']
    # )
    # models_results['multi_kmer'] = {'results': multi_kmer_results, 'evaluation': multi_kmer_eval}

    # print("Multi-K-mer Results:")
    # print_metrics_comparison(multi_kmer_eval['overall'])

    # CNN baseline
    print("\nTraining CNN Baseline...")
    cnn_results = train_cnn_baseline(X_train, X_test, y_train, y_test, epochs=10, use_improved=False)
    cnn_eval = evaluate_model_performance(
        y_test.values, cnn_results['y_test_pred'], test_data, cnn_results['y_test_proba']
    )
    models_results['cnn'] = {'results': cnn_results, 'evaluation': cnn_eval}

    print("CNN Results:")
    print_metrics_comparison(cnn_eval['overall'])

    # # Improved CNN baseline
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
    experiment_id, results = run_experiment(
        n_samples=20_000,
        engineering_fraction=0.02,
        methods=['random_replacement', 'codon_optimization', 'motif_insertion']
    )
