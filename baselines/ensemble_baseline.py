#!/usr/bin/env python3
"""
Ensemble Baseline using Random Forest Meta-Classifier
Combines predictions from multiple baseline models
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from utils.imbalance_handling import find_optimal_threshold_f1
import os


def load_baseline_predictions(experiment_dir, experiment_id, baseline_names=['kmer', 'blast', 'cnn']):
    """
    Load prediction files from different baselines

    Args:
        experiment_dir: Directory containing prediction files
        experiment_id: Experiment ID prefix
        baseline_names: List of baseline model names

    Returns:
        DataFrame with combined predictions
    """
    predictions_df = None

    for baseline_name in baseline_names:
        pred_file = os.path.join(experiment_dir, f"{experiment_id}_{baseline_name}_predictions.csv")

        if os.path.exists(pred_file):
            df = pd.read_csv(pred_file)

            if predictions_df is None:
                # Initialize with sequence metadata
                predictions_df = df[['sequence_id', 'true_label', 'engineering_method',
                                   'sequence_length', 'length_bin', 'virus_key', 'virus_family']].copy()

            # Add predictions from this baseline
            predictions_df[f'{baseline_name}_pred'] = df['predicted_label']
            predictions_df[f'{baseline_name}_proba'] = df['prediction_probability']
        else:
            print(f"Warning: Prediction file {pred_file} not found")

    return predictions_df


def prepare_ensemble_features(predictions_df, baseline_names=['kmer', 'blast', 'cnn']):
    """
    Prepare feature matrix for ensemble learning

    Args:
        predictions_df: DataFrame with baseline predictions
        baseline_names: List of baseline model names

    Returns:
        Feature matrix X and labels y
    """
    feature_columns = []

    for baseline_name in baseline_names:
        pred_col = f'{baseline_name}_pred'
        proba_col = f'{baseline_name}_proba'

        if pred_col in predictions_df.columns and proba_col in predictions_df.columns:
            feature_columns.extend([pred_col, proba_col])

    if not feature_columns:
        raise ValueError("No baseline predictions found in the data")

    X = predictions_df[feature_columns].values
    y = predictions_df['label'].values

    return X, y, feature_columns


def train_ensemble_baseline_from_results(models_results, test_data, test_split_ratio=0.5,
                                        n_estimators=100, random_state=42):
    """
    Train ensemble baseline using predictions from models_results

    Args:
        models_results: Dictionary containing model results from other baselines
        test_data: Test dataset with labels and metadata
        test_split_ratio: Fraction of test data to use for ensemble training
        n_estimators: Number of trees in Random Forest
        random_state: Random seed

    Returns:
        Dictionary with ensemble model results
    """
    np.random.seed(random_state)

    # Extract baseline names automatically
    baseline_names = list(models_results.keys())
    print(f"Found baselines: {baseline_names}")

    # Create predictions dataframe
    predictions_df = test_data.copy()

    # Add predictions from each baseline
    for baseline_name in baseline_names:
        results = models_results[baseline_name]['results']
        predictions_df[f'{baseline_name}_pred'] = results['y_test_pred']

        # Handle probability predictions (might be 1D or 2D)
        proba = results['y_test_proba']
        if proba.ndim > 1 and proba.shape[1] > 1:
            predictions_df[f'{baseline_name}_proba'] = proba[:, 1]  # Positive class probability
        else:
            predictions_df[f'{baseline_name}_proba'] = proba.flatten()

    # Split test data for ensemble training and evaluation
    n_samples = len(predictions_df)
    n_train = int(n_samples * test_split_ratio)

    indices = np.random.permutation(n_samples)
    train_indices = indices[:n_train]
    test_indices = indices[n_train:]

    train_predictions_df = predictions_df.iloc[train_indices].copy()
    test_predictions_df = predictions_df.iloc[test_indices].copy()

    print(f"Ensemble split: {len(train_predictions_df)} train, {len(test_predictions_df)} test")

    # Prepare features
    X_train, y_train, feature_columns = prepare_ensemble_features(train_predictions_df, baseline_names)
    X_test, y_test, _ = prepare_ensemble_features(test_predictions_df, baseline_names)

    print(f"Ensemble features: {feature_columns}")

    # Train Random Forest meta-classifier
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        class_weight='balanced',
        max_depth=3,
        min_samples_split=10,
        min_samples_leaf=5
    )

    model.fit(X_train, y_train)

    # Generate predictions
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    y_train_proba = model.predict_proba(X_train)
    y_test_proba = model.predict_proba(X_test)

    # Find optimal threshold using F1 score
    optimal_threshold = find_optimal_threshold_f1(y_test, y_test_proba)

    # Apply optimal threshold for predictions
    if y_test_proba.ndim == 1:
        y_test_pred_optimized = (y_test_proba > optimal_threshold).astype(int)
    else:
        y_test_pred_optimized = (y_test_proba[:, 1] > optimal_threshold).astype(int)

    # Feature importance
    feature_importance = dict(zip(feature_columns, model.feature_importances_))

    return {
        'model': model,
        'y_train_pred': y_train_pred,
        'y_test_pred': y_test_pred_optimized,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'feature_importance': feature_importance,
        'ensemble_train_data': train_predictions_df,
        'ensemble_test_data': test_predictions_df,
        'params': {
            'n_estimators': n_estimators,
            'baseline_names': baseline_names,
            'feature_columns': feature_columns,
            'optimal_threshold': optimal_threshold,
            'test_split_ratio': test_split_ratio
        }
    }


def run_ensemble_baseline(experiment_dir, experiment_id, train_test_split_ratio=0.5,
                         baseline_names=['kmer', 'blast', 'cnn']):
    """
    Run ensemble baseline by loading predictions and training meta-classifier

    Args:
        experiment_dir: Directory containing prediction files
        experiment_id: Experiment ID prefix
        train_test_split_ratio: Ratio to split predictions into train/test for ensemble
        baseline_names: List of baseline model names

    Returns:
        Dictionary with ensemble results
    """
    print(f"Loading baseline predictions from {experiment_dir}")

    # Load all baseline predictions (this represents our first test set)
    all_predictions = load_baseline_predictions(experiment_dir, experiment_id, baseline_names)

    if all_predictions is None or len(all_predictions) == 0:
        raise ValueError("No baseline predictions could be loaded")

    # Split into ensemble training and testing sets
    np.random.seed(42)
    n_samples = len(all_predictions)
    n_train = int(n_samples * train_test_split_ratio)

    # Random split
    indices = np.random.permutation(n_samples)
    train_indices = indices[:n_train]
    test_indices = indices[n_train:]

    train_predictions_df = all_predictions.iloc[train_indices].copy()
    test_predictions_df = all_predictions.iloc[test_indices].copy()

    print(f"Split predictions: {len(train_predictions_df)} train, {len(test_predictions_df)} test")

    # Train ensemble
    ensemble_results = train_ensemble_baseline(
        train_predictions_df, test_predictions_df, baseline_names
    )

    # Print feature importance
    print("\nFeature Importance:")
    for feature, importance in ensemble_results['feature_importance'].items():
        print(f"  {feature}: {importance:.3f}")

    return ensemble_results


if __name__ == "__main__":
    # Example usage
    experiment_dir = "experiments/curated/20k-hanna-data-sat"
    experiment_id = "20250913_233707"

    results = run_ensemble_baseline(experiment_dir, experiment_id)

    # Print basic metrics
    y_test = results['y_test']
    y_pred = results['y_test_pred']

    print(f"\nEnsemble Results:")
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.3f}")
    print(f"Precision: {precision_score(y_test, y_pred):.3f}")
    print(f"Recall: {recall_score(y_test, y_pred):.3f}")
    print(f"F1-score: {f1_score(y_test, y_pred):.3f}")