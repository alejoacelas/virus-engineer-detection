#!/usr/bin/env python3
"""
Threshold Calibration Utility
For adjusting prediction thresholds when test data has different positive label fraction than train data
"""

import numpy as np
from sklearn.metrics import precision_recall_curve


def calibrate_threshold(y_true, y_proba, metric='f1'):
    """
    Find optimal threshold for binary classification based on precision-recall curve.

    Args:
        y_true: True binary labels (0 or 1)
        y_proba: Predicted probabilities for positive class
        metric: Optimization metric ('f1', 'precision', 'recall', or 'balanced')

    Returns:
        float: Optimal threshold value
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_proba)

    # Skip the extra point at threshold = -inf
    precision = precision[:-1]
    recall = recall[:-1]

    if metric == 'f1':
        # Compute F1 scores for each threshold
        f1_scores = 2 * precision * recall / (precision + recall + 1e-12)
        best_idx = np.nanargmax(f1_scores)
    elif metric == 'precision':
        best_idx = np.nanargmax(precision)
    elif metric == 'recall':
        best_idx = np.nanargmax(recall)
    elif metric == 'balanced':
        # Balance between precision and recall
        balanced_scores = 2 * precision * recall / (precision + recall + 1e-12)
        best_idx = np.nanargmax(balanced_scores)
    else:
        raise ValueError(f"Unknown metric: {metric}. Use 'f1', 'precision', 'recall', or 'balanced'")

    return float(thresholds[best_idx])


def get_calibration_results(y_true, y_proba, metric='f1'):
    """
    Get detailed calibration results including optimal threshold and performance metrics.

    Args:
        y_true: True binary labels (0 or 1)
        y_proba: Predicted probabilities for positive class
        metric: Optimization metric ('f1', 'precision', 'recall', or 'balanced')

    Returns:
        dict: Dictionary containing optimal threshold and performance metrics
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_proba)

    # Skip the extra point at threshold = -inf
    precision = precision[:-1]
    recall = recall[:-1]

    # Compute F1 scores for all thresholds
    f1_scores = 2 * precision * recall / (precision + recall + 1e-12)

    if metric == 'f1':
        best_idx = np.nanargmax(f1_scores)
    elif metric == 'precision':
        best_idx = np.nanargmax(precision)
    elif metric == 'recall':
        best_idx = np.nanargmax(recall)
    elif metric == 'balanced':
        best_idx = np.nanargmax(f1_scores)  # F1 is already balanced
    else:
        raise ValueError(f"Unknown metric: {metric}. Use 'f1', 'precision', 'recall', or 'balanced'")

    optimal_threshold = float(thresholds[best_idx])

    return {
        "optimal_threshold": optimal_threshold,
        "optimal_f1": float(f1_scores[best_idx]),
        "precision": float(precision[best_idx]),
        "recall": float(recall[best_idx]),
        "metric_used": metric
    }


def calibrate_threshold_with_test_fraction(y_train_true, y_train_proba, y_test_true,
                                          metric='f1', random_state=42):
    """
    Calibrate threshold using a stratified sample from training data that matches
    the positive class fraction of the test data.

    Args:
        y_train_true: True binary labels from training data
        y_train_proba: Predicted probabilities for positive class from training
        y_test_true: True binary labels from test data (to compute fraction)
        metric: Optimization metric ('f1', 'precision', 'recall', 'balanced')
        random_state: Random seed for sampling

    Returns:
        float: Optimal threshold based on matched sample
    """
    from sklearn.model_selection import train_test_split

    # Compute test positive fraction
    test_pos_fraction = np.mean(y_test_true)

    # If test fraction matches train fraction, use standard calibration
    train_pos_fraction = np.mean(y_train_true)
    if abs(train_pos_fraction - test_pos_fraction) < 1e-6:
        return calibrate_threshold(y_train_true, y_train_proba, metric)

    # Create stratified sample from training data matching test fraction
    pos_indices = np.where(y_train_true == 1)[0]
    neg_indices = np.where(y_train_true == 0)[0]

    # Calculate sample sizes
    total_sample_size = min(len(y_train_true), 1000)  # Limit sample size
    n_pos_sample = int(total_sample_size * test_pos_fraction)
    n_neg_sample = total_sample_size - n_pos_sample

    # Sample with replacement if needed
    np.random.seed(random_state)
    sampled_pos_indices = np.random.choice(pos_indices, size=min(n_pos_sample, len(pos_indices)),
                                          replace=n_pos_sample > len(pos_indices))
    sampled_neg_indices = np.random.choice(neg_indices, size=min(n_neg_sample, len(neg_indices)),
                                          replace=n_neg_sample > len(neg_indices))

    # Combine samples
    sample_indices = np.concatenate([sampled_pos_indices, sampled_neg_indices])
    y_sample_true = y_train_true[sample_indices]
    y_sample_proba = y_train_proba[sample_indices]

    # Calibrate on matched sample
    return calibrate_threshold(y_sample_true, y_sample_proba, metric)


def apply_threshold(y_proba, threshold):
    """
    Apply threshold to probabilities to get binary predictions.

    Args:
        y_proba: Predicted probabilities for positive class
        threshold: Threshold value

    Returns:
        numpy.ndarray: Binary predictions (0 or 1)
    """
    return (y_proba >= threshold).astype(int)