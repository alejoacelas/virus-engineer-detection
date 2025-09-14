#!/usr/bin/env python3
"""
Imbalance Handling Utilities
For handling train-test imbalance scenarios (e.g., balanced train -> 2% positive test)
"""

import numpy as np


def adjust_probabilities_for_test_prior(y_proba, train_positive_rate, test_positive_rate):
    """
    Adjust predicted probabilities based on different train/test prior distributions.

    This method implements Bayes' theorem to correct for prior probability shifts:
    P(y=1|x, test) = P(y=1|x, train) * (P(y=1, test) / P(y=1, train)) / normalization

    Args:
        y_proba: Predicted probabilities from model, shape (n_samples, 2) or (n_samples,)
        train_positive_rate: Fraction of positive labels in training data
        test_positive_rate: Expected fraction of positive labels in test data

    Returns:
        numpy.ndarray: Adjusted probabilities with same shape as input
    """
    # Handle both 2D (n_samples, 2) and 1D (n_samples,) probability formats
    if y_proba.ndim == 1:
        # Convert 1D positive class probabilities to 2D format
        proba_pos = y_proba
        proba_neg = 1 - y_proba
        y_proba_2d = np.column_stack([proba_neg, proba_pos])
    else:
        y_proba_2d = y_proba.copy()

    # Calculate adjustment factor for positive class
    adjustment_factor = (test_positive_rate / (1 - test_positive_rate)) / (train_positive_rate / (1 - train_positive_rate))

    # Apply adjustment to positive class probabilities
    adjusted_proba = y_proba_2d.copy()
    adjusted_proba[:, 1] *= adjustment_factor

    # Renormalize to ensure probabilities sum to 1
    adjusted_proba = adjusted_proba / adjusted_proba.sum(axis=1, keepdims=True)

    # Return in same format as input
    if y_proba.ndim == 1:
        return adjusted_proba[:, 1]
    else:
        return adjusted_proba


def calculate_test_aware_class_weights(train_positive_rate, test_positive_rate):
    """
    Calculate class weights that account for expected test distribution.

    This enhances standard class weighting by considering the target test distribution,
    making the model more suitable for deployment on imbalanced test data.

    Args:
        train_positive_rate: Fraction of positive labels in training data
        test_positive_rate: Expected fraction of positive labels in test data

    Returns:
        numpy.ndarray: Class weights [weight_negative, weight_positive]
    """
    # Standard inverse frequency weights for training data
    train_neg_weight = train_positive_rate / (1 - train_positive_rate)

    # Adjustment factor based on test distribution
    test_adjustment = (1 - test_positive_rate) / test_positive_rate
    train_adjustment = train_positive_rate / (1 - train_positive_rate)

    # Combined weight for positive class
    pos_weight = test_adjustment * train_adjustment

    return np.array([1.0, pos_weight])


def apply_top_k_threshold(y_proba, expected_positive_rate):
    """
    Apply threshold to select top k% of samples as positive based on probability scores.

    This is useful when you know the expected positive rate in test data and want
    to select exactly that fraction of samples with highest predicted probabilities.

    Args:
        y_proba: Predicted probabilities, shape (n_samples, 2) or (n_samples,)
        expected_positive_rate: Expected fraction of positive samples (e.g., 0.02 for 2%)

    Returns:
        numpy.ndarray: Binary predictions (0 or 1)
    """
    # Extract positive class probabilities
    if y_proba.ndim == 1:
        pos_proba = y_proba
    else:
        pos_proba = y_proba[:, 1]

    # Calculate number of positive predictions needed
    n_positive = int(len(pos_proba) * expected_positive_rate)

    # Select top n_positive samples
    top_indices = np.argsort(pos_proba)[-n_positive:]

    predictions = np.zeros(len(pos_proba), dtype=int)
    predictions[top_indices] = 1

    return predictions


# Example usage functions for easy integration with baselines
def find_optimal_threshold_f1(y_true, y_proba, thresholds=None):
    """
    Find threshold that maximizes F1 score.

    Args:
        y_true: True labels
        y_proba: Predicted probabilities (1D array or 2D array with positive class in column 1)
        thresholds: Array of thresholds to test (default: 100 evenly spaced from 0 to 1)

    Returns:
        float: Optimal threshold
    """
    from sklearn.metrics import f1_score

    # Extract positive class probabilities
    if y_proba.ndim == 1:
        pos_proba = y_proba
    else:
        pos_proba = y_proba[:, 1]

    if thresholds is None:
        thresholds = np.linspace(0, 1, 101)

    best_f1 = 0
    best_threshold = 0.5

    for threshold in thresholds:
        y_pred = (pos_proba > threshold).astype(int)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        if f1 > best_f1:
            best_f1 = f1
            best_threshold = threshold

    return best_threshold

def apply_probability_adjustment(baseline_results, train_positive_rate=0.5, test_positive_rate=0.02):
    """
    Apply probability adjustment to baseline results.

    Args:
        baseline_results: Dictionary with 'y_test_proba' key
        train_positive_rate: Positive rate in training data (default 0.5 for balanced)
        test_positive_rate: Expected positive rate in test data (default 0.02)

    Returns:
        dict: Updated results with adjusted probabilities and predictions
    """
    results = baseline_results.copy()

    # Adjust probabilities
    adjusted_proba = adjust_probabilities_for_test_prior(
        baseline_results['y_test_proba'], train_positive_rate, test_positive_rate
    )

    # Update results
    results['y_test_proba_adjusted'] = adjusted_proba

    # Generate new predictions from adjusted probabilities
    if adjusted_proba.ndim == 1:
        results['y_test_pred_adjusted'] = (adjusted_proba > 0.5).astype(int)
    else:
        results['y_test_pred_adjusted'] = np.argmax(adjusted_proba, axis=1)

    return results


def get_test_aware_class_weights_for_pytorch(train_positive_rate=0.5, test_positive_rate=0.02):
    """
    Get class weights ready for PyTorch CrossEntropyLoss.

    Args:
        train_positive_rate: Positive rate in training data
        test_positive_rate: Expected positive rate in test data

    Returns:
        torch.FloatTensor: Class weights for CrossEntropyLoss
    """
    try:
        import torch
        weights = calculate_test_aware_class_weights(train_positive_rate, test_positive_rate)
        return torch.FloatTensor(weights)
    except ImportError:
        return calculate_test_aware_class_weights(train_positive_rate, test_positive_rate)