#!/usr/bin/env python3
"""
Simple Metrics Utility
Calculates precision, recall and naive baseline performance
"""

import numpy as np
from sklearn.metrics import precision_score, recall_score, precision_recall_curve

def calculate_precision_recall(y_true, y_pred, pos_label=1, test_pos_label_fraction=None):
    """
    Calculate precision and recall for binary classification, with an option to filter
    the dataset to have a given fraction of positive examples.
    
    Args:
        y_true: True labels
        y_pred: Predicted labels
        pos_label: Positive class label
        test_pos_label_fraction: Desired fraction of positive examples in the test set
    
    Returns:
        Dictionary with precision and recall
    """
    if test_pos_label_fraction is not None:
        y_true, y_pred = filter_by_positive_fraction(y_true, y_pred, test_pos_label_fraction, pos_label)
    
    precision = precision_score(y_true, y_pred, pos_label=pos_label, zero_division=0)
    recall = recall_score(y_true, y_pred, pos_label=pos_label, zero_division=0)
    
    return {
        'precision': precision,
        'recall': recall
    }

def naive_fraction_strategy(y_true, positive_fraction=None, test_pos_label_fraction=None):
    """
    Generate naive predictions by randomly predicting positive class
    with probability equal to the positive fraction, with an option to filter
    the dataset to have a given fraction of positive examples.
    
    Args:
        y_true: True labels 
        positive_fraction: Fraction to predict as positive (if None, uses actual fraction)
        test_pos_label_fraction: Desired fraction of positive examples in the test set
    
    Returns:
        Array of naive predictions and the positive fraction used
    """
    if test_pos_label_fraction is not None:
        y_true, _ = filter_by_positive_fraction(y_true, y_true, test_pos_label_fraction, 1)
    
    if positive_fraction is None:
        positive_fraction = np.mean(y_true)
    
    # Generate random predictions with the specified positive fraction
    np.random.seed(42)  # For reproducible results
    naive_pred = np.random.binomial(1, positive_fraction, size=len(y_true))
    
    return naive_pred, positive_fraction

def compare_to_naive_fraction(y_true, y_pred, test_pos_label_fraction=None):
    """
    Compare model predictions to naive fraction strategy, with an option to filter
    the dataset to have a given fraction of positive examples.
    
    Args:
        y_true: True labels
        y_pred: Model predictions  
        test_pos_label_fraction: Desired fraction of positive examples in the test set
    
    Returns:
        Dictionary with model and naive baseline metrics
    """
    if test_pos_label_fraction is not None:
        y_true, y_pred = filter_by_positive_fraction(y_true, y_pred, test_pos_label_fraction, 1)
    
    # Calculate model metrics
    model_metrics = calculate_precision_recall(y_true, y_pred)
    
    # Generate naive predictions
    naive_pred, fraction_used = naive_fraction_strategy(y_true, test_pos_label_fraction)
    naive_metrics = calculate_precision_recall(y_true, naive_pred)
    
    return {
        'model': model_metrics,
        'naive': naive_metrics,
        'positive_fraction': fraction_used
    }

def compare_with_calibration(y_train_true, y_train_proba, y_test_true, y_test_proba, 
                           test_pos_label_fraction=None):
    """
    Compare model predictions to naive baseline with threshold calibration.
    
    Args:
        y_train_true: True labels for training/validation (for threshold calibration)
        y_train_proba: Predicted probabilities for training/validation
        y_test_true: True labels for test set
        y_test_proba: Predicted probabilities for test set
        test_pos_label_fraction: Desired fraction of positive examples in test evaluation
    
    Returns:
        Dictionary with calibrated model and naive baseline metrics
    """
    # Calibrate threshold on training data considering test positive fraction
    calibration_results = calibrate_threshold(
        y_train_true, y_train_proba, test_pos_label_fraction
    )
    optimal_threshold = calibration_results['optimal_threshold']
    
    print(f"Optimal threshold: {optimal_threshold:.3f} (F1: {calibration_results['optimal_f1']:.3f})")
    
    # Apply calibrated threshold to test probabilities
    y_test_pred_calibrated = apply_threshold(y_test_proba, optimal_threshold)
    
    # Compare with naive baseline
    comparison = compare_to_naive_fraction(
        y_test_true, y_test_pred_calibrated, test_pos_label_fraction
    )
    
    # Add calibration info to results
    comparison['calibration'] = calibration_results
    
    return comparison

def filter_by_positive_fraction(y_true, y_pred, desired_fraction, pos_label):
    """
    Filter the dataset to have a given fraction of positive examples, maximizing the test set size
    without repeating datapoints.
    
    Args:
        y_true: True labels
        y_pred: Predicted labels
        desired_fraction: Desired fraction of positive examples
        pos_label: Positive class label
    
    Returns:
        Filtered true and predicted labels
    """
    pos_indices = np.where(y_true == pos_label)[0]
    neg_indices = np.where(y_true != pos_label)[0]
    
    num_pos = min(int(desired_fraction * len(y_true)), len(pos_indices))
    num_neg = min(len(y_true) - num_pos, len(neg_indices))
    
    selected_pos_indices = np.random.choice(pos_indices, num_pos, replace=False)
    selected_neg_indices = np.random.choice(neg_indices, num_neg, replace=False)
    
    selected_indices = np.concatenate([selected_pos_indices, selected_neg_indices])
    np.random.shuffle(selected_indices)
    
    return y_true[selected_indices], y_pred[selected_indices]

def calibrate_threshold(y_true, y_proba, test_pos_label_fraction=None, pos_label=1):
    """
    Calibrate optimal threshold using precision-recall curve, with option to filter
    the dataset to match test positive fraction.
    
    Args:
        y_true: True labels for calibration (typically validation set)
        y_proba: Predicted probabilities for positive class
        test_pos_label_fraction: Desired fraction of positive examples in test set
        pos_label: Positive class label
    
    Returns:
        Dictionary with optimal_threshold, f1_score, and other metrics
    """
    # Filter to match test positive fraction if specified
    if test_pos_label_fraction is not None:
        y_true_filtered, y_proba_filtered = filter_by_positive_fraction(
            y_true, y_proba, test_pos_label_fraction, pos_label
        )
    else:
        y_true_filtered, y_proba_filtered = y_true, y_proba
    
    # Calculate precision-recall curve
    precision, recall, thresholds = precision_recall_curve(y_true_filtered, y_proba_filtered)
    
    # Calculate F1 scores for each threshold
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
    
    # Find optimal threshold
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else 0.5
    optimal_f1 = f1_scores[optimal_idx]
    
    return {
        'optimal_threshold': optimal_threshold,
        'optimal_f1': optimal_f1,
        'precision': precision[optimal_idx],
        'recall': recall[optimal_idx]
    }

def apply_threshold(y_proba, threshold):
    """
    Apply threshold to probabilities to get binary predictions.
    
    Args:
        y_proba: Predicted probabilities for positive class
        threshold: Classification threshold
    
    Returns:
        Binary predictions
    """
    return (y_proba >= threshold).astype(int)