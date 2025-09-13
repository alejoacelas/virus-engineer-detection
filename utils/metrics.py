#!/usr/bin/env python3
"""
Enhanced Metrics Module with Subgroup Analysis
Calculates performance metrics overall and by subgroups
"""

import numpy as np
import pandas as pd
from sklearn.metrics import precision_score, recall_score, accuracy_score

def calculate_metrics(y_true, y_pred, y_proba=None):
    """Calculate basic classification metrics"""
    return {
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0)
    }

def calculate_subgroup_metrics(y_true, y_pred, subgroup_col, y_proba=None):
    """Calculate metrics by subgroup"""
    results = {}

    for group in subgroup_col.unique():
        mask = subgroup_col == group
        if mask.sum() > 0:  # Only if group has samples
            group_metrics = calculate_metrics(
                y_true[mask], y_pred[mask],
                y_proba[mask] if y_proba is not None else None
            )
            group_metrics['n_samples'] = mask.sum()
            results[group] = group_metrics

    return results

def naive_baseline_predictions(y_true, positive_fraction=None):
    """Generate naive baseline predictions"""
    if positive_fraction is None:
        positive_fraction = np.mean(y_true)

    np.random.seed(42)  # Reproducible results
    return np.random.binomial(1, positive_fraction, size=len(y_true))

def compare_with_naive(y_true, y_pred, subgroup_col=None, y_proba=None):
    """Compare model predictions with naive baseline"""
    # Model metrics
    model_metrics = calculate_metrics(y_true, y_pred, y_proba)

    # Naive baseline
    naive_pred = naive_baseline_predictions(y_true)
    naive_metrics = calculate_metrics(y_true, naive_pred)

    results = {
        'model': model_metrics,
        'naive': naive_metrics,
        'positive_fraction': np.mean(y_true)
    }

    # Subgroup analysis if requested
    if subgroup_col is not None:
        results['subgroups'] = {
            'model': calculate_subgroup_metrics(y_true, y_pred, subgroup_col, y_proba),
            'naive': calculate_subgroup_metrics(y_true, naive_pred, subgroup_col)
        }

    return results

def print_metrics_comparison(results):
    """Print formatted metrics comparison"""
    print(f"Model    - Acc: {results['model']['accuracy']:.3f}, "
          f"Prec: {results['model']['precision']:.3f}, "
          f"Rec: {results['model']['recall']:.3f}")
    print(f"Naive    - Acc: {results['naive']['accuracy']:.3f}, "
          f"Prec: {results['naive']['precision']:.3f}, "
          f"Rec: {results['naive']['recall']:.3f}")

    if 'subgroups' in results:
        print("\nSubgroup Analysis:")
        for group_name, group_metrics in results['subgroups']['model'].items():
            print(f"  {group_name}: Prec={group_metrics['precision']:.3f}, "
                  f"Rec={group_metrics['recall']:.3f}, n={group_metrics['n_samples']}")

def evaluate_model_performance(y_true, y_pred, metadata_df=None, y_proba=None):
    """Comprehensive model evaluation"""
    results = {}

    # Overall metrics
    results['overall'] = compare_with_naive(y_true, y_pred, y_proba=y_proba)

    # Subgroup analysis if metadata available
    if metadata_df is not None:
        for col in ['engineering_method', 'length_bin', 'virus_key', 'virus_family']:
            if col in metadata_df.columns:
                subgroup_results = compare_with_naive(
                    y_true, y_pred, metadata_df[col], y_proba
                )
                results[f'by_{col}'] = subgroup_results

    return results