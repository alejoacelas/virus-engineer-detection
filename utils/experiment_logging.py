#!/usr/bin/env python3
"""
Experiment Logging System
Tracks training runs, parameters, and results
"""

import json
import pandas as pd
from datetime import datetime
from pathlib import Path

def setup_experiment_dirs(base_dir="experiments"):
    """Setup experiment directories"""
    base_path = Path(base_dir)
    logs_dir = base_path / "logs"
    predictions_dir = base_path / "predictions"

    logs_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir.mkdir(parents=True, exist_ok=True)

    return logs_dir, predictions_dir

def generate_experiment_id():
    """Generate unique experiment ID"""
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def log_experiment(experiment_id, model_name, results, data_params, model_params, user_notes=""):
    """Log experiment results to JSON file"""
    logs_dir, _ = setup_experiment_dirs()

    log_entry = {
        'experiment_id': experiment_id,
        'timestamp': datetime.now().isoformat(),
        'model_name': model_name,
        'data_params': data_params,
        'model_params': model_params,
        'results': results,
        'user_notes': user_notes
    }

    log_file = logs_dir / f"{experiment_id}_{model_name}.json"
    with open(log_file, 'w') as f:
        json.dump(log_entry, f, indent=2, default=str)

    return log_file

def save_predictions(experiment_id, model_name, test_data, y_true, y_pred, y_proba=None):
    """Save individual predictions to CSV file"""
    _, predictions_dir = setup_experiment_dirs()

    # Prepare prediction data
    pred_data = {
        'sequence_id': test_data['sequence_id'] if 'sequence_id' in test_data else range(len(y_true)),
        'true_label': y_true,
        'predicted_label': y_pred
    }

    if y_proba is not None:
        pred_data['prediction_probability'] = y_proba[:, 1] if y_proba.ndim > 1 else y_proba

    # Add metadata if available
    for col in ['engineering_method', 'sequence_length', 'length_bin', 'virus_key', 'virus_family']:
        if col in test_data:
            pred_data[col] = test_data[col]

    pred_df = pd.DataFrame(pred_data)
    pred_file = predictions_dir / f"{experiment_id}_{model_name}_predictions.csv"
    pred_df.to_csv(pred_file, index=False)

    return pred_file

def prompt_user_notes(interactive=True):
    """Prompt user for experiment notes"""
    if not interactive:
        return "Automated test run"

    print("\nExperiment Notes:")
    print("Enter any notes about this experiment (changes made, hypothesis, etc.)")
    print("Press Enter twice to finish:")

    notes = []
    while True:
        try:
            line = input()
            if line == "":
                break
            notes.append(line)
        except EOFError:
            break

    return "\n".join(notes)