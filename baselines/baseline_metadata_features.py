#!/usr/bin/env python3
"""
Metadata Features Baseline
Uses only the binary metadata features for classification
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, top_k_accuracy_score

def train_metadata_baseline(csv_path, use_random_forest=True):
    """Train baseline using only metadata features"""
    # Load data
    df = pd.read_csv(csv_path)
    
    # Extract metadata features (exclude sequence_id, sequence, label)
    feature_cols = [col for col in df.columns 
                   if col not in ['sequence_id', 'sequence', 'label']]
    
    X = df[feature_cols].values
    y = df['label'].values
    
    print(f"Using {X.shape[1]} metadata features")
    print(f"Features: {feature_cols[:5]}...")  # Show first 5
    
    # Encode labels
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y_encoded, test_size=0.2, random_state=42, stratify=y_encoded
    )
    
    # Train model
    if use_random_forest:
        model = RandomForestClassifier(
            n_estimators=100, 
            max_depth=10, 
            random_state=42,
            class_weight='balanced'
        )
        model_name = "Random Forest"
    else:
        model = LogisticRegression(
            max_iter=1000, 
            random_state=42,
            class_weight='balanced'
        )
        model_name = "Logistic Regression"
    
    model.fit(X_train, y_train)
    
    # Predictions
    y_train_pred = model.predict(X_train)
    y_val_pred = model.predict(X_val)
    
    # Probabilities for top-k accuracy
    y_train_proba = model.predict_proba(X_train)
    y_val_proba = model.predict_proba(X_val)
    
    # Calculate metrics
    train_acc = accuracy_score(y_train, y_train_pred)
    val_acc = accuracy_score(y_val, y_val_pred)
    
    train_top10 = top_k_accuracy_score(y_train, y_train_proba, k=10)
    val_top10 = top_k_accuracy_score(y_val, y_val_proba, k=10)
    
    print(f"\n{model_name} Results:")
    print(f"Training Accuracy: {train_acc:.3f}, Top-10: {train_top10:.3f}")
    print(f"Validation Accuracy: {val_acc:.3f}, Top-10: {val_top10:.3f}")
    
    # Feature importance for Random Forest
    if use_random_forest and hasattr(model, 'feature_importances_'):
        importance_df = pd.DataFrame({
            'feature': feature_cols,
            'importance': model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print(f"\nTop 5 most important features:")
        for i, (_, row) in enumerate(importance_df.head().iterrows()):
            print(f"{i+1}. {row['feature']}: {row['importance']:.3f}")
    
    return model, le, feature_cols

if __name__ == "__main__":
    # Try Random Forest
    rf_model, le, features = train_metadata_baseline("data/example.csv", use_random_forest=True)
    
    # Try Logistic Regression
    lr_model, _, _ = train_metadata_baseline("data/example.csv", use_random_forest=False)