#!/usr/bin/env python3
"""
Simple CNN Baseline - PyTorch
Basic convolutional neural network for DNA sequence classification
"""

import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split

def dna_to_onehot(sequence, max_len=5000):
    """Convert DNA sequence to one-hot encoding"""
    # Map nucleotides to integers
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 0}
    
    # Truncate or pad sequence
    sequence = sequence.upper()[:max_len]
    seq_encoded = [mapping.get(nuc, 4) for nuc in sequence]
    
    # Pad with N's if needed
    if len(seq_encoded) < max_len:
        seq_encoded.extend([4] * (max_len - len(seq_encoded)))
    
    # Convert to one-hot (exclude N channel for simplicity)
    onehot = np.zeros((max_len, 4))
    for i, nuc_idx in enumerate(seq_encoded):
        if nuc_idx < 4:  # Skip N's
            onehot[i, nuc_idx] = 1
    
    return onehot

class SimpleCNN(nn.Module):
    def __init__(self, num_classes):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 64, 15)
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.fc1 = nn.Linear(64, 128)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(128, num_classes)
        
    def forward(self, x):
        x = torch.relu(self.conv1(x))
        x = self.pool(x).squeeze(-1)
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

class ImprovedCNN(nn.Module):
    def __init__(self, num_classes, seq_len=200):
        super().__init__()
        
        # Multi-scale convolutional branches
        self.conv_branch1 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=3, padding=1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(32, 32, kernel_size=3, padding=1),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )
        
        self.conv_branch2 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=7, padding=3),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(32, 32, kernel_size=7, padding=3),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )
        
        self.conv_branch3 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=11, padding=5),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(32, 32, kernel_size=11, padding=5),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )
        
        # Residual connection and fusion
        self.fusion_conv = nn.Conv1d(96, 64, kernel_size=1)  # 32*3 = 96 channels
        self.fusion_bn = nn.BatchNorm1d(64)
        
        # Final conv layer
        self.final_conv = nn.Sequential(
            nn.Conv1d(64, 64, kernel_size=3, padding=1),
            nn.BatchNorm1d(64),
            nn.ReLU()
        )
        
        # Pooling and classification
        self.pool = nn.AdaptiveMaxPool1d(4)  # Keep some spatial information
        self.fc1 = nn.Linear(64 * 4, 128)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(128, num_classes)
        
    def forward(self, x):
        # Multi-scale feature extraction
        branch1 = self.conv_branch1(x)
        branch2 = self.conv_branch2(x)
        branch3 = self.conv_branch3(x)
        
        # Concatenate multi-scale features
        fused = torch.cat([branch1, branch2, branch3], dim=1)
        
        # Fusion and final conv
        fused = torch.relu(self.fusion_bn(self.fusion_conv(fused)))
        x = self.final_conv(fused)
        
        # Classification head
        x = self.pool(x).flatten(1)
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

def train_simple_cnn(csv_path, max_len=200, epochs=50, batch_size=16):
    """Train simple CNN baseline"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load and prepare data
    df = pd.read_csv(csv_path)
    sequences = df['sequence'].values
    labels = df['label'].values.astype(int)
    
    # Convert to one-hot encoding
    print("Converting sequences to one-hot encoding...")
    X = np.array([dna_to_onehot(seq, max_len) for seq in sequences])
    X = torch.FloatTensor(X).transpose(1, 2)  # (batch, channels, length)
    
    y = torch.LongTensor(labels)
    num_classes = len(np.unique(labels))
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=labels
    )
    
    # Create data loaders
    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size)
    
    # Calculate class weights for imbalanced data
    unique_labels, label_counts = np.unique(labels, return_counts=True)
    class_weights = len(labels) / (len(unique_labels) * label_counts)
    class_weights_tensor = torch.FloatTensor(class_weights).to(device)
    
    # Create model
    model = SimpleCNN(num_classes).to(device)
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = optim.Adam(model.parameters())
    
    print(f"Training CNN on {len(X_train)} sequences, {num_classes} classes")
    
    # Training loop with progress tracking
    train_accs, val_accs = [], []
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_correct = 0
        train_total = 0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            pred = outputs.argmax(dim=1)
            train_correct += (pred == batch_y).sum().item()
            train_total += batch_y.size(0)
        
        # Validation
        model.eval()
        val_correct = 0
        val_total = 0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X)
                pred = outputs.argmax(dim=1)
                val_correct += (pred == batch_y).sum().item()
                val_total += batch_y.size(0)
        
        train_acc = train_correct / train_total
        val_acc = val_correct / val_total
        train_accs.append(train_acc)
        val_accs.append(val_acc)
        
        print(f"Epoch {epoch+1}/{epochs}, Train Acc: {train_acc:.3f}, Val Acc: {val_acc:.3f}")
    
    # Get final predictions on full datasets
    model.eval()
    with torch.no_grad():
        # Training predictions (for threshold optimization)
        train_loader_full = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size)
        y_train_proba = []
        for batch_X, batch_y in train_loader_full:
            batch_X = batch_X.to(device)
            outputs = model(batch_X)
            proba = torch.softmax(outputs, dim=1).cpu().numpy()
            y_train_proba.extend(proba)
        
        # Test predictions
        test_loader_full = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size)
        y_test_proba = []
        for batch_X, batch_y in test_loader_full:
            batch_X = batch_X.to(device)
            outputs = model(batch_X)
            proba = torch.softmax(outputs, dim=1).cpu().numpy()
            y_test_proba.extend(proba)
    
    y_train_proba = np.array(y_train_proba)
    y_test_proba = np.array(y_test_proba)
    
    # Default predictions using 0.5 threshold for accuracy calculation
    y_train_pred = np.argmax(y_train_proba, axis=1)
    y_test_pred = np.argmax(y_test_proba, axis=1)
    
    print(f"\nSimple CNN Results:")
    print(f"Final Training Accuracy: {train_accs[-1]:.3f}")
    print(f"Final Validation Accuracy: {val_accs[-1]:.3f}")
    
    return {
        'model': model,
        'y_train_true': y_train.numpy(),
        'y_train_pred': y_train_pred,
        'y_test_true': y_val.numpy(),
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'train_accs': train_accs,
        'val_accs': val_accs,
        'train_acc': train_accs[-1],
        'test_acc': val_accs[-1]
    }

def train_improved_cnn(csv_path, max_len=200, epochs=50, batch_size=16, use_improved=True):
    """Train improved multi-scale CNN for short sequences"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load and prepare data
    df = pd.read_csv(csv_path)
    sequences = df['sequence'].values
    labels = df['label'].values.astype(int)
    
    # Convert to one-hot encoding
    print("Converting sequences to one-hot encoding...")
    X = np.array([dna_to_onehot(seq, max_len) for seq in sequences])
    X = torch.FloatTensor(X).transpose(1, 2)  # (batch, channels, length)
    
    y = torch.LongTensor(labels)
    num_classes = len(np.unique(labels))
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=labels
    )
    
    # Create data loaders
    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size)
    
    # Calculate class weights for imbalanced data
    unique_labels, label_counts = np.unique(labels, return_counts=True)
    class_weights = len(labels) / (len(unique_labels) * label_counts)
    class_weights_tensor = torch.FloatTensor(class_weights).to(device)
    
    # Create model (choose between simple and improved)
    if use_improved:
        model = ImprovedCNN(num_classes, seq_len=max_len).to(device)
        model_name = "Improved CNN"
    else:
        model = SimpleCNN(num_classes).to(device)
        model_name = "Simple CNN"
    
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = optim.Adam(model.parameters())
    
    print(f"Training {model_name} on {len(X_train)} sequences, {num_classes} classes")
    
    # Training loop with progress tracking
    train_accs, val_accs = [], []
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_correct = 0
        train_total = 0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            pred = outputs.argmax(dim=1)
            train_correct += (pred == batch_y).sum().item()
            train_total += batch_y.size(0)
        
        # Validation
        model.eval()
        val_correct = 0
        val_total = 0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X)
                pred = outputs.argmax(dim=1)
                val_correct += (pred == batch_y).sum().item()
                val_total += batch_y.size(0)
        
        train_acc = train_correct / train_total
        val_acc = val_correct / val_total
        train_accs.append(train_acc)
        val_accs.append(val_acc)
        
        print(f"Epoch {epoch+1}/{epochs}, Train Acc: {train_acc:.3f}, Val Acc: {val_acc:.3f}")
    
    # Get final predictions on full datasets
    model.eval()
    with torch.no_grad():
        # Training predictions
        train_loader_full = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size)
        y_train_proba = []
        for batch_X, batch_y in train_loader_full:
            batch_X = batch_X.to(device)
            outputs = model(batch_X)
            proba = torch.softmax(outputs, dim=1).cpu().numpy()
            y_train_proba.extend(proba)
        
        # Test predictions
        test_loader_full = DataLoader(TensorDataset(X_val, y_val), batch_size=batch_size)
        y_test_proba = []
        for batch_X, batch_y in test_loader_full:
            batch_X = batch_X.to(device)
            outputs = model(batch_X)
            proba = torch.softmax(outputs, dim=1).cpu().numpy()
            y_test_proba.extend(proba)
    
    y_train_proba = np.array(y_train_proba)
    y_test_proba = np.array(y_test_proba)
    
    # Default predictions
    y_train_pred = np.argmax(y_train_proba, axis=1)
    y_test_pred = np.argmax(y_test_proba, axis=1)
    
    print(f"\n{model_name} Results:")
    print(f"Final Training Accuracy: {train_accs[-1]:.3f}")
    print(f"Final Validation Accuracy: {val_accs[-1]:.3f}")
    
    return {
        'model': model,
        'y_train_true': y_train.numpy(),
        'y_train_pred': y_train_pred,
        'y_test_true': y_val.numpy(),
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'train_accs': train_accs,
        'val_accs': val_accs,
        'train_acc': train_accs[-1],
        'test_acc': val_accs[-1]
    }

if __name__ == "__main__":
    results = train_simple_cnn("data/dummy_train.csv")
