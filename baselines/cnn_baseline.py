#!/usr/bin/env python3
"""
CNN Baseline for DNA Sequence Classification
Standalone PyTorch implementation
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

def dna_to_onehot(sequence, max_len=200):
    """Convert DNA sequence to one-hot encoding"""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    sequence = sequence.upper()[:max_len]

    # Pad if necessary
    if len(sequence) < max_len:
        sequence += 'N' * (max_len - len(sequence))

    # One-hot encode
    onehot = np.zeros((max_len, 4))
    for i, nuc in enumerate(sequence):
        if nuc in mapping:
            onehot[i, mapping[nuc]] = 1

    return onehot

class SimpleCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 64, kernel_size=15)
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.fc1 = nn.Linear(64, 128)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(128, 2)

    def forward(self, x):
        x = torch.relu(self.conv1(x))
        x = self.pool(x).squeeze(-1)
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x)

class ImprovedCNN(nn.Module):
    def __init__(self):
        super().__init__()
        # Multi-scale branches
        self.branch1 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=3, padding=1),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )
        self.branch2 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=7, padding=3),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )
        self.branch3 = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=11, padding=5),
            nn.BatchNorm1d(32),
            nn.ReLU()
        )

        # Fusion and classification
        self.fusion = nn.Conv1d(96, 64, kernel_size=1)
        self.pool = nn.AdaptiveMaxPool1d(4)
        self.fc1 = nn.Linear(64 * 4, 128)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(128, 2)

    def forward(self, x):
        b1 = self.branch1(x)
        b2 = self.branch2(x)
        b3 = self.branch3(x)
        fused = torch.relu(self.fusion(torch.cat([b1, b2, b3], dim=1)))
        x = self.pool(fused).flatten(1)
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x)

def train_cnn_baseline(X_train, X_test, y_train, y_test, max_len=200, epochs=20,
                      batch_size=32, use_improved=False):
    """Train CNN baseline"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Convert sequences to one-hot
    train_sequences = X_train['sequence']
    test_sequences = X_test['sequence']

    X_train_onehot = np.array([dna_to_onehot(seq, max_len) for seq in train_sequences])
    X_test_onehot = np.array([dna_to_onehot(seq, max_len) for seq in test_sequences])

    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train_onehot).transpose(1, 2)
    X_test_tensor = torch.FloatTensor(X_test_onehot).transpose(1, 2)
    y_train_tensor = torch.LongTensor(y_train.values)
    y_test_tensor = torch.LongTensor(y_test.values)

    # Data loaders
    train_loader = DataLoader(
        TensorDataset(X_train_tensor, y_train_tensor),
        batch_size=batch_size, shuffle=True
    )
    test_loader = DataLoader(
        TensorDataset(X_test_tensor, y_test_tensor),
        batch_size=batch_size
    )

    # Model
    model = ImprovedCNN() if use_improved else SimpleCNN()
    model = model.to(device)

    # Class weights for imbalanced data
    class_counts = np.bincount(y_train)
    class_weights = len(y_train) / (len(class_counts) * class_counts)
    criterion = nn.CrossEntropyLoss(weight=torch.FloatTensor(class_weights).to(device))
    optimizer = optim.Adam(model.parameters())

    # Training loop
    for epoch in range(epochs):
        model.train()
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()

        if (epoch + 1) % 5 == 0:
            print(f"Epoch {epoch+1}/{epochs}")

    # Get predictions
    model.eval()
    with torch.no_grad():
        # Training predictions
        y_train_proba_list = []
        for batch_X, _ in DataLoader(TensorDataset(X_train_tensor, y_train_tensor), batch_size=batch_size):
            batch_X = batch_X.to(device)
            outputs = torch.softmax(model(batch_X), dim=1)
            y_train_proba_list.append(outputs.cpu().numpy())
        y_train_proba = np.vstack(y_train_proba_list)

        # Test predictions
        y_test_proba_list = []
        for batch_X, _ in DataLoader(TensorDataset(X_test_tensor, y_test_tensor), batch_size=batch_size):
            batch_X = batch_X.to(device)
            outputs = torch.softmax(model(batch_X), dim=1)
            y_test_proba_list.append(outputs.cpu().numpy())
        y_test_proba = np.vstack(y_test_proba_list)

    y_train_pred = np.argmax(y_train_proba, axis=1)
    y_test_pred = np.argmax(y_test_proba, axis=1)

    return {
        'model': model,
        'y_train_pred': y_train_pred,
        'y_test_pred': y_test_pred,
        'y_train_proba': y_train_proba,
        'y_test_proba': y_test_proba,
        'params': {'max_len': max_len, 'epochs': epochs, 'batch_size': batch_size, 'use_improved': use_improved}
    }