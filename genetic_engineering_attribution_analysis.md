# Genetic Engineering Attribution Challenge - Winning Solutions Analysis

## Competition Overview

The Genetic Engineering Attribution Challenge aimed to develop tools for identifying the lab-of-origin from genetically engineered DNA sequences. The competition had two tracks:

- **Prediction Track**: Attribute DNA samples to lab-of-origin with highest top-10 accuracy
- **Innovation Track**: Demonstrate superior models beyond raw accuracy

**Evaluation Metric**: Top-10 accuracy (whether the correct lab appears in the top 10 predictions)

## Input Data Format

The challenge used the following data structure:

### Files Required:
- `train_values.csv`: Training sequences with metadata features
- `train_labels.csv`: One-hot encoded lab labels (1314 labs total)  
- `test_values.csv`: Test sequences for inference
- `holdout_values.csv`: New sequences for final evaluation

### Input Features:
- **Primary sequence**: DNA sequence (variable length, up to 16,000+ bp)
- **Metadata features (39 total)**:
  - Bacterial resistance markers (ampicillin, chloramphenicol, kanamycin, etc.)
  - Copy number (high/low/unknown)
  - Growth conditions (strain, temperature)
  - Selectable markers (blasticidin, hygromycin, neomycin, etc.)
  - Species information (human, mouse, yeast, synthetic, etc.)

## Model Types and Architectures

### 1st Place Solution (sorrge) - Score: 0.949

**Hybrid K-mer + CNN Ensemble Approach**

#### K-mer Models:
- **K-mer sizes**: 19, 21, 23, 25, 27 base pairs
- **Implementation**: C++ optimized extraction with reverse-complement normalization
- **Features**: Maximum 18 labs per k-mer for efficiency
- **Processing**: Near disk-reading speed performance

#### CNN Ensemble (7 models):
- **Architecture types**: SimpleCNN, ResNet1D variants
- **Input sequences**: 8,000-16,384 bp (longer sequences crucial)
- **Kernel width**: 18 for main CNN layers
- **Key components**:
  - JoinNet architecture combining sequence CNN with binary features MLP
  - Concatenated max-pooling and average-pooling outputs
  - 5 direct prediction models + 2 self-supervised (NT-Xent loss)

#### Simple Naive Baseline:
- Single CNN with basic architecture (2 conv layers + dense)
- K-mer counting with scikit-learn Naive Bayes
- 8-10 mer counts as features for logistic regression

### 2nd Place Solution (fit_dna) - Score: 0.944

**BLAST-based Feature Engineering + Deep Networks**

#### Model Architecture:
- **Core approach**: BLAST subsequence matching + PCA dimensionality reduction
- **Neural networks**: 5 different architectures with varying complexity
- **Key features**:
  - BLAST similarity matrices (reduced to 1024 PCA components)
  - Multiple feature types: alignment length, gap, bitscore, mismatch
  - Ensemble of DNNs with different input combinations

#### Architecture Variants:
- **v1-v3**: Multi-input networks combining sequence embeddings + engineered features
- **v4**: Simplified categorical features only
- **v5**: Pure sequence-only model
- **Sequence processing**: Embedding layers (8 channels, 24 dimensions) + Conv1D + GlobalMaxPooling

#### Simple Naive Baseline:
- BLAST similarity search + k-NN classification
- PCA on BLAST features + simple MLP (2-3 layers)
- Direct BLAST hit counting with basic aggregation

### 3rd Place Solution (eakmail) - Score: 0.934

**Fast N-gram Kernel Methods**

#### Model Approach:
- **Core method**: Various n-gram (k-mer) kernels with Naive Bayes
- **Optimization**: Fast 4-minute training solutions
- **Features**: Kernel-based similarity with rank model merging

#### Simple Naive Baseline:
- 6-8 mer counting with basic Naive Bayes classifier
- Simple frequency-based features
- Linear SVM on k-mer counts

## Minimal Python Dependencies

### Core Requirements:
```python
# Essential packages for baseline implementation
numpy>=1.18.0
pandas>=1.0.3
scikit-learn>=0.22.1
tensorflow>=2.1.0  # or pytorch>=1.6.0

# For BLAST-based approaches
biopython>=1.78

# For k-mer analysis
collections (built-in)
itertools (built-in)

# Optional for advanced features
scipy>=1.4.1
tqdm>=4.44.1
```

### For Advanced Solutions:
```python
# 1st place requirements
cudatoolkit==10.2
pytorch==1.6.0
torchvision==0.7.0

# 2nd place requirements  
tensorflow-gpu>=2.1.0
```

## Training Setup and Hyperparameters

### 1st Place Configuration:
- **Epochs**: 100 (fixed schedule, no early stopping)
- **Learning rate**: 0.001 with ExponentialLR (gamma=0.97)
- **Batch size**: 64-128 (depends on sequence length)
- **Weight decay**: 0.0001
- **Optimizer**: Adam
- **Data augmentation**: 5% base dropout, reverse complement, random subsequences

### 2nd Place Configuration:
- **Epochs**: 21
- **Learning rate**: 0.001 with step decay (0.5x after epoch 15)
- **Batch size**: 24-32
- **Optimizer**: Adam
- **Loss**: Binary crossentropy with label smoothing (0.7-0.9)
- **Sequence length**: 10,000-16,000 bp
- **Embedding**: 8 vocab, 24 dimensions
- **Conv kernels**: [15, 17], [19, 21, 23] depending on model

### Simple Baseline Hyperparameters:
- **K-mer size**: 6-8 (computational efficiency)
- **Learning rate**: 0.01 (higher for faster convergence)
- **Batch size**: 128-256
- **Epochs**: 10-20
- **Simple architecture**: 2 conv layers + 2 dense layers

## Simple Naive Baselines Recommendations

Based on the winning solutions, here are practical baselines for 2-3 datapoints:

### Option 1: K-mer Frequency Baseline
```python
# Extract 6-8 mers, count frequencies
# Use Naive Bayes or Logistic Regression
# ~10-50 lines of code, runs in <1 minute
```

### Option 2: Simple CNN Baseline
```python
# Basic CNN: Embedding -> Conv1D -> GlobalMaxPool -> Dense
# Input: One-hot encoded sequences (truncated to 5000 bp)
# ~100 lines of code, trains in 5-10 minutes
```

### Option 3: BLAST Similarity Baseline
```python
# Use BLAST to find similar sequences
# Simple k-NN on BLAST scores
# Requires BLAST installation but very effective
```

### Option 4: Hybrid Minimal Approach
```python
# Combine 6-mer counts + basic metadata features
# Simple ensemble of Naive Bayes + Logistic Regression  
# Most robust for small datasets
```

## Key Insights for Quick Implementation

1. **Sequence Length**: Even 5,000 bp can work well for baselines
2. **K-mer Size**: 6-8 mers provide good balance of specificity and computational efficiency
3. **Metadata**: Binary features (resistance, growth conditions) are highly informative
4. **Ensemble**: Even simple averaging of 2-3 models improves performance significantly
5. **Data Augmentation**: Reverse complement doubling is essential for DNA sequences

## Evaluation Metric
- **Primary**: Top-10 accuracy (whether correct lab appears in top 10 predictions)
- **Secondary**: Standard accuracy, per-class recall for rare labs
- **Competition specific**: Normalized submissions to match training distribution

This analysis provides a foundation for implementing both sophisticated solutions and practical baselines for genetic engineering attribution tasks.