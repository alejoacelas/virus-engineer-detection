#!/usr/bin/env python3
"""
Quick Test Script for Baseline Models
Tests that accuracy improves during training on dummy data
"""

import sys
sys.path.append('utils')
sys.path.append('baselines')

from utils.generate_engineering_data import save_engineering_detection_data
from utils.metrics import compare_with_calibration

def test_data_creation():
    """Test engineering detection data generation"""
    print("=" * 50)
    print("TESTING ENGINEERING DETECTION DATA GENERATION")
    print("=" * 50)
    try:
        # Create engineering detection dataset for testing
        detection_data = save_engineering_detection_data(
            csv_path="data/viral_genome_df_prelim.csv",
            output_path="data/engineering_detection_train.csv",
            n_samples=20000,  # Large dataset for robust testing
            # skip_if_exists=True,
            mean_length=187,
            sd_length=50,
            engineering_rate=0.5
        )

        print(f"✓ Successfully created {len(detection_data)} samples")
        print(f"✓ Classes: {sorted(detection_data['label'].unique())}")
        print(f"✓ Sequence lengths: {detection_data['sequence'].str.len().describe()}")
        print(f"✓ Engineering rate: {detection_data['label'].mean():.3f}")

        return True

    except Exception as e:
        print(f"✗ Data creation failed: {e}")
        return False

def test_kmer_baseline():
    """Test k-mer Logistic Regression baseline"""
    print("\n" + "=" * 50)
    print("TESTING K-MER LOGISTIC REGRESSION BASELINE")
    print("=" * 50)

    try:
        from baselines.baseline_kmer_logistic import train_kmer_baseline, train_multi_kmer_baseline

        results = train_kmer_baseline(
            csv_path="data/engineering_detection_train.csv",
            k=6,  # Larger k-mers for better pattern detection
            max_kmers=1000,  # More features for complex patterns
            test_size=0.2
        )
        # results = train_multi_kmer_baseline(
        #     csv_path="data/engineering_detection_train.csv",
        #     k_values=[4, 6, 8],  # Larger k-mers for better pattern detection
        #     max_kmers=1000,  # More features for complex patterns
        #     test_size=0.2
        # )

        # Compare to naive fraction baseline with calibration
        comparison = compare_with_calibration(
            results['y_train_true'], results['y_train_proba'][:, 1],
            results['y_test_true'], results['y_test_proba'][:, 1],
            test_pos_label_fraction=0.02
        )
        
        model_precision = comparison['model']['precision']
        model_recall = comparison['model']['recall']
        naive_precision = comparison['naive']['precision']
        naive_recall = comparison['naive']['recall']
        pos_fraction = comparison['positive_fraction']
        
        print(f"\nK-mer Logistic Regression Results:")
        print(f"Model    - Precision: {model_precision:.3f}, Recall: {model_recall:.3f}")
        print(f"Naive    - Precision: {naive_precision:.3f}, Recall: {naive_recall:.3f}")
        print(f"Positive fraction: {pos_fraction:.3f}")
        
        # Check if model outperforms naive baseline in both metrics
        precision_better = model_precision > naive_precision
        recall_better = model_recall > naive_recall
        
        if precision_better and recall_better:
            print("✓ K-mer baseline outperforms naive baseline in both precision and recall")
            return True
        else:
            print("✗ K-mer baseline does not outperform naive baseline in both metrics")
            return False

    except Exception as e:
        print(f"✗ K-mer baseline failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_cnn_baseline():
    """Test CNN baseline"""
    print("\n" + "=" * 50)
    print("TESTING CNN BASELINE")
    print("=" * 50)

    try:
        from baselines.baseline_simple_cnn import train_simple_cnn, train_improved_cnn

        # results = train_simple_cnn(
        #     csv_path="data/engineering_detection_train.csv",
        #     max_len=500,   # Match typical segment length
        #     epochs=5,    # More epochs for learning complex patterns
        #     batch_size=128  # Larger batches for stability
        # )
        results = train_improved_cnn(
            csv_path="data/engineering_detection_train.csv",
            max_len=500,   # Match typical segment length
            epochs=5,    # More epochs for learning complex patterns
            batch_size=128  # Larger batches for stability
        )

        # Compare to naive fraction baseline with calibration
        comparison = compare_with_calibration(
            results['y_train_true'], results['y_train_proba'][:, 1],
            results['y_test_true'], results['y_test_proba'][:, 1],
            test_pos_label_fraction=0.02
        )
        
        model_precision = comparison['model']['precision']
        model_recall = comparison['model']['recall']
        naive_precision = comparison['naive']['precision']
        naive_recall = comparison['naive']['recall']
        pos_fraction = comparison['positive_fraction']
        
        print(f"\nSimple CNN Results:")
        print(f"Model    - Precision: {model_precision:.3f}, Recall: {model_recall:.3f}")
        print(f"Naive    - Precision: {naive_precision:.3f}, Recall: {naive_recall:.3f}")
        print(f"Positive fraction: {pos_fraction:.3f}")
        
        # Check if model outperforms naive baseline in both metrics
        precision_better = model_precision > naive_precision
        recall_better = model_recall > naive_recall
        
        if precision_better and recall_better:
            print("✓ CNN baseline outperforms naive baseline in both precision and recall")
            return True
        else:
            print("✗ CNN baseline does not outperform naive baseline in both metrics")
            return False

    except Exception as e:
        print(f"✗ CNN baseline failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("Testing Genetic Engineering Detection Baselines")
    print("This will create engineering detection data and test that models can distinguish engineered sequences\n")

    results = []

    # Test data creation
    results.append(("Data Creation", test_data_creation()))

    # Test baselines if data creation succeeded
    if results[0][1]:
        results.append(("K-mer Baseline", test_kmer_baseline()))
        results.append(("CNN Baseline", test_cnn_baseline()))

    # Summary
    print("\n" + "=" * 50)
    print("TEST SUMMARY")
    print("=" * 50)

    passed = 0
    for test_name, test_success in results:
        status = "✓ PASS" if test_success else "✗ FAIL"
        print(f"{test_name:<20} {status}")
        if test_success:
            passed += 1

    print(f"\nOverall: {passed}/{len(results)} tests passed")

    if passed == len(results):
        print("\n🎉 All tests passed! Baselines are working correctly.")
    else:
        print("\n⚠️  Some tests failed. Check the error messages above.")

    return passed == len(results)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
