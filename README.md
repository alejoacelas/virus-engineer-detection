# Comprehensive Virus Engineering System

A more realistic virus engineering pipeline using 28 NCBI whole genome sequences that generates reproducible datasets with 10 different genetic engineering methods targeting restriction sites.

## 🧬 Features

- **28 NCBI Viruses**: Coronavirus, Adenovirus, Herpesvirus, Influenza, and more
- **10 Engineering Methods**: Region inversion, duplication, deletion, GFP/Tn5 insertion, frameshift mutations, SNPs, CRISPR simulation
- **Restriction Site Targeting**: Uses 29 different restriction enzymes for realistic engineering
- **Reproducible Results**: Comprehensive seed control for deterministic dataset generation
- **Full Genomes**: Complete viral sequences (15K-200K bp) with complete metadata
- **Balanced Datasets**: Configurable examples per virus with balanced method distribution

## 📁 Files

- `enhanced_virus_engineer_restriction.py` - Main engineering system with seed control
- `virus_sequences.py` - NCBI sequence loader with 29 restriction enzymes
- `NCBI_sequences/` - 28 real viral genome sequences in FASTA format
- `.gitignore` - Excludes large generated datasets

## 🚀 Usage

### Basic Usage
```bash
# Generate dataset with default seed (42) - reproducible results
python enhanced_virus_engineer_restriction.py

# Generate with custom seed
python enhanced_virus_engineer_restriction.py --seed 123

# Generate different dataset each time
python enhanced_virus_engineer_restriction.py --seed random

# Generate smaller dataset (100 per virus instead of 1000)
python enhanced_virus_engineer_restriction.py --num-per-virus 100

# Use environment variable for seed
VIRUS_SEED=999 python enhanced_virus_engineer_restriction.py

# No seed - truly random results each time
python enhanced_virus_engineer_restriction.py --no-seed
```

### Command Line Options
- `--seed` - Random seed for reproducibility (default: 42, or "random" for different results)
- `--num-per-virus` - Number of examples per virus (default: 1000)
- `--no-seed` - Don't set any seed (truly random results each time)

## 🦠 Available Viruses (28 total)

### Coronavirus Family
- SARS-CoV-2 (NC_045512.2) - 29,903 bp
- MERS virus (NC_019843.3) - 30,119 bp
- Human Coronavirus NL63 (NC_005831.2) - 27,553 bp
- Human Coronavirus OC43 (NC_006213.1) - 30,738 bp
- Human Coronavirus 229E (NC_002645.1) - 27,317 bp

### Adenovirus Family
- Human Adenovirus 1 (NC_001405.1) - 35,937 bp
- Human Adenovirus 5 (NC_001405.1) - 35,937 bp
- Human Adenovirus 54 (NC_012959.1) - 34,920 bp

### Herpesvirus Family
- Human Herpesvirus 4 Type 2 (NC_009334.1) - 171,823 bp
- Human Alphaherpesvirus 3 (NC_001348.1) - 124,884 bp

### Influenza Family
- Influenza A H1N1 (NC_026433.1) - 13,588 bp
- Influenza A H2N2 (NC_007357.1) - 13,588 bp
- Influenza A H3N2 (NC_007361.1) - 13,588 bp

### Other Families
- Ebola virus (NC_002549.1) - 18,959 bp
- Monkeypox virus (NC_063383.1) - 196,858 bp
- Chikungunya virus (NC_004162.2) - 11,805 bp
- Rubella virus (NC_001545.1) - 9,762 bp
- And 12 more...

## 🔬 Engineering Methods

### 1. Region Inversion
- Inverts DNA regions at restriction sites
- Realistic size variations (100-2000 bp)
- Biological randomness with flanking sequences

### 2. Region Duplication
- Duplicates genomic regions at restriction sites
- Maintains biological context
- Variable duplication sizes

### 3. Region Deletion
- Deletes regions using restriction enzymes
- Safe deletion sizes (50-1000 bp)
- Preserves genome integrity

### 4. GFP Insertion
- Inserts Green Fluorescent Protein reporter (720 bp)
- Uses restriction sites for precise integration
- Includes proper flanking sequences

### 5. Tn5 Insertion
- Inserts Tn5 transposon sequence (1,200 bp)
- Commonly used in synthetic biology
- Realistic integration patterns

### 6. iGEM GFP Insertion
- Inserts iGEM-specific GFP variant (750 bp)
- Popular in synthetic biology competitions
- Optimized for bacterial expression

### 7. Frameshift Mutations
- Creates insertion/deletion mutations at restriction sites
- Causes reading frame shifts
- Realistic mutation sizes (1-10 bp)

### 8. Single Nucleotide Changes
- Point mutations at restriction sites
- Transition and transversion mutations
- Biologically relevant changes

### 9. Random Substitutions
- Replaces regions with random sequences
- Maintains GC content patterns
- Variable substitution lengths

### 10. CRISPR Simulation
- Simulates CRISPR-Cas9 editing
- Includes PAM sequences (NGG)
- Realistic gRNA targeting

## 🧪 Restriction Enzymes (29 total)

### 6-Cutters (Common)
- **BamHI**: GGATCC
- **EcoRI**: GAATTC
- **HindIII**: AAGCTT
- **XbaI**: TCTAGA
- **SalI**: GTCGAC
- **NotI**: GCGGCCGC
- **SmaI**: CCCGGG
- **KpnI**: GGTACC
- **SacI**: GAGCTC
- **PstI**: CTGCAG
- **XhoI**: CTCGAG
- **NdeI**: CATATG
- **BglII**: AGATCT
- **NcoI**: CCATGG
- **SpeI**: ACTAGT
- **AgeI**: ACCGGT
- **ApaI**: GGGCCC
- **MfeI**: CAATTG
- **EcoRV**: GATATC
- **ScaI**: AGTACT

### 4-Cutters (Frequent)
- **AluI**: AGCT
- **DdeI**: CTNAG
- **HaeIII**: GGCC
- **HinfI**: GANTC
- **MboI**: GATC
- **RsaI**: GTAC
- **TaqI**: TCGA
- **HpaII**: CCGG
- **MspI**: CCGG

## 📊 Output

### Individual FASTA Files
Each engineered genome is saved as a separate FASTA file with:
- Complete modified sequence
- NCBI accession numbers in headers
- Engineering method and restriction enzyme used
- Position and size information

### Metadata JSON
Complete results stored in `restriction_virus_engineered.json`:
- All engineering parameters
- Original and modified sequences
- Restriction site information
- Method-specific details
- Reproducibility information

## 🎲 Reproducibility

The system includes comprehensive seed control:

- **Default seed (42)**: Ensures reproducible results by default
- **Custom seeds**: Use any integer for different but deterministic datasets
- **Random seeds**: Use `--seed random` for different results each time
- **Environment variables**: Set `VIRUS_SEED` for automated workflows
- **No seed option**: Use `--no-seed` for truly random results

### Example Reproducibility
```bash
# These will generate identical datasets
python enhanced_virus_engineer_restriction.py --seed 42
python enhanced_virus_engineer_restriction.py --seed 42

# This will generate a different but reproducible dataset
python enhanced_virus_engineer_restriction.py --seed 123
```

## 🔬 Biological Realism

- **Restriction site targeting**: 70% of operations use actual restriction sites
- **Biological randomness**: Size variations, flanking sequences, offsets
- **Real viral sequences**: All 28 viruses from NCBI GenBank
- **Proper FASTA format**: Standard headers with accession numbers
- **Complete genomes**: Full-length sequences, not fragments

## 📈 Dataset Generation

### Default Configuration
- **28 viruses** × **1000 examples** = **28,000 engineered genomes**
- **10 methods** = **100 examples per method per virus**
- **Balanced distribution** across all engineering methods
- **~1GB total** in FASTA files + metadata

### Customizable Parameters
- `--num-per-virus`: Control total dataset size
- Seed control: Ensure reproducibility
- Method balancing: Automatic distribution across all methods

## 🚫 Excluded Files

The following files are excluded from the repository (via `.gitignore`):
- `restriction_genomes/` - Directory with 28,000+ FASTA files
- `restriction_virus_engineered.json` - Large metadata file (~1GB)

These files are generated locally when running the pipeline and should not be committed to version control due to their size.

## 🔧 Requirements

- Python 3.6+
- Standard library only (no external dependencies)
- ~2GB free disk space for full dataset generation

## 📝 Example Output

```
>SARS-CoV-2 (NC_045512.2) | insert_gfp | 29,903 -> 30,623 bp | BamHI at 1234
[Full 30,623 bp engineered genome sequence...]

>Human Adenovirus 5 (NC_001405.1) | region_deletion | 35,937 -> 35,437 bp | EcoRI at 5678
[Full 35,437 bp engineered genome sequence...]
```

## 🎯 Use Cases

- **Machine Learning**: Training classifiers for engineered vs. natural sequences
- **Bioinformatics Research**: Studying genetic engineering signatures
- **Synthetic Biology**: Benchmarking engineering methods
- **Biosecurity**: Developing detection algorithms for engineered sequences
- **Education**: Teaching genetic engineering concepts

---

*This system generates realistic viral genome engineering datasets for research and educational purposes. All sequences are based on publicly available NCBI data.*