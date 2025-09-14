# Virus engineering simulation pipeline

A more realistic virus engineering pipeline using 25 NCBI whole genome sequences that generates reproducible datasets with 4 different genetic engineering methods.

##  Features

- **25 NCBI Viruses**: Coronavirus, Adenovirus, Herpesvirus, and more (complete genomes only)
- **4 Engineering Methods**: Region inversion, deletion, gene deletion, frameshift mutations
- **ORF-Aware Engineering**: Uses GenBank annotations to target specific genes and regions
- **Full Genomes**: Complete viral sequences (6K-200K bp) with modification coordinates
- **Balanced Datasets**: 1000 examples per virus with balanced method distribution

## Files

- `enhanced_sequence_generator.py` - Main sequence generator with ORFs
- `viral_genetic_engineering.py` - Core engineering methods with partial gene deletion
- `virus_sequences_all_orfs.py` - NCBI sequence loader with ORF annotations
- `all_virus_annotations.json` - Parsed GenBank ORF annotations for all viruses
- `annotations/` - GenBank files with gene and ORF information
- `NCBI_sequences/` - 25 real viral genome sequences in FASTA format
- `engineered_sequences/` - Generated engineered sequences (25,000 total)
- `.gitignore` - Excludes large generated datasets

## Usage

### Basic Usage
```bash
# Generate 25,000 engineered sequences (1000 per virus)
python enhanced_sequence_generator.py

# The system will automatically:
# - Load 25 complete viral genomes with ORF annotations
# - Apply 13 different engineering methods
# - Generate 1000 sequences per virus with biological realism
# - Save results to engineered_sequences/
```

### Output Files
- **JSON files**: Detailed metadata with modification coordinates (start/end positions)
- **FASTA files**: Full engineered genome sequences ready for analysis
- **Summary**: Complete generation statistics and method distribution

## Available Viruses (25 total - complete genomes only)

### Coronavirus Family (5 viruses)
- SARS-CoV-2 (NC_045512.2) - 29,903 bp, 12 ORFs
- MERS virus (NC_019843.3) - 30,119 bp, 11 ORFs
- Human Coronavirus NL63 (NC_005831.2) - 27,553 bp, 7 ORFs
- Human Coronavirus OC43 (NC_006213.1) - 30,741 bp, 10 ORFs
- Human Coronavirus 229E (NC_002645.1) - 27,317 bp, 8 ORFs

### Adenovirus Family (3 viruses)
- Human Adenovirus 1 (NC_001405.1) - 35,937 bp, 38 ORFs
- Human Adenovirus 5 (NC_001405.1) - 35,937 bp, 38 ORFs
- Human Adenovirus 54 (NC_012959.1) - 34,920 bp, 36 ORFs

### Herpesvirus Family (2 viruses)
- Human Herpesvirus 4 Type 2 (NC_009334.1) - 172,764 bp, 80 ORFs
- Human Alphaherpesvirus 3 (NC_001348.1) - 124,884 bp, 73 ORFs

### Other Families (15 viruses)
- Ebola virus (NC_002549.1) - 18,959 bp, 9 ORFs
- Monkeypox virus (NC_063383.1) - 197,209 bp, 179 ORFs
- Chikungunya virus (NC_004162.2) - 11,826 bp, 2 ORFs
- Rubella virus (NC_001545.1) - 9,755 bp, 2 ORFs
- Astrovirus MLB1 (NC_011400) - 6,171 bp, 3 ORFs
- And 10 more complete viral genomes...

## Engineering Methods (4 total)

### 1. Region Inversion
- Inverts DNA regions at restriction sites avoiding ORFs
- Uses 29 different restriction enzymes (BamHI, EcoRI, HindIII, etc.)
- Variable size variations (50-500 bp) with biological randomness
- Falls back to random positioning if no safe restriction sites available

### 2. Region Deletion
- Deletes regions at restriction sites between ORFs
- Uses restriction enzyme cutting sites for precise targeting
- Safe deletion sizes (50-500 bp) with position tolerance
- Preserves genome integrity and gene function

### 3. Gene Deletion
- Deletes 10-50% of genes (biologically realistic)
- Random deletion positions (start, middle, end)

### 4. Frameshift Mutations
- Creates 1bp insertion/deletion mutations
- Causes reading frame shifts
- Targets ORF regions


##  Restriction Enzymes (29 total)

### 6-Cutters 
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

### 4-Cutters 
- **AluI**: AGCT
- **DdeI**: CTNAG
- **HaeIII**: GGCC
- **HinfI**: GANTC
- **MboI**: GATC
- **RsaI**: GTAC
- **TaqI**: TCGA
- **HpaII**: CCGG
- **MspI**: CCGG

## Output

### JSON Files (Detailed Metadata)
Each virus has a detailed JSON file with:
- Complete modification coordinates (start/end positions)
- Original gene information and ORF annotations
- Deletion percentages and positions for gene deletions
- Method-specific parameters and randomness factors
- Full sequence metadata

### FASTA Files (Full Genomes)
Each virus has a FASTA file containing:
- Complete modified genome sequences
- Headers with method names and length changes
- 1000 engineered sequences per virus
- Ready for analysis and classification

### Example Output Structure
```
engineered_sequences/
├── astrovirus_mlb1__nc_011400_detailed_sequences.json
├── astrovirus_mlb1__nc_011400_full_genomes.fasta
├── sars_cov_2_nc_0455122_detailed_sequences.json
├── sars_cov_2_nc_0455122_full_genomes.fasta
└── generation_summary.json
```



### Example: Astrovirus Gene Deletion
- **Before**: Full ORF1b deletion (1537 bp) - would make virus non-viable
- **After**: Partial deletion (154-766 bp, avg 30.6%) - creates attenuated virus
- **Result**: More realistic for classifier training with natural-like modifications

## Dataset Generation

### Current Configuration
- **25 viruses** × **1000 examples** = **25,000 engineered genomes**
- **4 methods** with balanced distribution
- **Randomness** ensures 1000 sequences per virus
- **Complete metadata** with modification coordinates



## Requirements

- Python 3.6+
- Standard library only (no external dependencies)
- ~2GB free disk space for full dataset generation

## Example Output

### FASTA Header Format
```
>gene_del_2314_3850 | Gene Deletion | -284bp
[Full engineered genome sequence with partial gene deletion...]

>region_inversion_500_1200 | Region Inversion | +0bp
[Full engineered genome sequence with inverted region...]
```

### JSON Metadata Example
```json
{
  "id": "gene_del_2314_3850",
  "method": "gene_deletion",
  "length_change": -284,
  "modifications": [{
    "start_position": 2314,
    "end_position": 2597,
    "details": {
      "deleted_gene": "ORF1b",
      "deletion_percentage": 18.5,
      "deletion_position": "start"
    }
  }]
}
```



---

*This system generates realistic viral genome engineering datasets for research and educational purposes. All sequences are based on publicly available NCBI data.*
