# Real Virus Engineering System

A realistic virus engineering pipeline using actual NCBI sequences that generates full genomes with modifications at restriction sites.

## Features

✅ **Real NCBI Viruses**: Coronavirus NL63 and Adenovirus 54  
✅ **Full Genomes**: Complete viral sequences (28K-35K bp)  
✅ **Realistic Engineering**: Uses actual restriction sites when available  
✅ **Safe Fallbacks**: Uses safe genomic regions when no restriction sites found  
✅ **Multiple Operations**: GFP insertion, deletion, substitution  
✅ **Individual Output**: Each genome saved as separate FASTA file  

## Files

- `virus_engineer.py` - Main engineering system
- `virus_sequences.py` - Real virus database loader
- `corona_real.txt` - Human Coronavirus NL63 (NC_005831.2) - 28,393 bp
- `adeno_real.txt` - Human Adenovirus 54 (NC_012959.1) - 34,920 bp
- `engineered_genomes/` - Directory with individual FASTA files
- `multi_virus_engineered.json` - Complete results with metadata

## Usage

```bash
# Run engineering simulation with real viral sequences
python virus_engineer.py
```

## Real Viral Sequences

### Coronavirus NL63 (NC_005831.2)
- **Length**: 28,393 bp
- **Family**: Coronavirus
- **Source**: NCBI GenBank

### Adenovirus 54 (NC_012959.1)  
- **Length**: 34,920 bp
- **Family**: Adenovirus
- **Source**: NCBI GenBank

## Engineering Operations

### 1. GFP Insertion
- Inserts GFP reporter gene (720 bp)
- Uses restriction sites when available
- Falls back to safe regions

### 2. Region Deletion  
- Deletes between restriction sites
- Realistic deletion sizes (50-500 bp)
- Safe region fallback

### 3. Region Substitution
- Replaces regions with synthetic DNA
- Maintains length or allows size changes
- Uses restriction sites for precision

## Restriction Enzymes

- **BamHI**: GGATCC
- **EcoRI**: GAATTC  
- **HindIII**: AAGCTT
- **XbaI**: TCTAGA
- **SalI**: GTCGAC
- **NotI**: GCGGCCGC
- **SmaI**: CCCGGG

## Output

Each engineered genome includes:
- **Full sequence** in FASTA format
- **NCBI accession numbers** in headers
- **Modification details** and enzyme information
- **Length changes** from original
- **Complete metadata** in JSON

## Example Output

```
>Human Coronavirus NL63 (NC_005831.2) | insert_gfp | 28,393 -> 29,113 bp | XbaI
[Full 29,113 bp engineered genome sequence...]
```

## Key Features

- **Realistic**: Uses actual NCBI viral sequences
- **Full genomes**: Complete viral sequences, not fragments  
- **Biologically accurate**: Proper genomic regions and restriction sites
- **NCBI compliant**: Includes accession numbers and proper headers
- **Clean**: Simple, focused codebase with real data