"""
Virus sequences and restriction sites for genetic engineering
"""

import os
import glob

# Load all NCBI virus sequences
def load_ncbi_sequences():
    """Load all virus sequences from NCBI_sequences directory"""
    def clean_sequence(content):
        """Extract sequence from FASTA content"""
        lines = content.strip().split('\n')
        # Skip header line (starts with >)
        sequence_lines = [line for line in lines if not line.startswith('>')]
        return ''.join(sequence_lines).replace('\n', '').replace(' ', '').upper()
    
    def extract_virus_info(filename):
        """Extract virus information from filename"""
        base_name = os.path.basename(filename).replace('.fasta', '')
        
        # Parse filename to extract virus name and accession
        parts = base_name.split('__')
        if len(parts) >= 2:
            virus_name = parts[0].replace('_', ' ')
            accession = parts[1]
        else:
            # Try to extract accession from the end
            if '_NC_' in base_name:
                virus_name = base_name.split('_NC_')[0].replace('_', ' ')
                accession = 'NC_' + base_name.split('_NC_')[1]
            elif '_MT_' in base_name:
                virus_name = base_name.split('_MT_')[0].replace('_', ' ')
                accession = 'MT_' + base_name.split('_MT_')[1]
            else:
                virus_name = base_name.replace('_', ' ')
                accession = "Unknown"
        
        # Determine virus family based on name
        virus_lower = virus_name.lower()
        if 'coronavirus' in virus_lower or 'sars' in virus_lower or 'mers' in virus_lower:
            family = "coronavirus"
        elif 'adenovirus' in virus_lower:
            family = "adenovirus"
        elif 'influenza' in virus_lower:
            family = "influenza"
        elif 'herpes' in virus_lower or 'alphaherpesvirus' in virus_lower:
            family = "herpesvirus"
        elif 'ebola' in virus_lower:
            family = "filovirus"
        elif 'chikungunya' in virus_lower:
            family = "alphavirus"
        elif 'rhinovirus' in virus_lower:
            family = "picornavirus"
        elif 'measles' in virus_lower or 'mumps' in virus_lower:
            family = "paramyxovirus"
        elif 'rubella' in virus_lower:
            family = "rubivirus"
        elif 'monkeypox' in virus_lower:
            family = "poxvirus"
        elif 'astrovirus' in virus_lower:
            family = "astrovirus"
        elif 'metapneumovirus' in virus_lower or 'parainfluenza' in virus_lower or 'respirovirus' in virus_lower:
            family = "paramyxovirus"
        elif 'henipavirus' in virus_lower:
            family = "henipavirus"
        else:
            family = "unknown"
        
        return virus_name, accession, family, base_name
    
    sequences = {}
    ncbi_dir = "NCBI_sequences"
    
    if not os.path.exists(ncbi_dir):
        print(f"NCBI_sequences directory not found. Please create it and add FASTA files.")
        return {}
    
    fasta_files = glob.glob(os.path.join(ncbi_dir, "*.fasta"))
    
    if not fasta_files:
        print(f"No FASTA files found in {ncbi_dir}")
        return {}
    
    print(f"Loading {len(fasta_files)} virus sequences from NCBI_sequences/")
    
    for fasta_file in fasta_files:
        try:
            with open(fasta_file, "r") as f:
                content = f.read()
            
            sequence = clean_sequence(content)
            virus_name, accession, family, base_name = extract_virus_info(fasta_file)
            
            # Create a clean key for the virus
            key = base_name.replace(' ', '_').replace('-', '_').lower()
            key = ''.join(c for c in key if c.isalnum() or c == '_')
            
            sequences[key] = {
                "name": f"{virus_name} ({accession})",
                "sequence": sequence,
                "length": len(sequence),
                "family": family,
                "accession": accession,
                "filename": os.path.basename(fasta_file)
            }
            
            print(f"  ✅ {virus_name} ({len(sequence):,} bp)")
            
        except Exception as e:
            print(f"  ❌ Error loading {fasta_file}: {e}")
    
    print(f"Successfully loaded {len(sequences)} virus sequences")
    return sequences

# Load all NCBI sequences
VIRUSES = load_ncbi_sequences()

# Extended restriction sites for realistic engineering
RESTRICTION_SITES = {
    # Common 6-cutters
    "BamHI": "GGATCC",
    "EcoRI": "GAATTC", 
    "HindIII": "AAGCTT",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "NotI": "GCGGCCGC",
    "SmaI": "CCCGGG",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "PstI": "CTGCAG",
    "XhoI": "CTCGAG",
    "NdeI": "CATATG",
    "BglII": "AGATCT",
    "NcoI": "CCATGG",
    "SpeI": "ACTAGT",
    "AgeI": "ACCGGT",
    "ApaI": "GGGCCC",
    "MfeI": "CAATTG",
    "EcoRV": "GATATC",
    "ScaI": "AGTACT",
    # 4-cutters for more targets
    "AluI": "AGCT",
    "HaeIII": "GGCC",
    "TaqI": "TCGA",
    "MspI": "CCGG",
    "HpaII": "CCGG",
    "RsaI": "GTAC",
    "DpnI": "GATC",
    "MboI": "GATC",
    "Sau3AI": "GATC"
}

# Safe regions for different virus families
SAFE_REGIONS = {
    "coronavirus": {
        "3_UTR": (24500, 27540),
        "intergenic": (20000, 21000)
    },
    "adenovirus": {
        "ITR": (1, 200),  # Inverted terminal repeats
        "3_UTR": (35000, 35937)
    },
    "influenza": {
        "3_UTR": (1300, 1500),  # Segment-specific
        "intergenic": (100, 300)
    },
    "herpesvirus": {
        "repeat_regions": (1, 1000),
        "3_UTR": (150000, 155000)
    },
    "filovirus": {
        "3_UTR": (18000, 19000),
        "intergenic": (1000, 2000)
    },
    "alphavirus": {
        "3_UTR": (11000, 12000),
        "intergenic": (500, 1000)
    },
    "picornavirus": {
        "3_UTR": (7200, 7400),
        "intergenic": (100, 300)
    },
    "paramyxovirus": {
        "3_UTR": (15000, 16000),
        "intergenic": (1000, 2000)
    },
    "rubivirus": {
        "3_UTR": (9700, 9800),
        "intergenic": (100, 200)
    },
    "poxvirus": {
        "inverted_repeats": (1, 500),
        "3_UTR": (190000, 200000)
    },
    "astrovirus": {
        "3_UTR": (6800, 6900),
        "intergenic": (100, 200)
    },
    "henipavirus": {
        "3_UTR": (18000, 19000),
        "intergenic": (1000, 2000)
    }
}