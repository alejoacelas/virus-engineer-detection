#!/usr/bin/env python3
"""
Complete virus sequences module with ORF annotations for all viruses
"""

import os
import glob
import json
import random

def load_all_virus_annotations():
    """Load all virus ORF annotations from JSON file"""
    annotation_file = "all_virus_annotations.json"
    
    if not os.path.exists(annotation_file):
        print(f"Warning: {annotation_file} not found. Run parse_all_annotations.py first.")
        return {}
    
    with open(annotation_file, 'r') as f:
        data = json.load(f)
    
    print(f"Loaded annotations for {len(data)} viruses")
    return data

def load_ncbi_sequences_with_all_orfs():
    """Load all virus sequences from NCBI_sequences directory with complete ORF support"""
    
    def clean_sequence(content):
        """Extract sequence from FASTA content"""
        lines = content.strip().split('\n')
        sequence_lines = [line for line in lines if not line.startswith('>')]
        return ''.join(sequence_lines).replace('\n', '').replace(' ', '').upper()
    
    def extract_virus_info(filename):
        """Extract virus information from filename"""
        base_name = os.path.basename(filename).replace('.fasta', '')
        
        parts = base_name.split('__')
        if len(parts) >= 2:
            virus_name = parts[0].replace('_', ' ')
            accession = parts[1]
        else:
            if '_NC_' in base_name:
                virus_name = base_name.split('_NC_')[0].replace('_', ' ')
                accession = 'NC_' + base_name.split('_NC_')[1]
            elif '_MT_' in base_name:
                virus_name = base_name.split('_MT_')[0].replace('_', ' ')
                accession = 'MT_' + base_name.split('_MT_')[1]
            elif ' NC_' in base_name:
                virus_name = base_name.split(' NC_')[0].replace('_', ' ')
                accession = 'NC_' + base_name.split(' NC_')[1]
            elif 'MT' in base_name and '.' in base_name:
                parts = base_name.split('_')
                for part in parts:
                    if part.startswith('MT') and '.' in part:
                        accession = part
                        virus_name = base_name.replace('_' + part, '').replace('_', ' ')
                        break
                else:
                    virus_name = base_name.replace('_', ' ')
                    accession = "Unknown"
            else:
                virus_name = base_name.replace('_', ' ')
                accession = "Unknown"
        
        # Determine family
        family = "Unknown"
        if "coronavirus" in virus_name.lower() or "sars-cov" in virus_name.lower() or "mers" in virus_name.lower():
            family = "Coronavirus"
        elif "influenza" in virus_name.lower() or "parainfluenza" in virus_name.lower():
            family = "Influenza"
        elif "adenovirus" in virus_name.lower():
            family = "Adenovirus"
        elif "herpesvirus" in virus_name.lower():
            family = "Herpesvirus"
        elif "rhinovirus" in virus_name.lower() or "picornavirus" in virus_name.lower():
            family = "Picornavirus"
        elif "ebola" in virus_name.lower() or "filovirus" in virus_name.lower():
            family = "Filovirus"
        elif "chikungunya" in virus_name.lower() or "alphavirus" in virus_name.lower():
            family = "Alphavirus"
        elif "rubella" in virus_name.lower() or "rubivirus" in virus_name.lower():
            family = "Rubivirus"
        elif "monkeypox" in virus_name.lower() or "poxvirus" in virus_name.lower():
            family = "Poxvirus"
        elif "astrovirus" in virus_name.lower():
            family = "Astrovirus"
        elif "henipavirus" in virus_name.lower() or "nipah" in virus_name.lower():
            family = "Henipavirus"
        elif "respirovirus" in virus_name.lower() or "mumps" in virus_name.lower() or "measles" in virus_name.lower() or "metapneumovirus" in virus_name.lower():
            family = "Paramyxovirus"

        return virus_name, accession, family

    # Load all virus annotations
    all_annotations = load_all_virus_annotations()
    
    viruses = {}
    ncbi_dir = "NCBI_sequences"
    if not os.path.exists(ncbi_dir):
        print(f"Error: NCBI_sequences directory not found at {ncbi_dir}")
        return {}

    print(f"Loading {len(glob.glob(os.path.join(ncbi_dir, '*.fasta')))} virus sequences from {ncbi_dir}/")
    
    for filepath in glob.glob(os.path.join(ncbi_dir, '*.fasta')):
        with open(filepath, "r") as f:
            content = f.read()
            seq = clean_sequence(content)
            virus_name, accession, family = extract_virus_info(filepath)
            
            key = os.path.basename(filepath).replace('.fasta', '').lower().replace('.', '').replace('-', '_')
            
            # Get ORF annotations for this virus
            orfs = None
            if accession in all_annotations:
                orfs = all_annotations[accession]['orfs']
                print(f"  {virus_name} ({len(seq):,} bp) - WITH ORF annotations ({len(orfs)} CDS)")
            else:
                print(f"  {virus_name} ({len(seq):,} bp) - No ORF annotations available")
            
            viruses[key] = {
                "name": f"{virus_name} ({accession})",
                "sequence": seq,
                "length": len(seq),
                "family": family,
                "accession": accession,
                "orfs": orfs
            }
    
    print(f"Successfully loaded {len(viruses)} virus sequences")
    return viruses

# ORF-aware helper functions
def get_orfs_for_virus(virus_key):
    """Get ORF annotations for a specific virus"""
    if virus_key in VIRUSES and VIRUSES[virus_key]['orfs']:
        return VIRUSES[virus_key]['orfs']
    return None

def find_orfs_in_region(start_pos, end_pos, orfs):
    """Find ORFs that overlap with a given region"""
    if not orfs:
        return []
    
    overlapping_orfs = []
    for orf in orfs:
        orf_start = orf['start']
        orf_end = orf['end']
        
        # Check for overlap
        if not (end_pos < orf_start or start_pos > orf_end):
            overlapping_orfs.append(orf)
    
    return overlapping_orfs

def get_orf_boundaries(orfs):
    """Get start and stop positions of all ORFs for safe insertion sites"""
    if not orfs:
        return []
    
    boundaries = []
    for orf in orfs:
        boundaries.extend([orf['start'], orf['end']])
    
    return sorted(boundaries)

def get_safe_insertion_sites(sequence_length, orfs, num_sites=10):
    """Get safe insertion sites (between ORFs or in non-coding regions)"""
    if not orfs:
        # No ORF data, return random sites
        return random.sample(range(100, sequence_length - 100), min(num_sites, sequence_length - 200))
    
    # Get ORF boundaries
    boundaries = get_orf_boundaries(orfs)
    
    # Find gaps between ORFs
    safe_sites = []
    
    # Check region before first ORF
    if boundaries and boundaries[0] > 100:
        safe_sites.extend(range(100, min(boundaries[0] - 50, 1000)))
    
    # Check gaps between ORFs
    for i in range(len(boundaries) - 1):
        gap_start = boundaries[i] + 50
        gap_end = boundaries[i + 1] - 50
        
        if gap_end - gap_start > 100:  # Only consider gaps > 100bp
            gap_size = min(500, gap_end - gap_start)  # Limit gap size
            sites_in_gap = list(range(gap_start, gap_start + gap_size))
            safe_sites.extend(random.sample(sites_in_gap, min(3, len(sites_in_gap))))
    
    # Check region after last ORF
    if boundaries:
        last_boundary = boundaries[-1]
        if sequence_length - last_boundary > 100:
            end_sites = list(range(last_boundary + 50, sequence_length - 100))
            safe_sites.extend(random.sample(end_sites, min(3, len(end_sites))))
    
    # If we don't have enough safe sites, add some random ones
    while len(safe_sites) < num_sites:
        random_site = random.randint(100, sequence_length - 100)
        if random_site not in safe_sites:
            safe_sites.append(random_site)
    
    return safe_sites[:num_sites]

def get_orf_for_position(position, orfs):
    """Find which ORF (if any) contains the given position"""
    if not orfs:
        return None
    
    for orf in orfs:
        if orf['start'] <= position <= orf['end']:
            return orf
    
    return None

def get_orf_insertion_strategies(orfs, sequence_length):
    """Get realistic insertion strategies based on ORF structure"""
    strategies = []
    
    if not orfs:
        return ["random_insertion"]
    
    # Strategy 1: Insert after stop codons
    for orf in orfs:
        if orf['end'] + 100 < sequence_length:
            strategies.append({
                "type": "after_stop_codon",
                "target_orf": orf['name'],
                "position": orf['end'] + random.randint(10, 50),
                "description": f"Insert after {orf['name']} stop codon"
            })
    
    # Strategy 2: Insert between genes
    boundaries = get_orf_boundaries(orfs)
    for i in range(len(boundaries) - 1):
        gap_start = boundaries[i] + 50
        gap_end = boundaries[i + 1] - 50
        if gap_end - gap_start > 200:
            strategies.append({
                "type": "intergenic",
                "position": random.randint(gap_start, gap_end),
                "description": "Insert between genes"
            })
    
    # Strategy 3: Insert in non-essential ORFs (smaller ones)
    small_orfs = [orf for orf in orfs if orf['end'] - orf['start'] < 1000]
    for orf in small_orfs[:3]:  # Limit to first 3 small ORFs
        insertion_pos = orf['start'] + random.randint(50, 200)
        if insertion_pos < orf['end'] - 50:
            strategies.append({
                "type": "within_non_essential",
                "target_orf": orf['name'],
                "position": insertion_pos,
                "description": f"Insert within {orf['name']} (non-essential)"
            })
    
    return strategies

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
    # 4-cutters (more frequent)
    "AluI": "AGCT",
    "DdeI": "CTNAG",
    "HaeIII": "GGCC",
    "HinfI": "GANTC",
    "MboI": "GATC",
    "RsaI": "GTAC",
    "TaqI": "TCGA",
    "HpaII": "CCGG",
    "MspI": "CCGG"
}

# Load the complete virus data with all ORFs
VIRUSES = load_ncbi_sequences_with_all_orfs()

# Print comprehensive summary
print("\nCOMPLETE VIRUS ENGINEERING SYSTEM WITH ORF SUPPORT")
print("=" * 60)
print(f"Total viruses: {len(VIRUSES)}")
print(f"Viruses with ORF annotations: {sum(1 for v in VIRUSES.values() if v['orfs'] is not None)}")
print(f"Restriction enzymes: {len(RESTRICTION_SITES)}")

# Show ORF summary by family
families_with_orfs = {}
for key, virus in VIRUSES.items():
    if virus['orfs']:
        family = virus['family']
        if family not in families_with_orfs:
            families_with_orfs[family] = []
        families_with_orfs[family].append({
            'name': virus['name'],
            'orfs': len(virus['orfs'])
        })

if families_with_orfs:
    print("\nORF ANNOTATIONS BY FAMILY:")
    for family, viruses in sorted(families_with_orfs.items()):
        print(f"  {family}: {len(viruses)} viruses with ORFs")
        for virus in viruses[:2]:  # Show first 2
            print(f"    - {virus['name']}: {virus['orfs']} ORFs")
        if len(viruses) > 2:
            print(f"    ... and {len(viruses) - 2} more")

print("\nSystem ready for comprehensive ORF-aware virus engineering!")
