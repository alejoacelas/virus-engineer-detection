"""
Enhanced Virus Engineering System with Restriction Sites
Uses actual restriction sites for realistic genetic engineering
"""

import random
import json
import os
from virus_sequences import VIRUSES, RESTRICTION_SITES

def get_clean_sequence(sequence):
    """Clean sequence string"""
    return sequence.replace('\n', '').replace(' ', '').upper()

def find_restriction_sites(sequence, enzyme):
    """Find all restriction sites in sequence"""
    site = RESTRICTION_SITES[enzyme]
    positions = []
    start = 0
    while True:
        pos = sequence.find(site, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions

def get_random_restriction_site(sequence, preferred_enzymes=None):
    """Get a random restriction site from the sequence"""
    if preferred_enzymes is None:
        preferred_enzymes = ["BamHI", "EcoRI", "HindIII", "XbaI", "SalI", "NotI", "KpnI", "SacI", "PstI"]
    
    available_sites = {}
    
    # Try preferred enzymes first
    for enzyme in preferred_enzymes:
        if enzyme in RESTRICTION_SITES:
            sites = find_restriction_sites(sequence, enzyme)
            if sites:
                available_sites[enzyme] = sites
    
    # If no preferred sites, try all enzymes
    if not available_sites:
        for enzyme in RESTRICTION_SITES:
            sites = find_restriction_sites(sequence, enzyme)
            if sites:
                available_sites[enzyme] = sites
    
    if not available_sites:
        return None, None
    
    # Choose random enzyme and site
    enzyme = random.choice(list(available_sites.keys()))
    position = random.choice(available_sites[enzyme])
    
    return enzyme, position

def add_biological_randomness():
    """Add biological randomness to sequences"""
    # Add random flanking sequences (common in cloning)
    flanking_options = [
        "",  # No flanking
        "GAATTC",  # EcoRI site
        "GGATCC",  # BamHI site  
        "AAGCTT",  # HindIII site
        "GATC",    # DpnI site
        "GGCC",    # HaeIII site
        "TCGA",    # TaqI site
    ]
    return random.choice(flanking_options)

def region_inversion_at_restriction(sequence, min_size=100, max_size=2000):
    """Invert a region starting/ending at restriction sites"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        # Fallback to random region if no restriction sites
        seq_len = len(sequence)
        invert_size = random.randint(min_size, min(max_size, seq_len // 10))
        start = random.randint(0, seq_len - invert_size)
        end = start + invert_size
        restriction_info = {"method": "random_region", "enzyme": None}
    else:
        site_len = len(RESTRICTION_SITES[enzyme])
        start = site_pos
        # Add biological randomness to size
        size_variation = random.randint(min_size, min(max_size, len(sequence) // 20))
        end = min(start + size_variation, len(sequence))
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    # Reverse complement the region
    region = sequence[start:end]
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    inverted = ''.join(complement.get(base, base) for base in reversed(region))
    
    modified = sequence[:start] + inverted + sequence[end:]
    
    return modified, {
        "method": "region_inversion",
        "start": start,
        "end": end,
        "inverted_length": end - start,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info
    }

def region_duplication_at_restriction(sequence, min_size=200, max_size=3000):
    """Duplicate region starting at restriction site"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        # Fallback to random region
        seq_len = len(sequence)
        dup_size = random.randint(min_size, min(max_size, seq_len // 20))
        source_start = random.randint(0, seq_len - dup_size)
        source_end = source_start + dup_size
        restriction_info = {"method": "random_region", "enzyme": None}
    else:
        site_len = len(RESTRICTION_SITES[enzyme])
        source_start = site_pos
        size_variation = random.randint(min_size, min(max_size, len(sequence) // 20))
        source_end = min(source_start + size_variation, len(sequence))
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    region_to_dup = sequence[source_start:source_end]
    
    # Add biological randomness - insert at another restriction site
    insert_enzyme, insert_pos = get_random_restriction_site(sequence)
    if insert_pos is None:
        insert_pos = random.randint(0, len(sequence))
        insert_restriction_info = {"method": "random_position", "enzyme": None}
    else:
        insert_restriction_info = {"method": "restriction_site", "enzyme": insert_enzyme, "site": RESTRICTION_SITES[insert_enzyme]}
    
    modified = sequence[:insert_pos] + region_to_dup + sequence[insert_pos:]
    
    return modified, {
        "method": "region_duplication",
        "source_start": source_start,
        "source_end": source_end,
        "insert_position": insert_pos,
        "duplicated_length": source_end - source_start,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "source_restriction": restriction_info,
        "insert_restriction": insert_restriction_info
    }

def region_deletion_at_restriction(sequence, min_size=50, max_size=1000):
    """Delete region starting at restriction site"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        # Fallback to random region
        seq_len = len(sequence)
        del_size = random.randint(min_size, min(max_size, seq_len // 50))
        start = random.randint(0, seq_len - del_size)
        end = start + del_size
        restriction_info = {"method": "random_region", "enzyme": None}
    else:
        site_len = len(RESTRICTION_SITES[enzyme])
        start = site_pos
        size_variation = random.randint(min_size, min(max_size, len(sequence) // 50))
        end = min(start + size_variation, len(sequence))
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    modified = sequence[:start] + sequence[end:]
    
    return modified, {
        "method": "region_deletion",
        "deleted_start": start,
        "deleted_end": end,
        "deleted_length": end - start,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info
    }

def insert_gfp_at_restriction(sequence):
    """Insert GFP at restriction site with biological randomness"""
    gfp = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
    
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        insert_pos = random.randint(0, len(sequence))
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        insert_pos = site_pos
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    # Add biological randomness
    flanking = add_biological_randomness()
    construct = flanking + gfp + flanking
    
    modified = sequence[:insert_pos] + construct + sequence[insert_pos:]
    
    return modified, {
        "method": "gfp_insertion",
        "insert_position": insert_pos,
        "inserted_sequence": "GFP_with_flanking",
        "inserted_length": len(construct),
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info,
        "flanking_sequence": flanking
    }

def insert_tn5_at_restriction(sequence):
    """Insert Tn5 puromycin marker at restriction site"""
    tn5_puro = "ATGACCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCCGCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAGCTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCCGCGGTGGCGGTCTGGACCACGCCGGAGAGCGTCGAAGCGGGGGCGGTGTTCGCCGAGATCGGCCCGCGCATGGCCGAGTTGAGCGGTTCCCGGCTGGCCGCGCAGCAACAGATGGAAGGCCTCCTGGCGCCGCACCGGCCCAAGGAGCCCGCGTGGTTCCTGGCCACCGTCGGCGTCTCGCCCGACCACCAGGGCAAGGGTCTGGGCAGCGCCGTCGTGCTCCCCGGAGTGGAGGCGGCCGAGCGCGCCGGGGTGCCCGCCTTCCTGGAGACCTCCGCGCCCCGCAACCTCCCCTTCTACGAGCGGCTCGGCTTCACCGTCACCGCCGACGTCGAGGTGCCCGAAGGACCGCGCACCTGGTGCATGACCCGCAAGCCCGGTGCCTGA"
    
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        insert_pos = random.randint(0, len(sequence))
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        insert_pos = site_pos
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    flanking = add_biological_randomness()
    construct = flanking + tn5_puro + flanking
    
    modified = sequence[:insert_pos] + construct + sequence[insert_pos:]
    
    return modified, {
        "method": "tn5_insertion",
        "insert_position": insert_pos,
        "inserted_sequence": "Tn5_puromycin_with_flanking",
        "inserted_length": len(construct),
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info,
        "flanking_sequence": flanking
    }

def insert_igem_gfp_at_restriction(sequence):
    """Insert iGEM GFP at restriction site"""
    igem_gfp = "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAATAA"
    
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        insert_pos = random.randint(0, len(sequence))
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        insert_pos = site_pos
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    flanking = add_biological_randomness()
    construct = flanking + igem_gfp + flanking
    
    modified = sequence[:insert_pos] + construct + sequence[insert_pos:]
    
    return modified, {
        "method": "igem_gfp_insertion",
        "insert_position": insert_pos,
        "inserted_sequence": "iGEM_GFP_with_flanking",
        "inserted_length": len(construct),
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info,
        "flanking_sequence": flanking
    }

def frameshift_mutation_at_restriction(sequence):
    """Create frameshift near restriction site"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        pos = random.randint(50, len(sequence) - 50)
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        # Add frameshift near restriction site (biologically realistic)
        offset = random.randint(-50, 50)
        pos = max(0, min(site_pos + offset, len(sequence) - 10))
        restriction_info = {"method": "near_restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme], "offset": offset}
    
    if random.choice([True, False]):
        # Insertion
        bases = ['A', 'T', 'G', 'C']
        insert_length = random.randint(1, 3)  # Biologically realistic sizes
        insertion = ''.join(random.choices(bases, k=insert_length))
        modified = sequence[:pos] + insertion + sequence[pos:]
        operation = f"insert_{insertion}"
    else:
        # Deletion
        del_length = random.randint(1, 3)
        end_pos = min(pos + del_length, len(sequence))
        modified = sequence[:pos] + sequence[end_pos:]
        operation = f"delete_{del_length}bp"
    
    return modified, {
        "method": "frameshift_mutation",
        "position": pos,
        "operation": operation,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info
    }

def single_nucleotide_change_at_restriction(sequence):
    """Create SNP near restriction site"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        pos = random.randint(0, len(sequence) - 1)
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        # Add SNP near restriction site
        offset = random.randint(-20, 20)
        pos = max(0, min(site_pos + offset, len(sequence) - 1))
        restriction_info = {"method": "near_restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme], "offset": offset}
    
    original_base = sequence[pos]
    bases = ['A', 'T', 'G', 'C']
    new_base = random.choice([b for b in bases if b != original_base])
    
    modified = sequence[:pos] + new_base + sequence[pos + 1:]
    
    return modified, {
        "method": "single_nucleotide_change",
        "position": pos,
        "original_base": original_base,
        "new_base": new_base,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info
    }

def random_substitution_at_restriction(sequence, min_size=50, max_size=500):
    """Substitute region near restriction site"""
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        seq_len = len(sequence)
        sub_size = random.randint(min_size, min(max_size, seq_len // 100))
        start = random.randint(0, seq_len - sub_size)
        end = start + sub_size
        restriction_info = {"method": "random_region", "enzyme": None}
    else:
        # Substitute region starting at restriction site
        site_len = len(RESTRICTION_SITES[enzyme])
        start = site_pos
        sub_size = random.randint(min_size, min(max_size, len(sequence) // 100))
        end = min(start + sub_size, len(sequence))
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    bases = ['A', 'T', 'G', 'C']
    replacement = ''.join(random.choices(bases, k=end - start))
    
    modified = sequence[:start] + replacement + sequence[end:]
    
    return modified, {
        "method": "random_substitution",
        "substitution_start": start,
        "substitution_end": end,
        "substituted_length": end - start,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info
    }

def insert_crispr_at_restriction(sequence):
    """Insert CRISPR construct at restriction site"""
    gfp = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
    pam_site = "NGG"  # Simplified PAM
    
    enzyme, site_pos = get_random_restriction_site(sequence)
    
    if enzyme is None:
        insert_pos = random.randint(0, len(sequence))
        restriction_info = {"method": "random_position", "enzyme": None}
    else:
        insert_pos = site_pos
        restriction_info = {"method": "restriction_site", "enzyme": enzyme, "site": RESTRICTION_SITES[enzyme]}
    
    flanking = add_biological_randomness()
    construct = flanking + gfp + pam_site + flanking
    modified = sequence[:insert_pos] + construct + sequence[insert_pos:]
    
    return modified, {
        "method": "crispr_insertion",
        "insert_position": insert_pos,
        "inserted_sequence": "GFP_with_PAM_and_flanking",
        "inserted_length": len(construct),
        "original_length": len(sequence),
        "modified_length": len(modified),
        "restriction_info": restriction_info,
        "flanking_sequence": flanking
    }

# Enhanced engineering methods using restriction sites
RESTRICTION_ENGINEERING_METHODS = [
    region_inversion_at_restriction,
    region_duplication_at_restriction, 
    region_deletion_at_restriction,
    insert_gfp_at_restriction,
    insert_tn5_at_restriction,
    insert_igem_gfp_at_restriction,
    frameshift_mutation_at_restriction,
    single_nucleotide_change_at_restriction,
    random_substitution_at_restriction,
    insert_crispr_at_restriction
]

def engineer_virus_with_restriction_sites(virus_key, method_func):
    """Engineer a virus using specified method with restriction sites"""
    virus_data = VIRUSES[virus_key]
    sequence = get_clean_sequence(virus_data["sequence"])
    
    try:
        modified_seq, details = method_func(sequence)
        details["virus"] = virus_data["name"]
        details["virus_key"] = virus_key
        details["virus_family"] = virus_data["family"]
        details["full_genome"] = modified_seq
        
        return details
    except Exception as e:
        # Fallback to simple insertion if method fails
        return engineer_virus_simple_fallback(virus_key, sequence)

def engineer_virus_simple_fallback(virus_key, sequence):
    """Simple fallback engineering"""
    virus_data = VIRUSES[virus_key]
    insert_pos = random.randint(0, len(sequence))
    gfp = "ATGGTGAGCAAGGGCGAGGAG"  # Short GFP fragment
    modified = sequence[:insert_pos] + gfp + sequence[insert_pos:]
    
    return {
        "virus": virus_data["name"],
        "virus_key": virus_key,
        "method": "fallback_insertion",
        "insert_position": insert_pos,
        "original_length": len(sequence),
        "modified_length": len(modified),
        "full_genome": modified,
        "virus_family": virus_data["family"],
        "restriction_info": {"method": "fallback", "enzyme": None}
    }

def generate_restriction_based_dataset(num_per_virus=1000):
    """Generate dataset using restriction sites"""
    print(f"🧬 Generating Restriction-Based Virus Engineering Dataset")
    print(f"📊 Target: {num_per_virus} examples per virus")
    print("=" * 60)
    
    results = []
    virus_keys = list(VIRUSES.keys())
    methods_per_virus = len(RESTRICTION_ENGINEERING_METHODS)
    examples_per_method = num_per_virus // methods_per_virus
    
    print(f"Methods available: {methods_per_virus}")
    print(f"Examples per method: {examples_per_method}")
    print()
    
    for virus_key in virus_keys:
        virus_results = []
        print(f"Processing {VIRUSES[virus_key]['name']}...")
        
        for i, method_func in enumerate(RESTRICTION_ENGINEERING_METHODS):
            method_name = method_func.__name__
            print(f"  {method_name}: ", end="")
            
            for j in range(examples_per_method):
                result = engineer_virus_with_restriction_sites(virus_key, method_func)
                result["example_id"] = f"{virus_key}_{method_name}_{j+1}"
                virus_results.append(result)
                
                if (j + 1) % 20 == 0:
                    print(f"{j+1}", end=" ")
            
            print(f"✅ {examples_per_method} examples")
        
        # Add extra examples to reach exactly num_per_virus
        remaining = num_per_virus - len(virus_results)
        if remaining > 0:
            print(f"  Adding {remaining} extra examples...")
            for j in range(remaining):
                method_func = random.choice(RESTRICTION_ENGINEERING_METHODS)
                result = engineer_virus_with_restriction_sites(virus_key, method_func)
                result["example_id"] = f"{virus_key}_extra_{j+1}"
                virus_results.append(result)
        
        results.extend(virus_results)
        print(f"✅ Total for {virus_key}: {len(virus_results)} examples\n")
    
    return results

def save_restriction_results(results, base_filename="restriction_virus_engineered"):
    """Save results in multiple formats"""
    # JSON with all metadata
    json_file = f"{base_filename}.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"💾 JSON results saved to {json_file}")
    
    # Individual FASTA files
    os.makedirs("restriction_genomes", exist_ok=True)
    print(f"💾 Saving individual FASTA files to restriction_genomes/")
    
    for i, result in enumerate(results):
        filename = f"restriction_genomes/{result['example_id']}.fasta"
        
        with open(filename, 'w') as f:
            header = f">{result['virus']} | {result['method']} | {result['original_length']:,} -> {result['modified_length']:,} bp"
            if 'restriction_info' in result and result['restriction_info']['enzyme']:
                header += f" | {result['restriction_info']['enzyme']}"
            f.write(header + "\n")
            
            genome = result['full_genome']
            for j in range(0, len(genome), 80):
                f.write(genome[j:j+80] + "\n")
        
        if (i + 1) % 200 == 0:
            print(f"  Saved {i+1}/{len(results)} files...")
    
    print(f"✅ Saved {len(results)} individual FASTA files")

if __name__ == "__main__":
    print("Available viruses:")
    for key, data in VIRUSES.items():
        print(f"  {key}: {data['name']} ({data['length']:,} bp)")
    print()
    
    print("Available restriction-based engineering methods:")
    for i, method in enumerate(RESTRICTION_ENGINEERING_METHODS, 1):
        print(f"  {i}. {method.__name__}")
    print()
    
    # Generate the restriction-based dataset
    results = generate_restriction_based_dataset(num_per_virus=1000)
    save_restriction_results(results)
    
    print(f"\n🎉 RESTRICTION-BASED DATASET COMPLETE!")
    print(f"✅ Generated {len(results)} total engineered genomes")
    print(f"✅ Uses actual restriction sites for realistic editing")
    print(f"✅ Biological randomness added (flanking sequences, size variations)")
    print(f"✅ Full genomes with complete metadata")
