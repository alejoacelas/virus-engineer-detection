"""
Multi-Virus Engineering System with Realistic Restriction Sites
"""

import random
import json
from virus_sequences import VIRUSES, RESTRICTION_SITES, SAFE_REGIONS

def get_clean_sequence(sequence):
    """Clean sequence string"""
    return sequence.replace('\n', '').replace(' ', '').upper()

def find_restriction_sites(sequence, enzyme="BamHI"):
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

def find_safe_region(sequence, virus_family):
    """Find a safe region for modification"""
    if virus_family not in SAFE_REGIONS:
        # Default: use last 10% of genome
        return (int(len(sequence) * 0.9), len(sequence))
    
    regions = SAFE_REGIONS[virus_family]
    region_name = random.choice(list(regions.keys()))
    return regions[region_name]

def insert_at_restriction_site(sequence, virus_family, enzyme="BamHI"):
    """Insert sequence at restriction site"""
    sites = find_restriction_sites(sequence, enzyme)
    
    if not sites:
        # No restriction sites found, use safe region
        safe_start, safe_end = find_safe_region(sequence, virus_family)
        pos = random.randint(safe_start, safe_end - 100)
        insertion_site = "safe_region"
    else:
        # Use restriction site
        pos = random.choice(sites)
        insertion_site = f"{enzyme}_site"
    
    # GFP sequence
    gfp = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
    
    modified = sequence[:pos] + gfp + sequence[pos:]
    
    return modified, {
        "position": pos,
        "site_type": insertion_site,
        "enzyme": enzyme if insertion_site.endswith("_site") else None,
        "inserted_length": len(gfp)
    }

def delete_via_restriction_sites(sequence, virus_family):
    """Delete sequence between two restriction sites"""
    # Try different enzymes
    enzymes = list(RESTRICTION_SITES.keys())
    random.shuffle(enzymes)
    
    for enzyme in enzymes:
        sites = find_restriction_sites(sequence, enzyme)
        if len(sites) >= 2:
            # Use two restriction sites
            site1, site2 = sorted(random.sample(sites, 2))
            cut_length = site2 - site1
            if 50 <= cut_length <= 500:  # Reasonable deletion size
                modified = sequence[:site1] + sequence[site2:]
                return modified, {
                    "deletion_start": site1,
                    "deletion_end": site2,
                    "deleted_length": cut_length,
                    "enzyme": enzyme,
                    "site_type": f"between_{enzyme}_sites"
                }
    
    # Fallback: delete in safe region
    safe_start, safe_end = find_safe_region(sequence, virus_family)
    delete_length = random.randint(50, 200)
    start = random.randint(safe_start, safe_end - delete_length)
    end = start + delete_length
    
    modified = sequence[:start] + sequence[end:]
    return modified, {
        "deletion_start": start,
        "deletion_end": end,
        "deleted_length": delete_length,
        "site_type": "safe_region"
    }

def substitute_via_restriction_sites(sequence, virus_family):
    """Substitute sequence using restriction sites"""
    enzymes = list(RESTRICTION_SITES.keys())
    random.shuffle(enzymes)
    
    for enzyme in enzymes:
        sites = find_restriction_sites(sequence, enzyme)
        if len(sites) >= 2:
            site1, site2 = sorted(random.sample(sites, 2))
            cut_length = site2 - site1
            if 30 <= cut_length <= 300:
                # Generate replacement sequence
                bases = ['A', 'T', 'G', 'C']
                replacement = ''.join(random.choices(bases, k=cut_length))
                
                modified = sequence[:site1] + replacement + sequence[site2:]
                return modified, {
                    "substitution_start": site1,
                    "substitution_end": site2,
                    "original_length": cut_length,
                    "replacement_length": len(replacement),
                    "enzyme": enzyme,
                    "site_type": f"between_{enzyme}_sites"
                }
    
    # Fallback: substitute in safe region
    safe_start, safe_end = find_safe_region(sequence, virus_family)
    sub_length = random.randint(30, 150)
    start = random.randint(safe_start, safe_end - sub_length)
    end = start + sub_length
    
    bases = ['A', 'T', 'G', 'C']
    replacement = ''.join(random.choices(bases, k=sub_length))
    
    modified = sequence[:start] + replacement + sequence[end:]
    return modified, {
        "substitution_start": start,
        "substitution_end": end,
        "original_length": sub_length,
        "replacement_length": len(replacement),
        "site_type": "safe_region"
    }

def engineer_virus(virus_key, operation, enzyme=None):
    """Engineer a specific virus"""
    virus_data = VIRUSES[virus_key]
    sequence = get_clean_sequence(virus_data["sequence"])
    
    if operation == "insert_gfp":
        if not enzyme:
            enzyme = random.choice(list(RESTRICTION_SITES.keys()))
        modified_seq, details = insert_at_restriction_site(sequence, virus_data["family"], enzyme)
    elif operation == "delete_region":
        modified_seq, details = delete_via_restriction_sites(sequence, virus_data["family"])
    elif operation == "substitute_region":
        modified_seq, details = substitute_via_restriction_sites(sequence, virus_data["family"])
    else:
        raise ValueError(f"Unknown operation: {operation}")
    
    return {
        "virus": virus_data["name"],
        "virus_key": virus_key,
        "operation": operation,
        "original_length": len(sequence),
        "modified_length": len(modified_seq),
        "full_genome": modified_seq,
        "modification_details": details,
        "virus_family": virus_data["family"]
    }

def run_engineering(num_examples=30):
    """Run multi-virus engineering"""
    print("🧬 Multi-Virus Engineering with Restriction Sites")
    print("=" * 60)
    
    results = []
    virus_keys = list(VIRUSES.keys())
    operations = ["insert_gfp", "delete_region", "substitute_region"]
    enzymes = list(RESTRICTION_SITES.keys())
    
    for i in range(num_examples):
        virus = random.choice(virus_keys)
        operation = random.choice(operations)
        enzyme = random.choice(enzymes) if operation == "insert_gfp" else None
        
        result = engineer_virus(virus, operation, enzyme)
        results.append(result)
        
        print(f"Example {i+1}: {result['virus']} - {operation}")
        print(f"  Length: {result['original_length']:,} -> {result['modified_length']:,} bp")
        print(f"  Site: {result['modification_details']['site_type']}")
        if result['modification_details'].get('enzyme'):
            print(f"  Enzyme: {result['modification_details']['enzyme']}")
        print()
    
    return results

def save_results(results, filename="multi_virus_engineered.json"):
    """Save results to JSON"""
    with open(filename, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"💾 Results saved to {filename}")

def save_individual_genomes(results, output_dir="engineered_genomes"):
    """Save each genome as separate FASTA file"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"💾 Saving individual genomes to {output_dir}/")
    
    for result in results:
        filename = f"{output_dir}/{result['virus_key']}_{result['operation']}_{result['modification_details']['site_type']}.fasta"
        
        with open(filename, 'w') as f:
            # FASTA header with details
            header = f">{result['virus']} | {result['operation']} | {result['original_length']:,} -> {result['modified_length']:,} bp"
            if result['modification_details'].get('enzyme'):
                header += f" | {result['modification_details']['enzyme']}"
            f.write(header + "\n")
            
            # Full genome sequence (80 chars per line)
            genome = result['full_genome']
            for i in range(0, len(genome), 80):
                f.write(genome[i:i+80] + "\n")
        
        print(f"✅ {filename}")

if __name__ == "__main__":
    print("Available viruses:")
    for key, data in VIRUSES.items():
        print(f"  {key}: {data['name']} ({data['length']:,} bp)")
    print()
    
    print("Available restriction enzymes:")
    for enzyme, site in RESTRICTION_SITES.items():
        print(f"  {enzyme}: {site}")
    print()
    
    results = run_engineering(30)
    save_results(results)
    save_individual_genomes(results)
    
    print(f"\n✅ Generated {len(results)} engineered genomes")
    print("✅ Full genomes saved as individual FASTA files")
    print("✅ Realistic engineering at restriction sites")
