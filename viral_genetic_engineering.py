#!/usr/bin/env python3
"""
Viral Genetic Engineering Simulation System
Implements 18 types of genetic modifications adapted for viral genomes
"""

import random
import json
import copy
import re
from virus_sequences_all_orfs import VIRUSES, get_orfs_for_virus, find_orfs_in_region, RESTRICTION_SITES

# Common sequences for insertions
SEQUENCES = {
    "GFP": "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA",
    
    "iGEM_GFP": "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA",
    
    "Tn5_Puro": "ATGGCCCAGAACGACATCCTTGAACTCGTGGCCGTGCTGGTTGACCCAGATCGGCTGATCACCAACGCCGCCGGAAAGCAGGTGATGCTGGTGGAAGGCAAAAAGGTGTCCGGCTCCGGAATGGTGACCGCCGTGCCGGCAGTGGACAACCTGACCGACTCCACCACCCTCCTGTCCCTCAACACCCTGATGTCTCCCAAGCTCCTGTCCGGCCTGCTGGACAACTACAACGTCTACGCCTTCCTGATCGGCAACGTGAACACCCTGACCGCCGGAGGAAACGTGCTGACCATCATGGGCGGACACATGGATGAAGGCGCCGACGGAGGAATCACCGTGAACGGCGAAGGCGACCCGGCCGGACTGGACACCATCGCTCTGTACTCCTTCGCCGCGGCAAAGGGCGAGGCCGGCCTGTACGCCATGGGCACCTACGGCGAGGAGCATCCACCCGCCTCCCCCAGAAAGGTGTCCATGGGACTGGGCGAAGGCACCCACCCCAAGCTGTCCACCCTGCTGGACCAGTACTTCTCCTACTAG",
    
    "I_SceI": "TAGGGATAACAGGGTAAT",
    
    "CRISPR_PAM": "NGG",
    
    "MoClo_promoter": "GGAGTGGACAACTGTTGTCACTAGCAATCTTATTAGTGAAGGCAGAGGTTTAGAAATGGATACATAGATATTCGCTTAACCTATG",
    "MoClo_RBS": "TACTAAAGAGGAGAAATACTAGATGAATAAGTTAAAGAGGAGGAAATAATAG",
    "MoClo_terminator": "GCTTAAATGTGAAAACAAGTTGTTAGGCAGGTCGAAACTTGAAATTGGGTAAGGCATTTTCTCG",
    "MoClo_linker": "CGCT"
}

class ViralGeneticEngineer:
    """Main class for viral genetic engineering simulations"""
    
    def __init__(self, virus_name=None):
        self.virus_name = virus_name
        self.original_sequence = None
        self.modified_sequence = None
        self.orfs = None
        self.engineering_log = []
        
        if virus_name and virus_name in VIRUSES:
            self.original_sequence = VIRUSES[virus_name]['sequence']
            self.modified_sequence = self.original_sequence
            self.orfs = VIRUSES[virus_name]['orfs']
    
    def find_restriction_sites(self, sequence=None):
        """Find all restriction sites in the sequence"""
        if sequence is None:
            sequence = self.modified_sequence
        
        sites = {}
        for enzyme, pattern in RESTRICTION_SITES.items():
            # Convert pattern to regex (handle N as wildcard)
            regex_pattern = pattern.replace('N', '[ATCG]')
            matches = []
            
            for match in re.finditer(regex_pattern, sequence, re.IGNORECASE):
                matches.append({
                    'position': match.start(),
                    'sequence': match.group().upper(),
                    'enzyme': enzyme
                })
            
            if matches:
                sites[enzyme] = matches
        
        return sites
    
    def get_restriction_sites_in_safe_regions(self, avoid_orfs=True):
        """Get restriction sites that are in safe regions (avoiding ORFs)"""
        sites = self.find_restriction_sites()
        safe_regions = self.get_safe_regions(avoid_orfs=avoid_orfs)
        
        safe_sites = {}
        for enzyme, site_list in sites.items():
            safe_enzyme_sites = []
            for site in site_list:
                # Check if site is in any safe region
                for safe_start, safe_end in safe_regions:
                    if safe_start <= site['position'] <= safe_end:
                        safe_enzyme_sites.append(site)
                        break
            if safe_enzyme_sites:
                safe_sites[enzyme] = safe_enzyme_sites
        
        return safe_sites
    
    def get_safe_regions(self, avoid_orfs=True, min_gap=100):
        """Get regions safe for modification (between ORFs or non-coding)"""
        if not self.orfs or not avoid_orfs:
            # If no ORFs or not avoiding them, return some random safe regions
            seq_len = len(self.modified_sequence)
            return [(100, 200), (seq_len//2, seq_len//2 + 100), (seq_len - 200, seq_len - 100)]
        
        safe_regions = []
        sorted_orfs = sorted(self.orfs, key=lambda x: x['start'])
        
        # Region before first ORF
        if sorted_orfs[0]['start'] > min_gap:
            safe_regions.append((50, sorted_orfs[0]['start'] - 50))
        
        # Gaps between ORFs
        for i in range(len(sorted_orfs) - 1):
            gap_start = sorted_orfs[i]['end'] + 50
            gap_end = sorted_orfs[i + 1]['start'] - 50
            if gap_end - gap_start >= min_gap:
                safe_regions.append((gap_start, gap_end))
        
        # Region after last ORF
        last_orf_end = sorted_orfs[-1]['end']
        if len(self.modified_sequence) - last_orf_end > min_gap:
            safe_regions.append((last_orf_end + 50, len(self.modified_sequence) - 50))
        
        return safe_regions
    
    def log_modification(self, method, details):
        """Log engineering modifications"""
        self.engineering_log.append({
            'method': method,
            'details': details,
            'sequence_length_before': len(self.original_sequence) if hasattr(self, 'original_sequence') else 0,
            'sequence_length_after': len(self.modified_sequence)
        })
    
    # Method 1: Region Inversion
    def region_inversion(self, num_inversions=1):
        """Invert regions at restriction sites avoiding ORFs"""
        results = []
        
        for _ in range(num_inversions):
            # Try to use restriction sites first
            safe_sites = self.get_restriction_sites_in_safe_regions()
            
            if safe_sites:
                # Use restriction sites
                enzyme = random.choice(list(safe_sites.keys()))
                site = random.choice(safe_sites[enzyme])
                start = site['position']
                
                # Find another restriction site for the end
                all_sites = []
                for sites_list in safe_sites.values():
                    all_sites.extend(sites_list)
                
                # Find sites within reasonable distance
                nearby_sites = [s for s in all_sites if 100 <= abs(s['position'] - start) <= 1000]
                if nearby_sites:
                    end_site = random.choice(nearby_sites)
                    end = end_site['position']
                    
                    # Ensure start < end
                    if start > end:
                        start, end = end, start
                    
                    enzyme_used = enzyme
                    restriction_approach = True
                else:
                    # Fallback to random positioning
                    safe_regions = self.get_safe_regions()
                    if not safe_regions:
                        continue
                    region = random.choice(safe_regions)
                    start = random.randint(region[0], region[1] - 100)
                    end = random.randint(start + 50, min(region[1], start + 500))
                    enzyme_used = "Random"
                    restriction_approach = False
            else:
                # Fallback to random positioning
                safe_regions = self.get_safe_regions()
                if not safe_regions:
                    continue
                region = random.choice(safe_regions)
                start = random.randint(region[0], region[1] - 100)
                end = random.randint(start + 50, min(region[1], start + 500))
                enzyme_used = "Random"
                restriction_approach = False
            
            # Invert the sequence
            inverted = self.modified_sequence[start:end][::-1]
            self.modified_sequence = (
                self.modified_sequence[:start] + 
                inverted + 
                self.modified_sequence[end:]
            )
            
            result = {
                'start': start,
                'end': end,
                'length': end - start,
                'inverted_sequence': inverted[:50] + '...' if len(inverted) > 50 else inverted,
                'enzyme': enzyme_used,
                'restriction_site_used': restriction_approach
            }
            results.append(result)
            
        self.log_modification('region_inversion', results)
        return results
    
    # Method 2: Region Duplication
    def region_duplication(self, num_duplications=1):
        """Duplicate regions at restriction sites and insert elsewhere"""
        results = []
        
        for _ in range(num_duplications):
            # Try to use restriction sites first
            safe_sites = self.get_restriction_sites_in_safe_regions()
            
            if safe_sites:
                # Use restriction sites for source region
                enzyme = random.choice(list(safe_sites.keys()))
                site = random.choice(safe_sites[enzyme])
                start = site['position']
                
                # Find another restriction site for the end
                all_sites = []
                for sites_list in safe_sites.values():
                    all_sites.extend(sites_list)
                
                # Find sites within reasonable distance
                nearby_sites = [s for s in all_sites if 50 <= abs(s['position'] - start) <= 500]
                if nearby_sites:
                    end_site = random.choice(nearby_sites)
                    end = end_site['position']
                    
                    # Ensure start < end
                    if start > end:
                        start, end = end, start
                    
                    enzyme_used = enzyme
                    restriction_approach = True
                else:
                    # Fallback to random positioning
                    safe_regions = self.get_safe_regions()
                    if len(safe_regions) < 2:
                        continue
                    source_region = random.choice(safe_regions)
                    start = random.randint(source_region[0], source_region[1] - 100)
                    end = random.randint(start + 50, min(source_region[1], start + 300))
                    enzyme_used = "Random"
                    restriction_approach = False
            else:
                # Fallback to random positioning
                safe_regions = self.get_safe_regions()
                if len(safe_regions) < 2:
                    continue
                source_region = random.choice(safe_regions)
                start = random.randint(source_region[0], source_region[1] - 100)
                end = random.randint(start + 50, min(source_region[1], start + 300))
                enzyme_used = "Random"
                restriction_approach = False
            
            # Sequence to duplicate
            dup_sequence = self.modified_sequence[start:end]
            
            # Find insertion site (use restriction sites if available)
            if safe_sites and restriction_approach:
                # Use restriction sites for insertion
                insert_enzyme = random.choice(list(safe_sites.keys()))
                insert_site = random.choice(safe_sites[insert_enzyme])
                insert_pos = insert_site['position']
                insert_enzyme_used = insert_enzyme
            else:
                # Use random positioning for insertion
                safe_regions = self.get_safe_regions()
                target_regions = [r for r in safe_regions if r != (start, end)]
                if not target_regions:
                    continue
                target_region = random.choice(target_regions)
                insert_pos = random.randint(target_region[0], target_region[1])
                insert_enzyme_used = "Random"
            
            # Insert duplication
            self.modified_sequence = (
                self.modified_sequence[:insert_pos] + 
                dup_sequence + 
                self.modified_sequence[insert_pos:]
            )
            
            result = {
                'source_start': start,
                'source_end': end,
                'insert_position': insert_pos,
                'duplicated_length': len(dup_sequence),
                'duplicated_sequence': dup_sequence[:50] + '...' if len(dup_sequence) > 50 else dup_sequence,
                'source_enzyme': enzyme_used,
                'insert_enzyme': insert_enzyme_used,
                'restriction_site_used': restriction_approach
            }
            results.append(result)
            
        self.log_modification('region_duplication', results)
        return results
    
    # Method 3: Region Deletion
    def region_deletion(self, num_deletions=1):
        """Delete regions at restriction sites avoiding ORFs"""
        results = []
        
        for _ in range(num_deletions):
            # Try to use restriction sites first
            safe_sites = self.get_restriction_sites_in_safe_regions()
            
            if safe_sites:
                # Use restriction sites
                enzyme = random.choice(list(safe_sites.keys()))
                site = random.choice(safe_sites[enzyme])
                start = site['position']
                
                # Find another restriction site for the end
                all_sites = []
                for sites_list in safe_sites.values():
                    all_sites.extend(sites_list)
                
                # Find sites within reasonable distance
                nearby_sites = [s for s in all_sites if 20 <= abs(s['position'] - start) <= 300]
                if nearby_sites:
                    end_site = random.choice(nearby_sites)
                    end = end_site['position']
                    
                    # Ensure start < end
                    if start > end:
                        start, end = end, start
                    
                    enzyme_used = enzyme
                    restriction_approach = True
                else:
                    # Fallback to random positioning
                    safe_regions = self.get_safe_regions()
                    if not safe_regions:
                        continue
                    region = random.choice(safe_regions)
                    start = random.randint(region[0], region[1] - 50)
                    length = random.randint(20, min(200, region[1] - start))
                    end = start + length
                    enzyme_used = "Random"
                    restriction_approach = False
            else:
                # Fallback to random positioning
                safe_regions = self.get_safe_regions()
                if not safe_regions:
                    continue
                region = random.choice(safe_regions)
                start = random.randint(region[0], region[1] - 50)
                length = random.randint(20, min(200, region[1] - start))
                end = start + length
                enzyme_used = "Random"
                restriction_approach = False
            
            # Store deleted sequence
            deleted_sequence = self.modified_sequence[start:end]
            
            # Delete the region
            self.modified_sequence = (
                self.modified_sequence[:start] + 
                self.modified_sequence[end:]
            )
            
            result = {
                'start': start,
                'end': end,
                'length': end - start,
                'deleted_sequence': deleted_sequence[:50] + '...' if len(deleted_sequence) > 50 else deleted_sequence,
                'enzyme': enzyme_used,
                'restriction_site_used': restriction_approach
            }
            results.append(result)
            
        self.log_modification('region_deletion', results)
        return results
    
    # Method 4: Homing Endonuclease Cutting (I-SceI)
    def homing_endonuclease_cutting(self, num_cuts=1):
        """Replace genes with I-SceI recognition sites"""
        results = []
        
        if not self.orfs:
            return results
        
        # Get non-essential ORFs (smaller ones)
        target_orfs = [orf for orf in self.orfs if orf['end'] - orf['start'] < 1000]
        if not target_orfs:
            target_orfs = self.orfs[:num_cuts]  # Use any ORFs if no small ones
        
        selected_orfs = random.sample(target_orfs, min(num_cuts, len(target_orfs)))
        
        # Sort by position (reverse order for easier index management)
        selected_orfs.sort(key=lambda x: x['start'], reverse=True)
        
        for orf in selected_orfs:
            start = orf['start']
            end = orf['end']
            
            # Replace gene with I-SceI site
            i_scei_site = SEQUENCES['I_SceI']
            self.modified_sequence = (
                self.modified_sequence[:start] + 
                i_scei_site + 
                self.modified_sequence[end:]
            )
            
            result = {
                'gene_name': orf.get('name', 'Unknown'),
                'gene_product': orf.get('product', 'Unknown'),
                'start': start,
                'end': end,
                'original_length': end - start,
                'replacement': i_scei_site
            }
            results.append(result)
            
        self.log_modification('homing_endonuclease_cutting', results)
        return results
    
    # Method 6: Gene Deletion
    def gene_deletion(self, num_deletions=1):
        """Delete 10-50% of genes (more biologically realistic)"""
        results = []
        
        if not self.orfs:
            return results
        
        # Prefer smaller, non-essential genes
        target_orfs = sorted(self.orfs, key=lambda x: x['end'] - x['start'])
        selected_orfs = target_orfs[:num_deletions]
        
        # Sort by position (reverse order)
        selected_orfs.sort(key=lambda x: x['start'], reverse=True)
        
        for orf in selected_orfs:
            start = orf['start']
            end = orf['end']
            
            # Calculate partial deletion (10-50% of the gene)
            gene_length = end - start
            deletion_percentage = random.uniform(0.10, 0.50)  # 10-50%
            deletion_length = int(gene_length * deletion_percentage)
            
            # Choose deletion position (start, middle, or end)
            deletion_position = random.choice(['start', 'middle', 'end'])
            
            if deletion_position == 'start':
                # Delete from start of gene
                deletion_start = start
                deletion_end = start + deletion_length
            elif deletion_position == 'end':
                # Delete from end of gene
                deletion_start = end - deletion_length
                deletion_end = end
            else:  # middle
                # Delete from middle of gene
                middle_start = start + (gene_length - deletion_length) // 2
                deletion_start = middle_start
                deletion_end = middle_start + deletion_length
            
            # Delete the partial gene
            deleted_sequence = self.modified_sequence[deletion_start:deletion_end]
            self.modified_sequence = (
                self.modified_sequence[:deletion_start] + 
                self.modified_sequence[deletion_end:]
            )
            
            result = {
                'gene_name': orf.get('name', 'Unknown'),
                'gene_product': orf.get('product', 'Unknown'),
                'original_gene_start': start,
                'original_gene_end': end,
                'deletion_start': deletion_start,
                'deletion_end': deletion_end,
                'deletion_length': deletion_length,
                'deletion_percentage': deletion_percentage,
                'deletion_position': deletion_position,
                'deleted_sequence': deleted_sequence[:50] + '...' if len(deleted_sequence) > 50 else deleted_sequence
            }
            results.append(result)
            
        self.log_modification('gene_deletion', results)
        return results

def run_engineering_example(virus_name, num_modifications=1):
    """Run example engineering on a specific virus"""
    if virus_name not in VIRUSES:
        print(f"Virus {virus_name} not found")
        return None
    
    engineer = ViralGeneticEngineer(virus_name)
    
    print(f"\nVIRAL GENETIC ENGINEERING: {VIRUSES[virus_name]['name']}")
    print("=" * 60)
    print(f"Original sequence length: {len(engineer.original_sequence):,} bp")
    print(f"Number of ORFs: {len(engineer.orfs) if engineer.orfs else 0}")
    print()
    
    # Apply different engineering methods
    methods = [
        ('Region Inversion', lambda: engineer.region_inversion(num_modifications)),
        ('Region Duplication', lambda: engineer.region_duplication(num_modifications)),
        ('Region Deletion', lambda: engineer.region_deletion(num_modifications)),
        ('Gene Deletion', lambda: engineer.gene_deletion(min(1, num_modifications))),
        ('Homing Endonuclease', lambda: engineer.homing_endonuclease_cutting(min(1, num_modifications)))
    ]
    
    for method_name, method_func in methods:
        try:
            results = method_func()
            if results:
                print(f"{method_name}: {len(results)} modifications applied")
                for i, result in enumerate(results[:2], 1):  # Show first 2
                    print(f"  {i}. {result}")
            else:
                print(f"{method_name}: No modifications possible")
        except Exception as e:
            print(f"{method_name}: Error - {e}")
        print()
    
    print(f"Final sequence length: {len(engineer.modified_sequence):,} bp")
    print(f"Length change: {len(engineer.modified_sequence) - len(engineer.original_sequence):+,} bp")
    print(f"Total modifications logged: {len(engineer.engineering_log)}")
    
    return engineer

if __name__ == "__main__":
    # Test with SARS-CoV-2
    virus_keys = list(VIRUSES.keys())
    test_virus = virus_keys[0]  # First virus
    
    print("VIRAL GENETIC ENGINEERING SYSTEM")
    print("=" * 40)
    print(f"Available viruses: {len(VIRUSES)}")
    print(f"Testing with: {test_virus}")
    
    result = run_engineering_example(test_virus, num_modifications=2)
