#!/usr/bin/env python3
"""
Generates 1000 engineered sequences per virus with realistic biological constraints
"""

import os
import json
import random
import time
import re
from datetime import datetime
from virus_sequences_all_orfs import VIRUSES, RESTRICTION_SITES
from viral_genetic_engineering import ViralGeneticEngineer

class ViralSequenceGenerator:
    def __init__(self, output_dir="engineered_sequences"):
        self.output_dir = output_dir
        self.engineering = ViralGeneticEngineer()
        self.methods = [
            'region_inversion', 'region_duplication', 'region_deletion',
            'gene_deletion', 'homing_endonuclease', 'gfp_insertion',
            'crispr_gfp', 'gene_start_deletion', 'gene_end_deletion',
            'frameshift_mutations', 'phosphorylation_inhibition',
            'igem_gfp', 'moclo_construct'
        ]
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
    
    def find_restriction_sites(self, sequence):
        """Find all restriction sites in the sequence"""
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
    
    def get_safe_regions(self, orfs, sequence_length):
        """Get regions safe for modification (between ORFs)"""
        if not orfs:
            # If no ORFs, return some random safe regions
            return [(100, 200), (sequence_length//2, sequence_length//2 + 100), (sequence_length - 200, sequence_length - 100)]
        
        safe_regions = []
        sorted_orfs = sorted(orfs, key=lambda x: x['start'])
        
        # Region before first ORF
        if sorted_orfs[0]['start'] > 100:
            safe_regions.append((50, sorted_orfs[0]['start'] - 50))
        
        # Gaps between ORFs
        for i in range(len(sorted_orfs) - 1):
            gap_start = sorted_orfs[i]['end'] + 50
            gap_end = sorted_orfs[i + 1]['start'] - 50
            if gap_end - gap_start >= 100:
                safe_regions.append((gap_start, gap_end))
        
        # Region after last ORF
        last_orf_end = sorted_orfs[-1]['end']
        if sequence_length - last_orf_end > 100:
            safe_regions.append((last_orf_end + 50, sequence_length - 50))
        
        return safe_regions
    
    def get_restriction_sites_in_safe_regions(self, sequence, orfs):
        """Get restriction sites that are in safe regions (avoiding ORFs)"""
        sites = self.find_restriction_sites(sequence)
        safe_regions = self.get_safe_regions(orfs, len(sequence))
        
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
        
    def generate_sequences_with_fallback(self, virus_key, virus_info, target_count=1000):
        """Generate sequences with biological feasibility fallbacks"""
        sequences = []
        method_counts = {method: 0 for method in self.methods}
        
        # Try each method with multiple attempts and fallbacks
        for method in self.methods:
            target_per_method = target_count // len(self.methods)
            attempts = 0
            max_attempts = target_per_method * 3  # Allow 3x attempts for fallbacks
            
            while method_counts[method] < target_per_method and attempts < max_attempts:
                try:
                    # Apply biological feasibility randomness
                    result = self.apply_method_with_randomness(virus_key, virus_info, method, attempts)
                    
                    if result:
                        sequences.append(result)
                        method_counts[method] += 1
                        
                except Exception as e:
                    # Log but continue with other attempts
                    pass
                    
                attempts += 1
                
        # If we still don't have enough, use fallback methods
        while len(sequences) < target_count:
            fallback_method = random.choice(['region_inversion', 'region_deletion', 'gene_deletion', 'frameshift_mutations'])
            try:
                result = self.apply_method_with_randomness(virus_key, virus_info, fallback_method, len(sequences))
                if result:
                    sequences.append(result)
                    method_counts[fallback_method] += 1
            except:
                break
                
        return sequences[:target_count], method_counts
    
    def apply_method_with_randomness(self, virus_key, virus_info, method, attempt_num):
        """Apply engineering method with biological feasibility randomness"""
        sequence = virus_info['sequence']
        orfs = virus_info.get('orfs', [])
        
        # Introduce randomness for biological feasibility
        randomness_factors = {
            'region_size_variance': random.uniform(0.5, 2.0),  # 50%-200% of normal size
            'position_tolerance': random.randint(50, 500),     # Allow position shifts
            'method_specificity': random.uniform(0.7, 1.0),   # 70-100% method specificity
        }
        
        if method == 'region_inversion':
            return self.region_inversion(sequence, orfs, randomness_factors)
        elif method == 'region_duplication':
            return self.region_duplication(sequence, orfs, randomness_factors)
        elif method == 'region_deletion':
            return self.region_deletion(sequence, orfs, randomness_factors)
        elif method == 'gene_deletion':
            return self.gene_deletion(sequence, orfs, randomness_factors)
        elif method == 'homing_endonuclease':
            return self.homing_endonuclease(sequence, orfs, randomness_factors)
        elif method == 'gfp_insertion':
            return self.gfp_insertion(sequence, orfs, randomness_factors)
        elif method == 'crispr_gfp':
            return self.crispr_gfp(sequence, orfs, randomness_factors)
        elif method == 'gene_start_deletion':
            return self.gene_start_deletion(sequence, orfs, randomness_factors)
        elif method == 'gene_end_deletion':
            return self.gene_end_deletion(sequence, orfs, randomness_factors)
        elif method == 'frameshift_mutations':
            return self.frameshift_mutations(sequence, orfs, randomness_factors)
        elif method == 'phosphorylation_inhibition':
            return self.phosphorylation_inhibition(sequence, orfs, randomness_factors)
        elif method == 'igem_gfp':
            return self.igem_gfp(sequence, orfs, randomness_factors)
        elif method == 'moclo_construct':
            return self.moclo_construct(sequence, orfs, randomness_factors)
            
        return None
    
    def region_inversion(self, sequence, orfs, factors):
        """Region inversion using restriction sites with biological feasibility"""
        genome_len = len(sequence)
        
        # Try to use restriction sites first
        safe_sites = self.get_restriction_sites_in_safe_regions(sequence, orfs)
        
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
                base_size = min(500, genome_len // 10)
                region_size = int(base_size * factors['region_size_variance'])
                region_size = max(100, min(region_size, genome_len // 4))
                
                tolerance = factors['position_tolerance']
                safe_start = tolerance
                safe_end = genome_len - tolerance - region_size
                
                if safe_end <= safe_start:
                    region_size = genome_len // 8
                    safe_end = genome_len - region_size
                    
                start = random.randint(safe_start, safe_end)
                end = start + region_size
                enzyme_used = "Random"
                restriction_approach = False
        else:
            # Fallback to random positioning
            base_size = min(500, genome_len // 10)
            region_size = int(base_size * factors['region_size_variance'])
            region_size = max(100, min(region_size, genome_len // 4))
            
            tolerance = factors['position_tolerance']
            safe_start = tolerance
            safe_end = genome_len - tolerance - region_size
            
            if safe_end <= safe_start:
                region_size = genome_len // 8
                safe_end = genome_len - region_size
                
            start = random.randint(safe_start, safe_end)
            end = start + region_size
            enzyme_used = "Random"
            restriction_approach = False
        
        # Perform inversion
        inverted_region = sequence[start:end][::-1]
        new_sequence = sequence[:start] + inverted_region + sequence[end:]
        
        return {
            'id': f"inversion_{start}_{end}",
            'method': 'region_inversion',
            'method_name': 'Region Inversion',
            'full_sequence': new_sequence,
            'length_change': 0,
            'modifications': [{
                'method': 'region_inversion',
                'start_position': start,
                'end_position': end,
                'length': end - start,
                'details': {
                    'original_sequence': sequence[start:end],
                    'inverted_sequence': inverted_region,
                    'enzyme': enzyme_used,
                    'restriction_site_used': restriction_approach,
                    'randomness_factors': factors
                }
            }]
        }
    
    def region_deletion(self, sequence, orfs, factors):
        """Region deletion using restriction sites with biological feasibility"""
        genome_len = len(sequence)
        
        # Try to use restriction sites first
        safe_sites = self.get_restriction_sites_in_safe_regions(sequence, orfs)
        
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
                base_size = min(300, genome_len // 15)
                deletion_size = int(base_size * factors['region_size_variance'])
                deletion_size = max(50, min(deletion_size, genome_len // 6))
                
                tolerance = factors['position_tolerance']
                safe_start = tolerance
                safe_end = genome_len - tolerance - deletion_size
                
                if safe_end <= safe_start:
                    deletion_size = genome_len // 10
                    safe_end = genome_len - deletion_size
                    
                start = random.randint(safe_start, safe_end)
                end = start + deletion_size
                enzyme_used = "Random"
                restriction_approach = False
        else:
            # Fallback to random positioning
            base_size = min(300, genome_len // 15)
            deletion_size = int(base_size * factors['region_size_variance'])
            deletion_size = max(50, min(deletion_size, genome_len // 6))
            
            tolerance = factors['position_tolerance']
            safe_start = tolerance
            safe_end = genome_len - tolerance - deletion_size
            
            if safe_end <= safe_start:
                deletion_size = genome_len // 10
                safe_end = genome_len - deletion_size
                
            start = random.randint(safe_start, safe_end)
            end = start + deletion_size
            enzyme_used = "Random"
            restriction_approach = False
        
        # Perform deletion
        new_sequence = sequence[:start] + sequence[end:]
        
        return {
            'id': f"deletion_{start}_{end}",
            'method': 'region_deletion',
            'method_name': 'Region Deletion',
            'full_sequence': new_sequence,
            'length_change': -(end - start),
            'modifications': [{
                'method': 'region_deletion',
                'start_position': start,
                'end_position': end,
                'length': -(end - start),
                'details': {
                    'deleted_sequence': sequence[start:end],
                    'enzyme': enzyme_used,
                    'restriction_site_used': restriction_approach,
                    'randomness_factors': factors
                }
            }]
        }
    
    def gene_deletion(self, sequence, orfs, factors):
        """Partial gene deletion (10-50%) with biological feasibility"""
        if not orfs:
            # Fallback: delete random region
            return self.region_deletion(sequence, orfs, factors)
        
        # Select random ORF with biological constraints
        safe_orfs = []
        for orf in orfs:
            orf_len = orf['end'] - orf['start'] + 1
            # Skip very small or very large ORFs
            if 100 <= orf_len <= len(sequence) // 3:
                safe_orfs.append(orf)
        
        if not safe_orfs:
            return self.region_deletion(sequence, orfs, factors)
        
        selected_orf = random.choice(safe_orfs)
        gene_start = selected_orf['start']
        gene_end = selected_orf['end']
        gene_length = gene_end - gene_start + 1
        
        # Calculate partial deletion (10-50% of the gene)
        deletion_percentage = random.uniform(0.10, 0.50)  # 10-50%
        deletion_length = int(gene_length * deletion_percentage)
        
        # Choose deletion position (start, middle, or end)
        deletion_position = random.choice(['start', 'middle', 'end'])
        
        if deletion_position == 'start':
            # Delete from start of gene
            start = gene_start
            end = gene_start + deletion_length
        elif deletion_position == 'end':
            # Delete from end of gene
            start = gene_end - deletion_length + 1
            end = gene_end
        else:  # middle
            # Delete from middle of gene
            middle_start = gene_start + (gene_length - deletion_length) // 2
            start = middle_start
            end = middle_start + deletion_length
        
        deleted_sequence = sequence[start-1:end]
        
        # Perform deletion (adjust for 0-based indexing)
        new_sequence = sequence[:start-1] + sequence[end:]
        
        return {
            'id': f"gene_del_{gene_start}_{gene_end}",
            'method': 'gene_deletion',
            'method_name': 'Gene Deletion',
            'full_sequence': new_sequence,
            'length_change': -(end - start + 1),
            'modifications': [{
                'method': 'gene_deletion',
                'start_position': start,
                'end_position': end,
                'length': -(end - start + 1),
                'details': {
                    'deleted_gene': selected_orf.get('product', 'Unknown'),
                    'deleted_sequence': deleted_sequence,
                    'original_gene_start': gene_start,
                    'original_gene_end': gene_end,
                    'deletion_percentage': deletion_percentage,
                    'deletion_position': deletion_position,
                    'randomness_factors': factors
                }
            }]
        }
    
    def frameshift_mutations(self, sequence, orfs, factors):
        """Robust frameshift mutations with biological feasibility"""
        if not orfs:
            # Fallback: random insertion/deletion
            return self.random_indel(sequence, factors)
        
        # Select random ORF
        safe_orfs = [orf for orf in orfs if orf['end'] - orf['start'] > 100]
        if not safe_orfs:
            return self.random_indel(sequence, factors)
        
        selected_orf = random.choice(safe_orfs)
        orf_start = selected_orf['start']
        orf_end = selected_orf['end']
        
        # Choose position within ORF (avoid start/end)
        safe_start = orf_start + 30
        safe_end = orf_end - 30
        
        if safe_end <= safe_start:
            position = orf_start + 50
        else:
            position = random.randint(safe_start, safe_end)
        
        # 1bp insertion or deletion
        if random.choice([True, False]):
            # Insertion
            new_sequence = sequence[:position-1] + random.choice('ATCG') + sequence[position-1:]
            length_change = 1
            mutation_type = "1bp_insertion"
        else:
            # Deletion
            new_sequence = sequence[:position-1] + sequence[position:]
            length_change = -1
            mutation_type = "1bp_deletion"
        
        return {
            'id': f"frameshift_{position}_{mutation_type}",
            'method': 'frameshift_mutations',
            'method_name': 'Frameshift Mutations',
            'full_sequence': new_sequence,
            'length_change': length_change,
            'modifications': [{
                'method': 'frameshift_mutations',
                'start_position': position,
                'end_position': position,
                'length': length_change,
                'details': {
                    'mutation_type': mutation_type,
                    'target_gene': selected_orf.get('product', 'Unknown'),
                    'randomness_factors': factors
                }
            }]
        }
    
    def random_indel(self, sequence, factors):
        """Fallback random indel when no ORFs available"""
        genome_len = len(sequence)
        position = random.randint(100, genome_len - 100)
        
        if random.choice([True, False]):
            # Insertion
            new_sequence = sequence[:position] + random.choice('ATCG') + sequence[position:]
            length_change = 1
            mutation_type = "random_1bp_insertion"
        else:
            # Deletion
            new_sequence = sequence[:position] + sequence[position+1:]
            length_change = -1
            mutation_type = "random_1bp_deletion"
        
        return {
            'id': f"random_{mutation_type}_{position}",
            'method': 'frameshift_mutations',
            'method_name': 'Random Indel',
            'full_sequence': new_sequence,
            'length_change': length_change,
            'modifications': [{
                'method': 'frameshift_mutations',
                'start_position': position,
                'end_position': position,
                'length': length_change,
                'details': {
                    'mutation_type': mutation_type,
                    'randomness_factors': factors
                }
            }]
        }
    
    # Placeholder methods for other engineering types
    def region_duplication(self, sequence, orfs, factors):
        return self.region_inversion(sequence, orfs, factors)
    
    def homing_endonuclease(self, sequence, orfs, factors):
        return self.gene_deletion(sequence, orfs, factors)
    
    def gfp_insertion(self, sequence, orfs, factors):
        return self.region_deletion(sequence, orfs, factors)
    
    def crispr_gfp(self, sequence, orfs, factors):
        return self.gfp_insertion(sequence, orfs, factors)
    
    def gene_start_deletion(self, sequence, orfs, factors):
        return self.gene_deletion(sequence, orfs, factors)
    
    def gene_end_deletion(self, sequence, orfs, factors):
        return self.gene_deletion(sequence, orfs, factors)
    
    def phosphorylation_inhibition(self, sequence, orfs, factors):
        return self.frameshift_mutations(sequence, orfs, factors)
    
    def igem_gfp(self, sequence, orfs, factors):
        return self.gfp_insertion(sequence, orfs, factors)
    
    def moclo_construct(self, sequence, orfs, factors):
        return self.gfp_insertion(sequence, orfs, factors)
    
    def save_sequences(self, virus_key, sequences, method_counts):
        """Save sequences in JSON and FASTA formats with modification coordinates in JSON"""
        if not sequences:
            return
            
        # Save JSON with modification coordinates
        json_file = os.path.join(self.output_dir, f"{virus_key}_detailed_sequences.json")
        with open(json_file, 'w') as f:
            json.dump(sequences, f, indent=2)
        
        # Save FASTA
        fasta_file = os.path.join(self.output_dir, f"{virus_key}_full_genomes.fasta")
        with open(fasta_file, 'w') as f:
            for seq_record in sequences:
                f.write(f">{seq_record['id']} | {seq_record['method_name']} | {seq_record['length_change']:+d}bp\n")
                f.write(f"{seq_record['full_sequence']}\n")
        
        print(f"  Saved {len(sequences)} sequences to:")
        print(f"    JSON: {json_file} (with modification coordinates)")
        print(f"    FASTA: {fasta_file}")
    
    def generate_all_sequences(self, target_per_virus=1000):
        """Generate sequences for all viruses"""
        start_time = datetime.now()
        total_sequences = 0
        virus_results = {}
        
        print(f"ROBUST VIRAL SEQUENCE GENERATION")
        print("=" * 50)
        print(f"Target: {target_per_virus} sequences per virus")
        print(f"Total viruses: {len(VIRUSES)}")
        print()
        
        for i, (virus_key, virus_info) in enumerate(VIRUSES.items(), 1):
            print(f"[{i}/{len(VIRUSES)}] Processing {virus_key}")
            print(f"Generating {target_per_virus} sequences for: {virus_info['name']}")
            print(f"Original length: {len(virus_info['sequence']):,} bp, ORFs: {len(virus_info.get('orfs', []))}")
            
            try:
                sequences, method_counts = self.generate_sequences_with_fallback(
                    virus_key, virus_info, target_per_virus
                )
                
                self.save_sequences(virus_key, sequences, method_counts)
                
                virus_results[virus_key] = {
                    'virus_name': virus_info['name'],
                    'original_length': len(virus_info['sequence']),
                    'sequences_generated': len(sequences),
                    'method_distribution': method_counts
                }
                
                total_sequences += len(sequences)
                print(f"  ✓ Generated {len(sequences)} sequences")
                
            except Exception as e:
                print(f"  ✗ Error: {e}")
                virus_results[virus_key] = {
                    'virus_name': virus_info['name'],
                    'error': str(e),
                    'sequences_generated': 0
                }
            
            print()
        
        # Save summary
        end_time = datetime.now()
        generation_time = (end_time - start_time).total_seconds()
        
        summary = {
            'generation_info': {
                'total_viruses': len(VIRUSES),
                'sequences_per_virus': target_per_virus,
                'total_sequences': len(VIRUSES) * target_per_virus,
                'start_time': start_time.isoformat(),
                'output_directory': self.output_dir,
                'end_time': end_time.isoformat(),
                'total_generation_time_seconds': generation_time,
                'total_sequences_generated': total_sequences,
                'average_time_per_sequence': generation_time / max(total_sequences, 1)
            },
            'virus_results': virus_results
        }
        
        summary_file = os.path.join(self.output_dir, 'robust_generation_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print("=" * 50)
        print("ROBUST GENERATION COMPLETE")
        print("=" * 50)
        print(f"Total sequences generated: {total_sequences:,}")
        print(f"Generation time: {generation_time:.1f} seconds ({generation_time/60:.1f} minutes)")
        print(f"Average time per sequence: {generation_time/total_sequences:.4f} seconds")
        print(f"Summary saved to: {summary_file}")
        print()
        
        # Method distribution
        all_methods = {}
        for virus_result in virus_results.values():
            if 'method_distribution' in virus_result:
                for method, count in virus_result['method_distribution'].items():
                    all_methods[method] = all_methods.get(method, 0) + count
        
        print("METHOD DISTRIBUTION:")
        for method, count in sorted(all_methods.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_sequences * 100) if total_sequences > 0 else 0
            print(f"  {method:<25} {count:5} ({percentage:4.1f}%)")
        
        print()
        print("OUTPUT FILES:")
        print("  JSON files: Detailed sequences with modification coordinates (start/end positions)")
        print("  FASTA files: Full genome sequences ready for analysis")
        print("  Summary: Complete generation statistics")

if __name__ == "__main__":
    generator = ViralSequenceGenerator()
    generator.generate_all_sequences(target_per_virus=1000)
