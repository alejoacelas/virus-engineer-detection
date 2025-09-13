# Real virus sequences for engineering
def load_real_sequences():
    """Load real sequences from FASTA files"""
    def clean_sequence(content):
        """Extract sequence from FASTA content"""
        lines = content.strip().split('\n')
        # Skip header line (starts with >)
        sequence_lines = [line for line in lines if not line.startswith('>')]
        return ''.join(sequence_lines).replace('\n', '').replace(' ', '')
    
    try:
        with open("corona_real.txt", "r") as f:
            corona_content = f.read()
            corona_seq = clean_sequence(corona_content)
        
        with open("adeno_real.txt", "r") as f:
            adeno_content = f.read()
            adeno_seq = clean_sequence(adeno_content)
            
    except FileNotFoundError:
        print("Real sequence files not found. Please save corona_real.txt and adeno_real.txt")
        return {}
    
    return {
        "coronavirus_nl63": {
            "name": "Human Coronavirus NL63 (NC_005831.2)",
            "sequence": corona_seq,
            "length": len(corona_seq),
            "family": "coronavirus"
        },
        
        "adenovirus_54": {
            "name": "Human Adenovirus 54 (NC_012959.1)", 
            "sequence": adeno_seq,
            "length": len(adeno_seq),
            "family": "adenovirus"
        }
    }

VIRUSES = load_real_sequences()

# Common restriction sites for realistic engineering
RESTRICTION_SITES = {
    "BamHI": "GGATCC",
    "EcoRI": "GAATTC", 
    "HindIII": "AAGCTT",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "NotI": "GCGGCCGC",
    "SmaI": "CCCGGG"
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
    }
}
