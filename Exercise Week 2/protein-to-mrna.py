amino_acid_codon_mapping = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    '*': ['UAA', 'UAG', 'UGA']  # End of codons
}

# Function to translate an amino acid sequence into mRNA
def translate_to_mRNA(amino_acid_sequence):
    mRNA_sequence = []
    codon_frequency = {}
    
    for amino_acid in amino_acid_sequence:
        if amino_acid == '-':
            continue  # Skip dashes
        codons = amino_acid_codon_mapping.get(amino_acid, [])
        if not codons:
            raise ValueError(f"No codons found for amino acid {amino_acid}")
        
        # Add the first codon to the mRNA sequence
        mRNA_sequence.append(codons[0])
        
        # Update the codon frequency
        codon_frequency[codons[0]] = codon_frequency.get(codons[0], 0) + 1
    
    return ''.join(mRNA_sequence), codon_frequency

# Testing and Result
amino_acid_sequence = "N-A-N"
amino_acid_sequence = amino_acid_sequence.replace('-', '')  # Remove dashes
mRNA, codon_frequency = translate_to_mRNA(amino_acid_sequence)

print(f"mRNA: {mRNA}")
print("Codon Frequency:")
for codon, frequency in codon_frequency.items():
    print(f"{codon}: {frequency}")
