def validate_dna(dna_seq):
    # Convert to uppercase
    seqm = dna_seq.upper()
    print(seqm)

    valid = seqm.count("A") + seqm.count("T") + seqm.count("C") + seqm.count("G")

    if valid == len(seqm):
        print("VALID!")
        return True
    else:
        print("INVALID!")
        return False

def calcfreq(seq):
    dictionary = {}

    for s in seq.upper():
        if s in dictionary:
            dictionary[s] += 1
        else:
            dictionary[s] = 1
    return dictionary

def complement_dna(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base.upper()] for base in seq)

def dna_to_rna(seq):
    rna_dict = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
    return ''.join(rna_dict[base.upper()] for base in seq)

def rna_to_aminoacid(mrna_seq):
    codon_table = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    amino_acids = []
    codon = ''

    for base in mrna_seq:
        codon += base
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'Unknown')
            amino_acids.append(amino_acid)
            codon = ''

    return amino_acids



# RESULT & TESTING 
dna = "ttacga"
print("Input DNA = ", end="")
validate_result = validate_dna(dna)

if validate_result == True:
    freq_data = calcfreq(dna)
    print(freq_data)

    comp_dna = complement_dna(dna)
    print("Complement = ", comp_dna)

    transcrib = dna_to_rna(comp_dna)
    print("mRNA = ", transcrib)

    amino_acids = rna_to_aminoacid(transcrib)
    print("Aminoacid = ", amino_acids)
