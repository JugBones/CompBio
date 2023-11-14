from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt

def get_user_sequences():
    sequences = {}
    while True:
        species = input("Enter species name (or 'done' to finish): ")
        if species.lower() == 'done':
            break
        sequence = input(f"Enter DNA sequence for {species}: ")
        sequences[species] = sequence
    return sequences

def build_phylogenetic_tree(sequences, method="identity", output_file="phylogenetic_tree.png"):
    sequence_records = [SeqRecord.SeqRecord(Seq(seq), id=species) for species, seq in sequences.items()]
    alignment = MultipleSeqAlignment(sequence_records)
    
    calculator = DistanceCalculator(method)
    dm = calculator.get_distance(alignment)
    
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    
    Phylo.draw(tree)
    plt.title("Phylogenetic Tree")
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    print("Phylogenetic Tree Builder")
    user_sequences = get_user_sequences()
    if not user_sequences:
        print("No sequences provided. Exiting.")
    else:
        method = input("Choose a distance calculation method (e.g., identity, kimura, etc.): ")
        output_file = input("Enter the output file name (e.g., phylogenetic_tree.png): ")
        build_phylogenetic_tree(user_sequences, method, output_file)
