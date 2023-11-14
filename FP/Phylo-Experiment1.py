from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt

# Sample DNA sequences (replace with your own data)
sequences = {
    "Species_A": "ATCGATCGATCG",
    "Species_B": "ATAGATCGATCG",
    "Species_C": "ATCGATCGCTCG",
    "Species_D": "ATGGATCGATCG",
}

# Create a list of sequence records
sequence_records = [SeqRecord.SeqRecord(Seq(seq), id=species) for species, seq in sequences.items()]

# Create a multiple sequence alignment object
alignment = MultipleSeqAlignment(sequence_records)

# Calculate the distance matrix
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

# Build the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

# Draw and save the phylogenetic tree as a PNG file
Phylo.draw(tree)
plt.title("Phylogenetic Tree")
plt.savefig("phylogenetic_tree.png")
plt.show()
