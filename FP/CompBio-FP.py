import os
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from dendropy import TreeList, DnaCharacterMatrix, treecalc
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, faces, AttrFace
import time
import sys

# Step 1: Data Collection (example doang, ganti dataset sendiri nanti)
input_file = "input_sequences.fasta"
output_file = "aligned_sequences.fasta"

# Sequence Alignment
mafft_cmd = f"mafft --auto {input_file} > {output_file}"
os.system(mafft_cmd)

# Phylogenetic Tree Construction
alignment = AlignIO.read(output_file, "fasta")

# Maximum Likelihood Tree Construction
calculator_ml = DistanceCalculator("identity")
constructor_ml = DistanceTreeConstructor(calculator_ml, method="upgma")
ml_tree = constructor_ml.build_tree(alignment)

# Bayesian Inference Tree Construction 
char_matrix = DnaCharacterMatrix.get(file=open(output_file), schema="fasta")
tree_list = TreeList.get_from_path(output_file, schema="fasta")

burnin = 1000
post_burnin_trees = tree_list[burnin:]
bi_consensus_tree = treecalc.dendropy_simple_majority_consensus_tree(post_burnin_trees)

# measurement (metrics) : Calculate RF Distance
rf_distance = treecalc.symmetric_difference(ml_tree, bi_consensus_tree)
print(f"Robinson-Foulds Distance: {rf_distance}")

# measurement (metrics) : Time & Space Complexity
start_time = time.time()
ml_tree = constructor_ml.build_tree(alignment)
end_time = time.time()
time_ml = end_time - start_time
space_ml = sys.getsizeof(ml_tree)

start_time = time.time()
bi_consensus_tree = treecalc.dendropy_simple_majority_consensus_tree(post_burnin_trees)
end_time = time.time()
time_bi = end_time - start_time
space_bi = sys.getsizeof(bi_consensus_tree)

print(f"ML Tree - Time: {time_ml}s, Space: {space_ml} bytes")
print(f"BI Tree - Time: {time_bi}s, Space: {space_bi} bytes")

# Clustering and Visualization
ml_distance_matrix = treecalc.symmetric_difference_matrix(ml_tree, tree_list)
bi_distance_matrix = treecalc.symmetric_difference_matrix(bi_consensus_tree, tree_list)

# Perform hierarchical clustering
clustering_ml = AgglomerativeClustering(n_clusters=3).fit(ml_distance_matrix)
clustering_bi = AgglomerativeClustering(n_clusters=3).fit(bi_distance_matrix)

# Visualize clustering results
sns.clustermap(ml_distance_matrix, row_cluster=False, col_cluster=False, cmap="viridis", method="average")
plt.title("ML Tree Clustering")
plt.show()

sns.clustermap(bi_distance_matrix, row_cluster=False, col_cluster=False, cmap="viridis", method="average")
plt.title("BI Tree Clustering")
plt.show()

# Visualization with ETE Toolkit
ml_tree_ete = Tree(ml_tree.format("newick"))
ml_tree_ete.show(tree_style=TreeStyle())

bi_consensus_tree_ete = Tree(bi_consensus_tree.as_string(schema="newick"))
bi_consensus_tree_ete.show(tree_style=TreeStyle())
