from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _Matrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import ParsimonyScorer
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor
import matplotlib.pyplot as plt
import os
import time
import tracemalloc
import tkinter as tk
from tkinter import filedialog


class PhylogeneticApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Phylogenetic Analysis App")

        # Variables
        self.alignment_file = ""
        self.alignment_label = tk.StringVar()
        self.ml_tree = None
        self.file_labels = [None]

        # Load Alignment Button and Label
        align_frame = tk.Frame(root)
        align_frame.grid(row=0, column=0, padx=10, pady=5, sticky="w")

        align_btn = tk.Button(
            align_frame, text="Load Alignment File", command=self.load_alignment_file
        )
        align_btn.grid(row=0, column=0)

        align_label = tk.Label(
            align_frame, textvariable=self.alignment_label, wraplength=200
        )
        align_label.grid(row=0, column=1, padx=10, pady=5, sticky="w")

        # Phylogenetic Analysis Button
        nj_btn = tk.Button(
            root,
            text="NJ",
            command=self.nj_analysis,
        )
        nj_btn.grid(row=1, column=0, pady=10, columnspan=2)

        upgma_btn = tk.Button(
            root,
            text="UPGMA",
            command=self.upgma_analysis,
        )
        upgma_btn.grid(row=1, column=1, pady=10, columnspan=2)

        # Adjust window size
        self.root.geometry("400x200")  # Set the width and height of the window

    def update_labels(self, index=None):
        if index is not None and 0 <= index < 5:
            # Update file label
            self.file_labels[index].config(text=self.sequence_files[index])

        # Update alignment label
        align_label = self.root.nametowidget(".!frame.!label")
        align_label.config(text=self.alignment_file)

    def load_alignment_file(self):
        file_path = filedialog.askopenfilename(
            title="Select Alignment File", filetypes=[("Alignment Files", "*.aln")]
        )
        if file_path:
            self.alignment_file = file_path
            print(f"Alignment file loaded: {file_path}")

            # Update the label
            self.alignment_label.set(file_path)  # Set the text variable

    def nj_analysis(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fasta_folder = os.path.join(script_dir, "falcons-fasta")
        # Import sequences
        t1 = SeqIO.read(os.path.join(fasta_folder, "falco-alopex.fasta"), "fasta")
        t2 = SeqIO.read(os.path.join(fasta_folder, "falco-amurensis.fasta"), "fasta")
        t3 = SeqIO.read(os.path.join(fasta_folder, "falco-biarmicus.fasta"), "fasta")
        t4 = SeqIO.read(os.path.join(fasta_folder, "falco-cherrug.fasta"), "fasta")
        t5 = SeqIO.read(os.path.join(fasta_folder, "falco-chicquera.fasta"), "fasta")
        t6 = SeqIO.read(os.path.join(fasta_folder, "falco-columbarius.fasta"), "fasta")
        t7 = SeqIO.read(os.path.join(fasta_folder, "falco-cuvierii.fasta"), "fasta")
        t8 = SeqIO.read(os.path.join(fasta_folder, "falco-dickinsoni.fasta"), "fasta")
        t9 = SeqIO.read(os.path.join(fasta_folder, "falco-eleonorae.fasta"), "fasta")
        t10 = SeqIO.read(os.path.join(fasta_folder, "falco-fasciinucha.fasta"), "fasta")
        t11 = SeqIO.read(os.path.join(fasta_folder, "falco-jugger.fasta"), "fasta")
        t12 = SeqIO.read(os.path.join(fasta_folder, "falco-longipennis.fasta"), "fasta")
        t13 = SeqIO.read(os.path.join(fasta_folder, "falco-mexicanus.fasta"), "fasta")
        t14 = SeqIO.read(os.path.join(fasta_folder, "falco-zoniventris.fasta"), "fasta")
        t15 = SeqIO.read(os.path.join(fasta_folder, "sooty-falcon.fasta"), "fasta")

        # Rename sequences
        sequences = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
        sequence_ids = [
            "Kestrel",
            "Amur",
            "Lanner",
            "Saker",
            "RedNeck",
            "Merlin",
            "AfrHobby",
            "Dickinson",
            "Eleonora",
            "Taita",
            "Laggar",
            "AusHobby",
            "Prairie",
            "Banded",
            "Sooty",
        ]

        # Combine sequences into a new file
        combined_file = SeqIO.write(sequences, "typicalfalcons.fasta", "fasta")

        # Open alignment file
        with open(self.alignment_file, "r") as aln:
            self.alignment = AlignIO.read(aln, "clustal")

        tracemalloc.start()
        start_time = time.time()

        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(self.alignment)
        print(distance_matrix)

        # NJ Tree Construction using Bio.Phylo with labels for internal nodes
        nj_constructor = DistanceTreeConstructor(calculator, "nj")
        self.ml_tree = nj_constructor.build_tree(self.alignment)
        self.ml_tree.rooted = True
        end_time = time.time()

        memory_usage = tracemalloc.get_traced_memory()

        tracemalloc.stop()
        # Create a new Toplevel window for metrics
        metrics_popup = tk.Toplevel(self.root)
        metrics_popup.title("Phylogenetic NJ Analysis Metrics")

        metrics_label = tk.Label(
            metrics_popup,
            text=f"Time taken: {end_time - start_time} seconds\nMemory usage: {memory_usage}",
        )
        metrics_label.pack()
        print(self.ml_tree)

        Phylo.write(self.ml_tree, "nj_falcons_tree.xml", "phyloxml")

        def get_species_name(clade):
            # For leaf nodes, return the species name
            if clade.is_terminal():
                return clade.name.split("_")[
                    0
                ]  # Assuming the species name is the part before the underscore
            else:
                return ""  # For internal nodes, return an empty string

        # Visualize ML tree
        fig, ax = plt.subplots(figsize=(18, 8), dpi=100)
        Phylo.draw(
            self.ml_tree,
            axes=ax,
            do_show=False,
            label_func=get_species_name,
            branch_labels=None,
        )

        # Add labels and customize layout
        ax.set_title("NJ Phylogenetic Tree")
        ax.set_xlabel("Evolutionary Distance")
        ax.set_ylabel("Species")

        # Save the figure
        fig.savefig("nj_falcons_phylogenetic_tree.png")

        # Convert the tree to a different format (optional)
        Phylo.convert("nj_falcons_tree.xml", "phyloxml", "nj-falcons_tree.nex", "nexus")

        falcon_nex = Phylo.read("nj-falcons_tree.nex", "nexus")

        fig = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
        plt.rc("font", size=12)  # fontsize of the leaf and node labels
        plt.rc("xtick", labelsize=10)  # fontsize of the tick labels
        plt.rc("ytick", labelsize=10)  # fontsize of the tick labels
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(falcon_nex, axes=axes, branch_labels=None)
        fig.savefig("typicalfalcons2_cladogram")

    def upgma_analysis(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fasta_folder = os.path.join(script_dir, "falcons-fasta")
        # Import sequences
        t1 = SeqIO.read(os.path.join(fasta_folder, "falco-alopex.fasta"), "fasta")
        t2 = SeqIO.read(os.path.join(fasta_folder, "falco-amurensis.fasta"), "fasta")
        t3 = SeqIO.read(os.path.join(fasta_folder, "falco-biarmicus.fasta"), "fasta")
        t4 = SeqIO.read(os.path.join(fasta_folder, "falco-cherrug.fasta"), "fasta")
        t5 = SeqIO.read(os.path.join(fasta_folder, "falco-chicquera.fasta"), "fasta")
        t6 = SeqIO.read(os.path.join(fasta_folder, "falco-columbarius.fasta"), "fasta")
        t7 = SeqIO.read(os.path.join(fasta_folder, "falco-cuvierii.fasta"), "fasta")
        t8 = SeqIO.read(os.path.join(fasta_folder, "falco-dickinsoni.fasta"), "fasta")
        t9 = SeqIO.read(os.path.join(fasta_folder, "falco-eleonorae.fasta"), "fasta")
        t10 = SeqIO.read(os.path.join(fasta_folder, "falco-fasciinucha.fasta"), "fasta")
        t11 = SeqIO.read(os.path.join(fasta_folder, "falco-jugger.fasta"), "fasta")
        t12 = SeqIO.read(os.path.join(fasta_folder, "falco-longipennis.fasta"), "fasta")
        t13 = SeqIO.read(os.path.join(fasta_folder, "falco-mexicanus.fasta"), "fasta")
        t14 = SeqIO.read(os.path.join(fasta_folder, "falco-zoniventris.fasta"), "fasta")
        t15 = SeqIO.read(os.path.join(fasta_folder, "sooty-falcon.fasta"), "fasta")

        # Rename sequences
        sequences = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
        sequence_ids = [
            "Kestrel",
            "Amur",
            "Lanner",
            "Saker",
            "RedNeck",
            "Merlin",
            "AfrHobby",
            "Dickinson",
            "Eleonora",
            "Taita",
            "Laggar",
            "AusHobby",
            "Prairie",
            "Banded",
            "Sooty",
        ]

        # Combine sequences into a new file
        combined_file = SeqIO.write(sequences, "typicalfalcons.fasta", "fasta")

        # Open alignment file
        with open(self.alignment_file, "r") as aln:
            self.alignment = AlignIO.read(aln, "clustal")

        tracemalloc.start()
        start_time = time.time()

        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(self.alignment)
        print(distance_matrix)

        # NJ Tree Construction using Bio.Phylo with labels for internal nodes
        upgma_constructor = DistanceTreeConstructor(calculator, "upgma")
        self.ml_tree = upgma_constructor.build_tree(self.alignment)
        self.ml_tree.rooted = True

        end_time = time.time()

        memory_usage = tracemalloc.get_traced_memory()

        tracemalloc.stop()

        # Create a new Toplevel window for metrics
        metrics_popup = tk.Toplevel(self.root)
        metrics_popup.title("Phylogenetic UPGMA Analysis Metrics")

        metrics_label = tk.Label(
            metrics_popup,
            text=f"Time taken: {end_time - start_time} seconds\nMemory usage: {memory_usage}",
        )
        metrics_label.pack()

        print(self.ml_tree)

        Phylo.write(self.ml_tree, "upgma_falcons_tree.xml", "phyloxml")

        # Visualize ML tree
        def get_species_name(clade):
            # For leaf nodes, return the species name
            if clade.is_terminal():
                return clade.name.split("_")[
                    0
                ]  # Assuming the species name is the part before the underscore
            else:
                return ""  # For internal nodes, return an empty string

        # Visualize ML tree
        fig, ax = plt.subplots(figsize=(18, 8), dpi=100)
        Phylo.draw(
            self.ml_tree,
            axes=ax,
            do_show=False,
            label_func=get_species_name,
            branch_labels=None,
        )

        # Add labels and customize layout
        ax.set_title("UPGMA Phylogenetic Tree")
        ax.set_xlabel("Evolutionary Distance")
        ax.set_ylabel("Species")

        # Save the figure
        fig.savefig("upgma_falcons_phylogenetic_tree.png")

        # Convert the tree to a different format (optional)
        Phylo.convert(
            "upgma_falcons_tree.xml", "phyloxml", "upgma-falcons_tree.nex", "nexus"
        )

        falcon_nex = Phylo.read("upgma-falcons_tree.nex", "nexus")

        fig = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
        plt.rc("font", size=12)  # fontsize of the leaf and node labels
        plt.rc("xtick", labelsize=10)  # fontsize of the tick labels
        plt.rc("ytick", labelsize=10)  # fontsize of the tick labels
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(falcon_nex, axes=axes)
        fig.savefig("typicalfalcons2_cladogram")


if __name__ == "__main__":
    root = tk.Tk()
    app = PhylogeneticApp(root)
    root.mainloop()
