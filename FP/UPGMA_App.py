import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import time
import tracemalloc


class PhylogeneticApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Phylogenetic Analysis App")

        # Variables
        self.sequence_files = [""] * 5
        self.alignment_file = ""
        self.alignment = None
        self.ml_tree = None
        self.file_labels = [None] * 5  # To store the labels dynamically

        # Load File Buttons and Labels
        for i in range(5):
            frame = tk.Frame(root)
            frame.grid(row=i, column=0, padx=10, pady=5, sticky="w")

            btn = tk.Button(
                frame, text=f"Load File {i+1}", command=lambda i=i: self.load_file(i)
            )
            btn.grid(row=0, column=0)

            # Create a label to the right of the button
            label = tk.Label(frame, text="", wraplength=200)
            label.grid(row=0, column=1, padx=10, pady=5, sticky="w")
            self.file_labels[i] = label

        # Load Alignment Button and Label
        align_frame = tk.Frame(root)
        align_frame.grid(row=5, column=0, padx=10, pady=5, sticky="w")

        align_btn = tk.Button(
            align_frame, text="Load Alignment File", command=self.load_alignment_file
        )
        align_btn.grid(row=0, column=0)

        align_label = tk.Label(align_frame, text="", wraplength=200)
        align_frame.grid(row=5, column=0, padx=10, pady=5, sticky="w")

        # Phylogenetic Analysis Button
        analyze_btn = tk.Button(
            root,
            text="Perform Phylogenetic Analysis",
            command=self.perform_phylogenetic_analysis,
        )
        analyze_btn.grid(row=6, column=0, pady=10, columnspan=2)

    def load_file(self, index):
        file_path = filedialog.askopenfilename(
            title=f"Select File {index+1}", filetypes=[("FASTA Files", "*.fasta")]
        )
        if file_path:
            self.sequence_files[index] = file_path
            print(f"File {index+1} loaded: {file_path}")

            # Update the label
            self.update_labels(index)

    def update_labels(self, index=None):
        if index is not None and 0 <= index < 5:
            # Update file label
            self.file_labels[index].config(text=self.sequence_files[index])

        # Update alignment label
        align_label = self.root.nametowidget(".!frame6.!label")
        align_label.config(text=self.alignment_file)

    def load_alignment_file(self):
        file_path = filedialog.askopenfilename(
            title="Select Alignment File", filetypes=[("Alignment Files", "*.aln")]
        )
        if file_path:
            self.alignment_file = file_path
            print(f"Alignment file loaded: {file_path}")

            # Update the label
            self.update_labels()

    def perform_phylogenetic_analysis(self):
        if (
            any(file_path == "" for file_path in self.sequence_files)
            or self.alignment_file == ""
        ):
            print("Please load all five sequence files and an alignment file.")
            return

        # Read and combine sequences
        sequences = [
            SeqIO.read(file_path, "fasta") for file_path in self.sequence_files
        ]
        combined_file = "combined.fasta"
        SeqIO.write(sequences, combined_file, "fasta")

        # Load the combined sequences into MUSCLE
        # ... (as in your original code)

        # Open the alignment file as a MultipleSeqAlignment object
        with open(self.alignment_file, "r") as aln:
            self.alignment = AlignIO.read(aln, "clustal")

        # Calculate Distance Matrix
        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(self.alignment)

        # ML Tree Construction using Bio.Phylo
        ml_constructor = DistanceTreeConstructor(calculator, "nj")
        self.ml_tree = ml_constructor.build_tree(self.alignment)

        # Save ML tree to a new file (optional)
        Phylo.write(self.ml_tree, "nj_waterfrog_tree.newick", "newick")

        # Display visualization
        self.visualize_tree()

        # Display metrics
        self.display_metrics()

    def visualize_tree(self):
        fig, ax = plt.subplots(figsize=(15, 8), dpi=100)
        Phylo.draw(self.ml_tree, axes=ax, do_show=False)

        # Add labels and customize layout
        ax.set_title("NJ Phylogenetic Tree")
        ax.set_xlabel("Evolutionary Distance")
        ax.set_ylabel("Species")

        # Save the figure
        fig.savefig("nj_waterfrog_phylogenetic_tree.png")
        plt.show()

    def display_metrics(self):
        tracemalloc.start()

        start_time = time.time()

        # Create a DistanceCalculator object
        calculator = DistanceCalculator("identity")

        # Calculate the distance matrix
        distance_matrix = calculator.get_distance(self.alignment)

        # Create a DistanceTreeConstructor object
        ml_constructor = DistanceTreeConstructor(calculator, "upgma")

        # Build the ML tree
        self.ml_tree = ml_constructor.build_tree(self.alignment)

        end_time = time.time()

        memory_usage = tracemalloc.get_traced_memory()

        tracemalloc.stop()

        # Create a new Toplevel window for metrics
        metrics_popup = tk.Toplevel(self.root)
        metrics_popup.title("Phylogenetic Analysis Metrics")

        metrics_label = tk.Label(
            metrics_popup,
            text=f"Time taken: {end_time - start_time} seconds\nMemory usage: {memory_usage}",
        )
        metrics_label.pack()

        # You can add more widgets or customize the layout in the metrics_popup window


if __name__ == "__main__":
    root = tk.Tk()
    app = PhylogeneticApp(root)
    root.mainloop()
