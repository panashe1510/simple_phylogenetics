# Import required modules from Biopython and Matplotlib
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import matplotlib.pyplot as plt

# Function to read sequences from a FASTA file


def read_sequences(fasta_file):
    return list(SeqIO.parse(fasta_file, "fasta"))

# Function to calculate distance between two sequences
# Uses global alignment to compute % identity, then returns (1 - identity) as a distance


def calculate_identity(seq1, seq2):
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    matches = sum(res1 == res2 for res1, res2 in zip(
        best_alignment.target, best_alignment.query))
    identity = matches / max(len(seq1), len(seq2))  # normalized identity
    return 1 - identity  # distance = 1 - identity

# Builds a distance matrix from all pairwise sequence comparisons


def build_distance_matrix(records):
    names = [record.id for record in records]  # Extract species names
    matrix = [[0] * (i + 1)
              for i in range(len(records))]  # Lower triangle format

    for i in range(len(records)):
        for j in range(i):
            dist = calculate_identity(str(records[i].seq), str(records[j].seq))
            matrix[i][j] = round(dist, 4)  # Fill only the lower triangle

    return DistanceMatrix(names, matrix)

# Function to show alignments between each sequence pair (optional output)


def show_alignments(records):
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            print(f"\nAlignment: {records[i].id} vs {records[j].id}")
            alignments = pairwise2.align.globalxx(
                records[i].seq, records[j].seq, one_alignment_only=True)
            print(pairwise2.format_alignment(*alignments[0]))

# Main function that ties everything together


def main():
    fasta_file = "mini_genes.fasta"
    records = read_sequences(fasta_file)

    print(" Building distance matrix...")
    matrix = build_distance_matrix(records)

    print(" Constructing phylogenetic tree (UPGMA)...")
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(matrix)

    print(" Visualizing phylogenetic tree...")
    Phylo.draw(tree)
    plt.show(block=False)
    Phylo.write(tree, "simple_tree.nwk", "newick")
    print("Tree saved as simple_tree.nwk")

    #  Bonus: Ask user if they want to inspect alignments
    user_input = input(
        "\n Do you want to see alignment scores for each pair? (y/n): ").lower()
    if user_input == 'y':
        show_alignments(records)


#  Run the script
if __name__ == "__main__":
    main()
