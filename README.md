#  Simple Phylogenetic Tree Constructor

 A Python-based tool for constructing **phylogenetic trees** from DNA sequences using pairwise alignment and UPGMA clustering. Ideal for bioinformatics students, researchers, and computational biology enthusiasts.

##  Features
- Reads DNA sequences from a FASTA file  
- Performs **pairwise sequence alignment** using Biopython  
- Builds a **distance matrix** based on percent identity  
- Constructs a **phylogenetic tree using UPGMA**  
-  Visualizes the tree using Matplotlib  
-  Saves the tree in **Newick format** for further analysis  
- Allows users to **view sequence alignment scores**  

---

##  Repository Contents

📁 simple-phylogenetics/ ├── simple_phylo.py      # Python script for tree construction ├── mini_genes.fasta          # Output phylogenetic tree (Newick format) ├── README.md


---

## 🔧 Installation & Setup

Ensure you have **Python 3.7+** installed, then install the required dependencies:

bash
pip install biopython matplotlib

---
##  Running the code
1️⃣ Ensure mini_genes.fasta exists in the project folder.
2️⃣ Run the script

bash
python simple_phylo.py

##  Example output

      ┌── Species_A
      │      
      │    ┌── Species_B
      │    │      
      │    │    ┌── Species_C
      │    │    │      
      └────└────└── Species_D
                  │      
                  └── Species_E


