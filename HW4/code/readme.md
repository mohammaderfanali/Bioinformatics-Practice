# Bioinformatics Algorithms: Evolutionary Trees & Sequence Analysis 🧬

This repository contains essential bioinformatics algorithms implemented in Python and C++. These tools focus on reconstructing evolutionary relationships and identifying conserved motifs in genetic data.

---

## 📂 Projects Overview

| Project | Language | Description |
| :--- | :--- | :--- |
| **UPGMA Reconstructor** | Python | Reconstructs phylogenetic trees and calculates MRCA times from distance matrices. |
| **Greedy Motif Search** | C++ | Finds conserved regulatory patterns (motifs) in DNA using Laplace Pseudocounts. |

---

## 1. UPGMA Phylogenetic Reconstructor (Python) 🌳

This tool implements the **Unweighted Pair Group Method with Arithmetic Mean**. It is designed to take an evolutionary distance matrix and rebuild the most likely ultrametric tree.



### Key Features:
- **MRCA Identification**: Calculates the exact "Most Recent Common Ancestor" (MRCA) height for any pair of species.
- **Matrix Updating**: Uses weighted averages to maintain the ultrametric property of the tree throughout the clustering process.
- **Precision**: High-accuracy float processing with 6-decimal place output for scientific consistency.


## 2. Greedy Motif Search (C++) 🔍

A robust implementation of **motif discovery** in DNA sequences. This algorithm is essential for identifying **Transcription Factor Binding Sites (TFBS)**. It searches for hidden patterns across multiple DNA strings by iteratively optimizing a **profile matrix**.

---

### ✨ Key Features

- **Laplace Pseudocounts**  
  Uses the *Rule of Succession* (initializing the Profile Matrix with 1s) to handle the **zero-probability problem**, ensuring that rare nucleotides do not eliminate potential motifs.

- **Greedy Optimization**  
  Efficiently explores the search space by building the motif set one sequence at a time, making it significantly faster than exhaustive search methods.

- **Consensus Scoring**  
  Evaluates motifs based on the frequency of the most common nucleotide at each position, producing a score that reflects the strength of the motif.

---

### ⚙️ Technical Details

- **Profile Matrix**  
  A probability matrix \( P \) where each element \( P_{i,j} \) represents the likelihood of nucleotide \( i \) appearing at position \( j \).

- **Computational Complexity**  
  Optimized for high performance in **C++**, enabling efficient processing of large datasets where Python implementations may become slow.
