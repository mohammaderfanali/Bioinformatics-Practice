import sys
def get_score(motifs):
  
    k = len(motifs[0])
    t = len(motifs)
    score = 0
    for j in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for i in range(t):
            counts[motifs[i][j]] += 1
        score += max(counts.values())
    return score

def create_profile(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = {'A': [1.0] * k, 'C': [1.0] * k, 'G': [1.0] * k, 'T': [1.0] * k}
    
    for j in range(k):
        for i in range(t):
            profile[motifs[i][j]][j] += 1
        
        for char in 'ACGT':
            profile[char][j] /= (t + 4)
    return profile

def profile_most_probable(text, k, profile):
    max_prob = -1.0
    best_kmer = text[0:k]
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1.0
        for j in range(k):
            prob *= profile[kmer[j]][j]
        
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

def greedy_motif_search(dna, k, t):
    best_motifs = [seq[0:k] for seq in dna]
    best_score = get_score(best_motifs)
    
    first_seq = dna[0]
    for i in range(len(first_seq) - k + 1):
        motifs = [first_seq[i:i+k]]
        
        for j in range(1, t):
            current_profile = create_profile(motifs)
            next_motif = profile_most_probable(dna[j], k, current_profile)
            motifs.append(next_motif)
            
        current_score = get_score(motifs)
        if current_score > best_score:
            best_score = current_score
            best_motifs = motifs
            
    return best_motifs

def main():
    line = sys.stdin.readline().split()
    if not line:
        return
    
    k = int(line[0])
    t = int(line[1])
    
    dna_sequences = []
    for _ in range(t):
        seq = sys.stdin.readline().strip()
        if seq:
            dna_sequences.append(seq)
    
    result = greedy_motif_search(dna_sequences, k, t)
    
    for motif in result:
        print(motif)

if __name__ == "__main__":
    main()