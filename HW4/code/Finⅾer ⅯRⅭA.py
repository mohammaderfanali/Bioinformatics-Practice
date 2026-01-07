import sys

def get_names_and_n():
    try:
        n = int(input().strip())
        names = input().split()
        return n, names
    except EOFError:
        return None, None

def get_distance_matrix(n):
    matrix = []
    for _ in range(n):
        row = list(map(float, input().split()))
        matrix.append(row)
    return matrix

def update_mrca_results(mrca_results, c1_members, c2_members, time):
    for m1 in c1_members:
        for m2 in c2_members:
            mrca_results[m1][m2] = mrca_results[m2][m1] = time

def run_upgma_process(n, dist_matrix):
    clusters = [[i] for i in range(n)]
    mrca_results = [[0.0] * n for _ in range(n)]
    active_indices = list(range(n))
    
    while len(active_indices) > 1:
        min_d = float('inf')
        u_idx, v_idx = -1, -1
        
        for i in range(len(active_indices)):
            for j in range(i + 1, len(active_indices)):
                idx1, idx2 = active_indices[i], active_indices[j]
                if dist_matrix[idx1][idx2] < min_d:
                    min_d = dist_matrix[idx1][idx2]
                    u_idx, v_idx = i, j
        
        c1 = active_indices[u_idx]
        c2 = active_indices[v_idx]
        mrca_time = min_d / 2.0
        
        update_mrca_results(mrca_results, clusters[c1], clusters[c2], mrca_time)
        
        size1, size2 = len(clusters[c1]), len(clusters[c2])
        total_size = size1 + size2
        
        for k in range(n):
            if k != c1 and k != c2:
                # D(C_new, k) = (size1 * D(c1, k) + size2 * D(c2, k)) / total_size
                new_dist = (size1 * dist_matrix[c1][k] + size2 * dist_matrix[c2][k]) / total_size
                dist_matrix[c1][k] = dist_matrix[k][c1] = new_dist
        
        clusters[c1].extend(clusters[c2])
        active_indices.pop(v_idx)
        
    return mrca_results

def process_queries(mrca_results, names):
    name_to_idx = {name: i for i, name in enumerate(names)}
    try:
        q_count = int(input().strip())
        for _ in range(q_count):
            u_name, v_name = input().split()
            idx1, idx2 = name_to_idx[u_name], name_to_idx[v_name]
            
            print(f"{mrca_results[idx1][idx2]:.6f}")
    except EOFError:
        pass

def main():
    n, names = get_names_and_n()
    if n is None: return
    
    dist_matrix = get_distance_matrix(n)
    
    mrca_results = run_upgma_process(n, dist_matrix)
    
    process_queries(mrca_results, names)

if __name__ == "__main__":
    main()