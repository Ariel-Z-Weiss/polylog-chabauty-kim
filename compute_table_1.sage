# The authors thank Martin Lüdtke for for sharing a preliminary version of this code, which we modified slightly.

"""
    Given a tuple of weights (w_1, ..., w_k) and a bound N, return a list of length N+1
    whose n-th entry is the number of vectors (a_1, ..., a_k) with sum(a_i * w_i) = n
"""
def count_weighted_partitions(N, weights):
    count = [0] * (N + 1)
    
    # Base case: There is exactly one way to get a total of 0:
    # by taking 0 copies of each weight
    count[0] = 1
    
    for w in weights:
        for n in range(w, N + 1):
            # Given a vector that sums to (n-w), we can add 1 in w's position to form n. 
            count[n] += count[n-w]
    return count

max_depth = 30
max_weight = 3000

for d in range(2,max_depth+1):
    shuffle_dims = count_weighted_partitions(max_weight, (1,)*4 + tuple(range(3, d+1, 2)))
    polylog_dims = count_weighted_partitions(max_weight, (1,) + tuple(range(1, d+1)))

    # Find the first weight where the dim_PL(d, v) > dim_\Phi(d, v)
    for w in range(2, max_weight + 1):
        if shuffle_dims[w] < polylog_dims[w]:
            print(f"d={d}, w={w} ({shuffle_dims[w]} < {polylog_dims[w]})")
            break
