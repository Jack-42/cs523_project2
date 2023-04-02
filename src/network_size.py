import math
from decimal import Decimal


# note: r_s is .8523 for the RBD

def synonymous_size(r_s, length, subs, k):
    """
    calculate the size of a k-step neutral network
    params:
    r_s -> the ratio of synonymous mutations
    length -> length of the sequence to evaluate
    subs -> the number of possible substitutions (e.g. 3 for nucleotide sequences, 19 for amino acids)
    k -> k steps
    """

    sum = 0
    for step in range(1, k + 1):
        N_k = math.pow(r_s, step) * math.pow(subs,step) * (length - step + 1)
        sum += N_k

    return sum


def kplus1(r_s, length, subs, k, k_max=0):
    """
    calcluate the number of sequences 1 step away from a k-step
    neutral network.
    params:
    r_s -> the ratio of synonymous mutations
    length -> length of the sequence to evaluate
    subs -> the number of possible substitutions (e.g. 3 for nucleotide sequences, 19 for amino acids)
    k -> start at k steps
    k_max -> evaluate k+k_max rows and output as csv
    """

    results_path = "../results/kplus1.csv"

    if k < 1:
        raise Exception('k must be at least 1')
    
    results = []

    for k_step in range(k, k + k_max):
        N_k = math.pow(r_s, k_step) * math.pow(subs,k_step) * (length - k_step + 1)
        N_syn = r_s * N_k
        N_syn_str = str(N_syn)
        if (N_syn >= 10000000):
            N_syn_str = '%.2E' % Decimal(N_syn)
        N_nonsyn = N_k - N_syn
        N_nonsyn_str = str(N_nonsyn)
        if (N_nonsyn >= 10000000):
            N_nonsyn_str = '%.2E' % Decimal(N_nonsyn)
        
        results.append([k_step, N_syn_str, N_nonsyn_str])
    with open(results_path, 'w') as f:
        for result in results:
            f.write(str(result[0]) + ',' + str(result[1]) + ',' + str(result[2]))
            f.write("\n")

def qb4():
    r_s = .8523
    length = 200
    subs = 19
    k = 4

    syn = synonymous_size(r_s, length, subs, k)
    print('%.3E' % Decimal(syn))
    tot = (1/r_s) * syn
    print('%.3E' % Decimal(tot))

    return syn, tot

def qb5():
    syn, tot = qb4()
    diff = tot - syn
    print('%.3E' % Decimal(diff))

def qb7():
    r_s = .8523
    length = 200
    subs = 19
    k = 1
    k_max = 8

    kplus1(r_s, length, subs, k, k_max)


if __name__ == "__main__":
    #synonymous_size(.8523, 200, 19, 3)
    #qb4()
    #qb5()
    qb7()
