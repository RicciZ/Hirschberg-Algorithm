import argparse
import numpy as np

def gen_delta(c,d):
    if c == d:
        if c == '-':
            return -float('inf')
        else:
            return args.m
    elif c == '-':
        return args.i
    elif d == '-':
        return args.d
    else:
        return args.s

def traceback_global(v, w, pointers):
    i,j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    return ''.join(new_v[::-1])+'\n'+''.join(new_w[::-1])
    

def global_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as 
    computed by traceback_global. 
    
    :param: v
    :param: w
    :param: delta
    """
    print("The length of the input sequences are: n = ",len(v),", m = ",len(w))

    M = np.zeros((len(v)+1,len(w)+1))
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
    score, alignment = None, None
    for i in range(len(v)+1):
      for j in range(len(w)+1):
        if i == 0 and j == 0:
          continue
        if i == 0 and j > 0:
          M[i][j] = M[i][j-1] + delta['-'][w[j-1]]
          pointers[i][j] = LEFT
          continue
        if j == 0 and i > 0:
          M[i][j] = M[i-1][j] + delta[v[i-1]]['-']
          pointers[i][j] = UP
          continue
        up = M[i-1][j] + delta[v[i-1]]['-']
        left = M[i][j-1] + delta['-'][w[j-1]]
        topleft = M[i-1][j-1] + delta[v[i-1]][w[j-1]]
        if left >= up and left >= topleft:
          M[i][j] = left
          pointers[i][j] = LEFT
        elif up >= left and up >= topleft:
          M[i][j] = up
          pointers[i][j] = UP
        else:
          M[i][j] = topleft
          pointers[i][j] = TOPLEFT
    score = M[i][j]
    alignment = traceback_global(v,w, pointers)
    return score, alignment

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Hirschberg')
    parser.add_argument('--i', type=float, default=-1, help='Insertion score')
    parser.add_argument('--d', type=float, default=-1, help='Deletion score')
    parser.add_argument('--s', type=float, default=-1, help='Substitution score')
    parser.add_argument('--m', type=float, default=1, help='Match score')
    parser.add_argument('--seqA', type=str, default=None, help='The first input sequence')
    parser.add_argument('--seqB', type=str, default=None, help='The second input sequence')
    parser.add_argument('--lenA', type=int, default=100, help='The length of the first input sequence')
    parser.add_argument('--lenB', type=int, default=100, help='The length of the second input sequence')
    args = parser.parse_args()

    if args.seqA == None and args.seqB == None:
        seqA = np.random.choice(['A', 'C', 'T', 'G'], args.lenA)
        seqB = np.random.choice(['A', 'C', 'T', 'G'], args.lenB)
    elif args.seqA == None:
        seqA = np.random.choice(['A', 'C', 'T', 'G'], args.lenA)
        seqB = args.seqB.upper()
        for i in seqB:
            if i not in ['A', 'C', 'T', 'G']:
                raise ValueError('Invalid characters in the input sequence.')
    elif args.seqB == None:
        seqA = args.seqA.upper()
        seqB = np.random.choice(['A', 'C', 'T', 'G'], args.lenB)
        for i in seqA:
            if i not in ['A', 'C', 'T', 'G']:
                raise ValueError('Invalid characters in the input sequence.')
    else:
        seqA = args.seqA.upper()
        seqB = args.seqB.upper()
        for i in seqA:
            if i not in ['A', 'C', 'T', 'G']:
                raise ValueError('Invalid characters in the input sequence.')
        for i in seqB:
            if i not in ['A', 'C', 'T', 'G']:
                raise ValueError('Invalid characters in the input sequence.')
    
    keys = ['A', 'C', 'T', 'G', '-']
    delta = {}
    for i in range(len(keys)):
        delta[keys[i]] = {k: v for (k, v) in zip(
            keys, [gen_delta(keys[i],keys[j]) for j in range(len(keys))])}

    # pointers
    UP = (-1, 0)
    LEFT = (0, -1)
    TOPLEFT = (-1, -1)
    ORIGIN = (0, 0)

    score, align = global_align(seqA, seqB, delta)
    print(score)
    print(align)
    print()
    

