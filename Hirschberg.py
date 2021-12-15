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

def last_col(v, w, delta):
    """
    Returns the scores of the last column of alignment of the strings v and w

    :param: v
    :param: w
    :param: delta
    """
    pre_col = np.zeros(len(v)+1, int)
    cur_col = np.zeros(len(v)+1, int)

    for i in range(1, len(v)+1):
        pre_col[i] = pre_col[i-1] + delta[v[i-1]]['-']

    for j in range(1, len(w)+1):
        cur_col[0] = pre_col[0] + delta['-'][w[j-1]]
        for i in range(1, len(v)+1):
            up = cur_col[i-1] + delta[v[i-1]]['-']
            left = pre_col[i] + delta['-'][w[j-1]]
            topleft = pre_col[i-1] + delta[v[i-1]][w[j-1]]
            cur_col[i] = max(up, left, topleft)
        # print(pre_col)
        # print(cur_col)
        pre_col = cur_col
        cur_col = np.zeros(len(v)+1)
    return pre_col


def Hirschberg(v, w, delta):
    def HB(i, j, k, l):
        if l-j <= 1:
            return
        mid = j + (l-j) // 2
        prefix = last_col(v[i:k], w[j:mid], delta)
        suffix = last_col(v[i:k][::-1], w[mid:l][::-1], delta)
        i_star = i + np.argmax(prefix + suffix[::-1])
        align.append((i_star, mid))
        HB(i, j, i_star, mid)
        HB(i_star, mid, k, l)

    print("The length of the input sequences are: n = ",len(v),", m = ",len(w))

    align = [(len(v), len(w))]
    HB(0, 0, len(v), len(w))
    align.sort()

    score = 0.0
    new_v = []
    new_w = []
    prev = (0, 0)
    for (i, j) in align:
        if tuple(map(lambda x, y: x-y, prev, (i, j))) == LEFT:
            score += delta['-'][w[j-1]]
            new_v.append('-')
            new_w.append(w[j-1])
        elif tuple(map(lambda x, y: x-y, prev, (i, j))) == UP:
            score += delta[v[i-1]]['-']
            new_v.append(v[i-1])
            new_w.append('-')
        elif tuple(map(lambda x, y: x-y, prev, (i, j))) == TOPLEFT:
            score += delta[v[i-1]][w[j-1]]
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        else:
            while tuple(map(lambda x, y: x-y, prev, (i, j))) != TOPLEFT:
                score += delta[v[prev[0]]]['-']
                new_v.append(v[prev[0]])
                new_w.append('-')
                prev = (prev[0]+1, prev[1])
                # print(score)
            score += delta[v[i-1]][w[j-1]]
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        prev = (i, j)
        # print(score)
    return score, ''.join(new_v)+'\n'+''.join(new_w)

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

    score, align = Hirschberg(seqA, seqB, delta)
    print(score)
    print(align)
    print()
    
