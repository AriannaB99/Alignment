'''This was a development file, the final work should all be in the
Genome.py file. This was a file for experimentation. '''

import timeit
import numpy as np
'''Dictionary to store the weights of the different changes we are going to need to make
key : value -> (1st nucleotide, 2nd nucleotide)'''
S = {}

'''List to store the traceback routine, which ways we've 
gone through the cache
D -> backwards diagonal, L -> back to the left j-1, U -> above i-1'''
X = []

'''Method to create and return our cache'''
def CreateCache(A, B):
    m = len(B)
    n = len(A)
    '''Creating our cache
    2D array, the major axis is the place in A and the minor axis is the place in B '''
    C = np.zeros((n,m))
    return C

'''Method to construct the dictionary containing edit weights from the provided table'''
def EditCosts():
    '''Easiest steps = most positive point/effort value'''
    S[('a', 'a')] = 5
    S[('a', 'c')] = -1
    S[('a', 'g')] = -2
    S[('a', 't')] = -1
    S[('a', '-')] = -3
    S[('c', 'a')] = -1
    S[('c',  'c')] = 5
    S[('c', 'g')] = -3
    S[('c', 't')] = -2
    S[('c', '-')] = -4
    S[('g', 'a')] = -2
    S[('g', 'c')] = -3
    S[('g', 'g')] = 5
    S[('g', 't')] = -2
    S[('g', '-')] = -2
    S[('t', 'a')] = -1
    S[('t', 'c')] = -2
    S[('t', 'g')] = -2
    S[('t', 't')] = 5
    S[('t', '-')] = -1
    S[('-', 'a')] = -3
    S[('-', 'c')] = -4
    S[('-', 'g')] = -2
    S[('-', 't')] = -1

'''Dynamic method to compute the minimum edit distance between two words
A -> 1st word
B -> 2nd word (misspelled word)'''
def MaxAlign(A, B, C):
    m = len(B)
    n = len(A)
    '''Fill in our base cases, when we are at the _ of either word'''
    C[0][0] = 0

    for i in range(1, n):
        C[i][0] = S[(A[i], '-')] + C[i-1][0]
    for j in range(1,  m):
        C[0][j] = S[('-', B[j])] + + C[j-1][0]

    '''Compute the three options from our previous computations, and then fill in the maximum value'''
    for i in range(1, n):
        for j in range(1, m):
            x = C[i][j-1] + S[('-', B[j])]
            y = C[i-1][j] + S[(A[i], '-')]
            z = C[i-1][j-1] + S[(A[i], B[j])]
            C[i][j] = max(x, y, z)
    '''Return our final solution, -1 because of the _ included in the length of the string'''
    return C[n-1][m-1]

'''Method to go back  through the cache C and find the solution that 
we chose, then apply it, store it, and send back to our main function'''
def Traceback(i, j, A, B, C):
    if j == 0 or i == 0:
        return 0

    if C[i][j] == C[i][j - 1] + S[('-', B[j])]:
        X.append(['_', '->', B[j]])
        return Traceback(i, j-1, A, B, C)
    elif C[i][j] == C[i - 1][j] + S[(A[i], '-')]:
        X.append([A[i], '->', '_'])
        return Traceback(i-1, j, A, B, C)
    else:
        if A[i] == B[j]:
            X.append([A[i], '=', B[j]])
        else:
            X.append([A[i], '->', B[j]])
        return Traceback(i-1, j-1, A, B, C)

def experiment():
    SETUP_CODE = '''
from __main__ import Traceback, CreateCache, MaxAlign, results,  EditCosts'''

    TEST_CODE= '''
A = '_AGTATC'
B = '_GTACA'
EditCosts()
C = CreateCache(A, B)
print(MaxAlign(A, B, C))
Traceback(len(A)-1, len(B)-1, A, B, C)
results(A, B)'''

    times = timeit.timeit(TEST_CODE, SETUP_CODE, number=1)
    print(times)

def main():
    #experiment()
    A = '_aaagct'
    B = '_cgtacg'
    EditCosts()
    C = CreateCache(A, B)
    print(MaxAlign(A, B, C))
    Traceback(len(A)-1, len(B)-1, A, B, C)
    moves = X[::-1]
    for i in moves:
        print(i[0], end = "  ")
        print(i[1], end = "  ")
        print(i[2])



if __name__ == '__main__':
    main()