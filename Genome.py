import timeit
import numpy as np

'''Dictionary to store the weights of the different changes we are going to need to make
key : value -> (1st nucleotide, 2nd nucleotide)'''
S = {}

'''List to store the traceback routine, which ways we've 
gone through the cache'''
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

'''Dynamic method to compute the max alignment between two sequences'''
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

'''Recursive method to go back  through the cache C and find the solution that 
we chose, then apply it, store it, and send back to our main function'''
def Traceback(i, j, A, B, C):
    if j == 0 or i == 0:
        return 0
    '''If we went down in the cache'''
    if C[i][j] == C[i][j - 1] + S[('-', B[j])]:
        X.append(['_', '->', B[j]])
        return Traceback(i, j-1, A, B, C)
    elif C[i][j] == C[i - 1][j] + S[(A[i], '-')]:
        '''If we went to the right in the cache '''
        X.append([A[i], '->', '_'])
        return Traceback(i-1, j, A, B, C)
    else:
        '''If we moved diagonally, and we have to check whether we did a match 
        or a change to another character'''
        if A[i] == B[j]:
            X.append([A[i], '=', B[j]])
        else:
            X.append([A[i], '->', B[j]])
        return Traceback(i-1, j-1, A, B, C)

'''Dynamic method to go back through the cache C and find the solution 
that produced the maximum alignment score'''
def TracebackDP(A, B, C):
    i = len(A)-1
    j = len(B)-1

    '''While we are not at the start (the throwaway character) of either word'''
    while i > 0 or j > 0:
        if C[i][j] == C[i][j - 1] + S[('-', B[j])]:
            '''If we moved down in the cache'''
            X.append(['_', '->', B[j]])
            j = j-1
        elif C[i][j] == C[i - 1][j] + S[(A[i], '-')]:
            '''If we moved to the right in the cache'''
            X.append([A[i], '->', '_'])
            i = i-1
        else:
            '''If we moved diagonally in the cache, and checking whether 
            the characters are equal, or we're doing a substitution '''
            if A[i] == B[j]:
                X.append([A[i], '=', B[j]])
            else:
                X.append([A[i], '->', B[j]])
            i = i-1
            j = j-1

'''Method to run the timing experiment'''
def experiment():
    SETUP_CODE = '''
from __main__ import Traceback, CreateCache, MaxAlign, results,  EditCosts'''

    TEST_CODE= '''
A = '_agctgag'
B = '_aacgagg'
EditCosts()
C = CreateCache(A, B)
print(MaxAlign(A, B, C))
Traceback(len(A)-1, len(B)-1, A, B, C)
results(A, B)'''

    times = timeit.timeit(TEST_CODE, SETUP_CODE, number=1)
    print(times)

def main():
    '''Set up the dictionary with the edit costs'''
    EditCosts()

    '''Opening the relevant files, names can be changed here depending on which 
    two sequences you want to compare. Human.txt, Gorilla.txt, and Neanderthal.txt'''
    with open('Neanderthal.txt', 'r') as f:
        A = '_' + f.read().strip('\n')
    with open('Gorilla.txt', 'r') as f1:
        B = '_' + f1.read().strip('\n')

    '''Creating our cache with the dimensions of the two sequences we're working with'''
    C = CreateCache(A, B)

    '''Calculating the maximum alignment between the two strings, and printing it'''
    print(MaxAlign(A, B, C), flush=True)

    '''Tracing back through our cache to print out the changes that we made to the 
    original alignment'''
    TracebackDP(A, B, C)

    '''When we compute the traceback, we do it in reverse order. So, we need to 
    turn it around to make sense'''
    moves = X[::-1]

    '''Open the file that we want to write the results to 
    NA.txt -> Neanderthal-Gorilla comparison
    HN.txt -> Human-Neanderthal comparison
    HA.txt -> Human-Gorilla comparison'''
    with open('NA.txt', 'w') as fp:
        for i in moves:
            '''For every entry in the list of moves we made, 
            we're going to print them out, one per line'''
            fp.write(i[0])
            fp.write(i[1])
            fp.write(i[2])
            fp.write('\n')
            fp.flush()

if __name__ == '__main__':
    main()