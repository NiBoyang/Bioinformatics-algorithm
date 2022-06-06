global match
global mismatch
global gap_opening
global gap_ext
global pointer

pointer = 0
match = 5
mismatch = -4
gap_opening = -10
gap_ext = -0.5

# !!unmute to allow command line parsing
# import argparse
# parser = argparse.ArgumentParser(description='Aligning sequences...')
# parser.add_argument('seq1',action="store",help="First sequence")
# parser.add_argument('seq2',action="store",help="Second sequence")
# args = parser.parse_args()
#========================================================================
#!!Unmute to read fasta using Biopython
# from Bio import SeqIO
# def read_fasta_filename(filename):
#     with open(filename, 'r') as handle:
#         for myseq in SeqIO.parse(handle, "fasta"):
#             return str(myseq.seq)
#========================================================================

def create_matrix(rows, cols):
    """Three matrices for calculating affine gap penalties"""
    def matrix_m():
        """Main Matrix, responsible for match/mismatch"""
        my_matrix = [[0.0 for col in range(cols+1)] for row in range(rows+1)]
        return my_matrix
    def matrix_x():
        """For gap in X"""
        my_matrix = matrix_m()
        for row in range(rows+1):
            my_matrix[row][0] = -float("inf")
        return my_matrix
    def matrix_y():
        """For gap in Y"""
        my_matrix = matrix_m()
        for col in range(cols+1):
            my_matrix[0][col] = -float("inf")
        return my_matrix
    return matrix_m(), matrix_x(), matrix_y()

def calc_score(x, y):
    return match if seq1[y-1] == seq2[x-1] else mismatch

def cal_three_matrices(seq1, seq2):
    M, X, Y = create_matrix(len(seq2), len(seq1))
    for i in range(1, len(seq2)+1):
        for j in range(1, len(seq1)+1):
            X[i][j] = max(M[i-1][j] + gap_opening, X[i-1][j] + gap_ext, 0)
            Y[i][j] = max(M[i][j-1] + gap_opening, Y[i][j-1] + gap_ext, 0)
            M[i][j] = max(calc_score(i, j) + M[i - 1][j - 1], calc_score(i, j) + X[i - 1][j - 1],
                          calc_score(i, j) + Y[i - 1][j - 1], 0)
    return[M, X, Y]


def get_max(mymatrix):

    max = mymatrix[0][0]
    mrow, mcol = 0, 0
    rows = len(mymatrix)
    cols = len(mymatrix[0])

    for i in range(1, rows):
        for j in range(1, cols):
            if mymatrix[i][j] > max:
                max = mymatrix[i][j]
                mrow = i
                mcol = j
    print("max score: ", max)
    return [mrow, mcol]

def traceback(mymatrix,maxv,matrix_pointer):

    global pointer
    x=maxv[0]
    y=maxv[-1]
    sc = match if seq1[y - 1] == seq2[x - 1] else mismatch

    if matrix_pointer == 0:
        if mymatrix[0][x-1][y-1] == mymatrix[1][x-1][y-1] == mymatrix[2][x-1][y-1] == 0:
            return False

        elif mymatrix[0][x-1][y-1] + sc == mymatrix[0][x][y]:
            pointer = 0
            return [x-1, y-1]

        elif mymatrix[1][x-1][y-1] +sc == mymatrix[0][x][y]:
            pointer = 1
            return [x-1, y-1]

        elif mymatrix[2][x-1][y-1] + sc == mymatrix[0][x][y]:
            pointer = 2
            return [x-1, y-1]

    elif matrix_pointer == 1:
        if mymatrix[1][x-1][y] + gap_ext == mymatrix[1][x][y]:
            pointer = 1
            return [x-1, y]

        elif mymatrix[0][x-1][y] + gap_opening == mymatrix[1][x][y]:
            pointer = 0
            return [x-1, y]

    elif matrix_pointer == 2:
        if mymatrix[2][x][y-1] + gap_ext == mymatrix[2][x][y]:
            pointer = 2
            return [x, y-1]

        elif mymatrix[0][x][y-1] + gap_opening == mymatrix[2][x][y]:
            pointer = 0
            return [x, y-1]
    else:
        return False


def print_traceback(mymatrix):
    # this will print as expected with internal gaps
    print("Building traceback...")
    maxv = get_max(mymatrix[0])

    # traverse the matrix to find the traceback elements
    # if more than one path just pick one
    topstring = ""
    midstring = ""
    bottomstring = ""

    # pad the sequences so indexes into the sequences match the matrix indexes
    asequence1 = "#" + seq1
    asequence2 = "#" + seq2

    # add first element and go to previous
    topstring += asequence1[maxv[-1]]
    bottomstring += asequence2[maxv[0]]
    if asequence1[maxv[-1]] == asequence2[maxv[0]]:
        midstring += "|"
    else:
        midstring += " "

    # here we need to store the old position so we can track if it is an insertion of deletion
    # code assumes highest score cannot be a gap!
    # old_maxv=maxv

    # add the rest of the elements
    search = True
    while (search):
        maxv = traceback(mymatrix, maxv, pointer)
        if maxv == False or ((maxv[0] == 0) or (maxv[-1] == 0)):
            break

        # deal with single gaps or matches
        # if(old_maxv[-1]==maxv[-1]):
        if pointer == 1:
            topstring += "-"
            bottomstring += asequence2[maxv[0]]
        # elif(old_maxv[0]==maxv[0]):
        elif pointer == 2:
            # insertion or deletion
            topstring += asequence1[maxv[-1]]
            bottomstring += "-"
        else:
            # add the next element and go to previous
            topstring += asequence1[maxv[-1]]
            bottomstring += asequence2[maxv[0]]

        if (asequence1[maxv[-1]] == asequence2[maxv[0]]) & (pointer == 0):
            midstring += "|"
        elif (asequence1[maxv[-1]] != asequence2[maxv[0]]) & (pointer == 0):
            midstring += "*"
        else:
            midstring += " "
        # old_maxv = maxv
    print(topstring[::-1])
    print(midstring[::-1])
    print(bottomstring[::-1])

if __name__ == "__main__":

    seq1 = "CAACGGCATGCGCAACTTGTTGCGTA"
    seq2 = "ATGGTAGGGATTATCAACGGDGDGGAAAATCTAGCCA"

    mymatrix = cal_three_matrices(seq1, seq2)
    print_traceback(mymatrix)

    # !! unmute to use command line parsing to input sequences
    # seq1 = read_fasta_filename(args.seq1)
    # seq2 = read_fasta_filename(args.seq2)

    # !!unmute to show three matrices
    # import pandas as pd
    # pd.set_option('display.expand_frame_repr', False)
    # print(pd.DataFrame(mymatrix[0]))
    # print(pd.DataFrame(mymatrix[1]))
    # print(pd.DataFrame(mymatrix[2]))
