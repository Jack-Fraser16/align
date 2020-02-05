# John Fraser
# Professor Slonim
# Comp 167 - Computational Biology, Tufts University
# February 4th, 2020

# align.py reads input in the form of a FASTA-formatted file and determines
# the global alignment of the first two sequences listed within that file.    
# It then prints out how the two sequences are best aligned using penalties 
# for gaps and mismatches between the sequences. 

#!/usr/bin/env python3
import argparse
import sys

def main():

    global match
    global mismatch
    global gap
    global query 
    global subject
   
    # Tool for parsing input from the command line.  
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('-m', '--manual', type=int, nargs=3, help="First comma\
                         nd line argument sets the value of a match, second \
                         sets the penalty for a mismatch, and the third sets \
                         the penalty for a gap.")
    args = parser.parse_args()

    # Assumes default scoring system unless optional alternatives are entered.
    if args.manual:
        match = args.manual[0]
        mismatch = args.manual[1]
        gap = args.manual[2]
    else:
        match = 4
        mismatch = -2
        gap = -2

    # Instantiating key variables that are declared globally but need 
    # placeholder values before being used. Query and subject will contain
    # the sequences being compared, while curr_string and counter simply help
    # parse input.
    query = ''
    subject = ''
    curr_string = ''
    counter = 0

    # Reading in the file line by line to retrieve the first two sequences 
    # that are contained within the file. 
    with open(args.filename) as file:
        next(file)
        contents = file.readlines()

    for line in contents:
        curr_string = curr_string[:-1]
        if line[0] == '>' or line[0] == '':
            query = curr_string
            curr_string = ''
            counter += 1
            if counter == 2:
                break
        else:
            for i in line:
                curr_string = curr_string + i

    subject = curr_string

    # Calling the helper functions that instantiate and populate the necessary
    # solutions matrix using data structures defined below.  
    matrix, print_order = create_matrix(query, subject)
    subject, query = '', ''
    full_matrix = populate_matrix(matrix)
    recursive_arrow_search(full_matrix, full_matrix.horizontal - 1, 
                                                     full_matrix.vertical - 1)

    # Formatting output to meet expectations.
    if print_order == True:
        print(subject)
        print(query, end = '')
        return

    print(query)
    print(subject, end = '')
    return

# create_matrix both instantiates the matrix struct that contains the backbone 
# 2D-Array representing the solutions matrix and populates each index with the
# struct cell.
def create_matrix(seq1, seq2):

    length_s1 = int(len(seq1))
    length_s2 = int(len(seq2))

    # Determines if the order of how sequences were handed is switched, it is
    # only important for formatting output correctly. 
    print_switch = True

    # Determines which sequence is longer and ensures it to be on the 
    # horizontal axis for consistency.
    if (length_s1 > length_s2 or length_s1 == length_s2):
        empty_matrix = make_matrix(length_s1, length_s2, seq1, seq2)
    else: 
        empty_matrix = make_matrix(length_s2, length_s1, seq2, seq1)
        print_switch = False

    for i in range(empty_matrix.vertical):
        for j in range(empty_matrix.horizontal):
            empty_matrix.board[i][j] = cell(0)

    return empty_matrix, print_switch

# populate_matrix populates each cell struct within each index of the matrix
# with the appropriate value based on both the Needleman-Wunsch algorithm and
# either the default or manually-entered scoring criteria. 
def populate_matrix(matrix):
     
    for i in range(matrix.vertical):
        matrix.board[i][0].value = gap * i
        matrix.board[i][0].up = True

    for j in range(matrix.horizontal):
        matrix.board[0][j].value = gap * j
        matrix.board[0][j].left = True
    
    # Iterate through the matrix in column-first order and determine the 
    # optimal score for each cell and assign arrow flags accordingly.
    for i in range(1, matrix.vertical):
        for j in range(1, matrix.horizontal):

            move_cost = [0] * 3
            diagonal_move = match

            if (matrix.vert_seq[i - 1] != matrix.horiz_seq[j - 1]):
                diagonal_move = mismatch

            move_cost[0] = matrix.board[i - 1][j].value + gap
            move_cost[1] = matrix.board[i - 1][j - 1].value + diagonal_move   
            move_cost[2] = matrix.board[i][j - 1].value + gap

            matrix.board[i][j].value = max(move_cost)

            for l in range(3):
                if move_cost[l] == max(move_cost):
                    if l == 0: 
                        matrix.board[i][j].up = True
                    if l == 1:
                        matrix.board[i][j].diagonal = True
                    if l == 2:
                        matrix.board[i][j].left = True

    return matrix

# recursive_arrow_search is a recursive function that begins at the bottom 
# right-hand corner of the solutions matrix and searches for the optimal
# pathway to the upper left-hand corner. It modifies the global variables
# subject and query with the letters of the alignment that correspond to
# each cell it visits along the optimal pathway.  
def recursive_arrow_search(matrix, h_index, v_index):

    global query
    global subject

    if h_index == 0 and v_index == 0:
        return
    # Note that the desicion of always checking the diagonal arrow first
    # and then the leftwards arrow and finally the upwards arrow was made more
    # or less arbitrarily as any choice in which arrow to follow a) becomes 
    # meaningless with user-defined scoring schemes and b) is ultimately
    # inconsequential as both arrows in a tie would result in globally aligned
    # sequences with the same socre.
    if matrix.board[v_index][h_index].diagonal == True:
        subject = matrix.horiz_seq[h_index - 1] + subject
        query = matrix.vert_seq[v_index - 1] + query
        recursive_arrow_search(matrix, h_index - 1, v_index - 1)

    elif matrix.board[v_index][h_index].left == True:
        subject = matrix.horiz_seq[h_index - 1] + subject
        query = '-' + query
        recursive_arrow_search(matrix, h_index - 1, v_index)

    else:
        subject = '-' + subject
        query = matrix.vert_seq[v_index - 1] + query
        recursive_arrow_search(matrix, h_index, v_index - 1)

    return

# cell is the data structure that populates each index of the solutions matrix
# and is designed to keep track of the arrows and numerical value associated
# with a given cell.
class cell:
    def __init__(self, value):
        self.up = False
        self.left = False
        self.diagonal = False
        self.value = value

# make_matrix is the struct that contains the solutions matrix and information 
# about its size and the sequences it is used to compare. Note that horizontal
# is reserved for the larger of the two sequences and this design choice was 
# made purely to help me consistently visualize the matrix.
class make_matrix:
    def __init__(self, horizontal, vertical, s1, s2):
        self.horizontal = horizontal + 1
        self.vertical = vertical + 1
        self.horiz_seq = s1
        self.vert_seq = s2
        self.board = [[0 for i in range(self.horizontal)] for j in range(self.vertical)]

if __name__ == "__main__":
    main()