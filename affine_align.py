
def read_init(file):
    """Reads the first line of the input file"""

    f = open(file, "r")
    first_line = [int(i) for i in f.readline().strip("\n").split(" ")]
    match = first_line[0]
    mismatch = first_line[1]
    gap = first_line[2]
    extend = first_line[3]
    f.close()
    return match, mismatch, gap, extend


def read_file(file):
    """"Reads the sequences to be aligned"""

    f = open(file, "r")
    output = []
    next(f)
    for line in f:
        sequences = [i for i in line.rstrip().split("\n")]
        output.append(sequences)
    f.close()
    return output


def init_match(i, j):
    """Initializes match matrix"""

    # M(0,0) = 0
    if j == 0 and i == 0:
        return 0
    else:
        # -inf for first row/col
        if j == 0 or i == 0:
            return float('-inf')
        else:
            return 0


def init_ix(i, j, gap, extend):
    """Initializes Ix matrix"""

    # -inf for first row
    if j > 0 and i == 0:
        return float('-inf')
    else:
        if i > 0:
            return -gap - ((i - 1) * extend)
        else:
            return float('-inf')


def init_iy(i, j, gap, extend):
    """Initializes Iy matrix"""

    # -inf for first col
    if j == 0 and i > 0:
        return float('-inf')
    else:
        if j > 0:
            return -gap - ((j - 1) * extend)
        else:
            return float('-inf')


def calc_match(seq_1, seq_2, match, mismatch, i, j):
    """Calculates match based on numbers specified in the input"""

    if seq_2[j-1] == seq_1[i-1]:
        return match
    else:
        return -mismatch


def create_matrices(seq_1, seq_2, gap, extend, match, mismatch):
    """Initializes and populates the three matrices with the proper values"""

    len_i = len(seq_1) + 1
    len_j = len(seq_2) + 1

    # Initializes Match, Ix, and Iy matrices
    M = [[init_match(i, j) for j in range(0, len_j)] for i in range(0, len_i)]
    Ix = [[init_ix(i, j, gap, extend) for j in range(0, len_j)] for i in range(0, len_i)]
    Iy = [[init_iy(i, j, gap, extend) for j in range(0, len_j)] for i in range(0, len_i)]

    # Populates the matrices based on the affine gap penalty function
    for j in range(1, len_j):
        for i in range(1, len_i):
            match_result = calc_match(seq_1, seq_2, match, mismatch, i, j)
            M[i][j] = max(M[i-1][j-1] + match_result, Ix[i-1][j-1] + match_result, Iy[i-1][j-1] + match_result)
            Ix[i][j] = max(M[i-1][j] - gap, Ix[i-1][j] - extend)
            Iy[i][j] = max(M[i][j-1] - gap, Iy[i][j-1] - extend)

    return M, Ix, Iy


def traceback(seq_1, seq_2, M, Ix, Iy):
    """Performs the traceback of the matrices to find the global alignment"""

    align_1 = ''
    align_2 = ''
    i = len(seq_1)
    j = len(seq_2)

    # Move through all of the matrices to find the alignment
    while i > 0 or j > 0:

        # Find the max value in each cell amongst the matrices
        max_align = max(M[i][j], Ix[i][j], Iy[i][j])

        # Follow priority order: M(i,j) > Ix(i,j) > Iy(i,j)
        # If we match, move diagonally
        if max_align == M[i][j]:
            align_1 += seq_1[i-1]
            align_2 += seq_2[j-1]
            i -= 1
            j -= 1
        # If we put a gap in the second sequence, stay in the same column but move up
        elif max_align == Ix[i][j]:
            align_1 += seq_1[i-1]
            align_2 += '_'
            i -= 1
        # If we put a gap in the first sequence, stay in the same row but move left
        elif max_align == Iy[i][j]:
            align_1 += '_'
            align_2 += seq_2[j-1]
            j -= 1

    # Convert the newly-aligned sequences into the proper format
    aligned_1 = ''.join([align_1[j] for j in range(-1, -(len(align_1)+1), -1)])
    aligned_2 = ''.join([align_2[j] for j in range(-1, -(len(align_2)+1), -1)])

    return aligned_1, aligned_2


# Obtain values from input file
match_value, mismatch_pen, gap_pen, extend_pen = read_init('align_input.txt')
sequences = read_file('align_input.txt')
seq_1 = sequences[0][0]
seq_2 = sequences[1][0]

# Create matrices and perform traceback
M, Ix, Iy = create_matrices(seq_1, seq_2, gap_pen, extend_pen, match_value, mismatch_pen)
align_1, align_2 = traceback(seq_1, seq_2, M, Ix, Iy)

# Print alignment
print(align_1)
print(align_2)