################################################################
#
# Needleman-Wunsch alignment algorithm, but enhanced by a
# reconstruction procedure that goes through all compatible solutions.
#
################################################################

penalties = { "sub" : -1,
              "gap" : -1,
              "equ" :  0 }

def sc(b1,b2):
    """Individual letter score."""
    if b1==b2:
        return penalties["equ"]
    else:
        return penalties["sub"]

def sc_gap():
    return penalties["gap"]

########################################################################

def score(g1, g2):
    f = [[0 for j in range(len(g2)+1)] for i in range(len(g1)+1)]

    for i in range(len(g1)+1):
        f[i][0] = sc_gap() * i
    for j in range(len(g2)+1):
        f[0][j] = sc_gap() * j

    for i in range(1, len(g1)+1):
        for j in range(1, len(g2)+1):
            score_both = f[i-1][j-1] + sc(g1[i-1],g2[j-1])
            f[i][j] = max(score_both,
                          f[i-1][j] + sc_gap(),
                          f[i][j-1] + sc_gap())
    
    for k in penalties:
        print(k, "penalty", "%2d" % penalties[k])
    print()

    print(g1)
    print(g2)
    print()

    # print table

    print(" " + "".join("%4c" % c for c in "_" + g2))

    for i, l in enumerate(f):
        print(("_"+g1)[i] +             # vertical print of alignment word
              "".join(["%4d" % n for n in l]))
    
    
    print()
    for align_a, align_b in reconstruct(f, g1, g2, len(g1), len(g2)):
        print(align_a)
        print(align_b)
        print()

########################################################################

def reconstruct(f, g1, g2, i, j):
    # recursive reconstruction of the alignment, going through all
    # cases

    if not i and not j:              # reached beginning of both words
        yield "", ""                 # base case

    if j-1 >= 0:
        if f[i][j] == f[i][j-1] + sc_gap(): # matches!
            for align_a, align_b in reconstruct(f, g1, g2, i, j-1):
                yield align_a + "-", align_b + g2[j-1]
    if i-1 >= 0:
        if f[i][j] == f[i-1][j] + sc_gap():
            for align_a, align_b in reconstruct(f, g1, g2, i-1, j):
                yield align_a + g1[i-1], align_b + "-"
        
    if i-1 >= 0 and j-1 >= 0:                     # diag
        if f[i][j] == f[i-1][j-1] + sc(g1[i-1],g2[j-1]):
            for align_a, align_b in reconstruct(f, g1, g2, i-1, j-1):
                yield align_a + g1[i-1], align_b + g2[j-1]

################################################################

def load_data(filename):
    f = open(filename)

    sequences = [seq for seq in f if len(seq) > 0]
    return sequences

########################################################################

def align():
    # get sequences from file
    sequences = load_data('basestring.dat')
    sequences = [s.strip() for s in sequences]

    # and run the alignment for all combinations (typically only 2)
    for i in range(1, len(sequences)):
        score(sequences[0], sequences[i])

# main prog

penalties["sub"] = -1
penalties["gap"] = -1
align()
