#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
    max_score = 0
    max_score_index = (0,0)
    alnseq1 = ""
    alnseq2 = ""

    scoring_matrix = [[None for _ in range(len(seq2)+1)] for _ in range(len(seq1)+1)]

    path_matrix = [[None for _ in range(len(seq2))] for _ in range(len(seq1))]
    for i in range(1, len(seq1)+1):
        scoring_matrix[i][0]=0
    for j in range(0, len(seq2)+1):
        scoring_matrix[0][j]=0
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            # left
            scoring_matrix[i][j] = max((scoring_matrix[i][j-1] - gapopen), 0)
            path_matrix[i-1][j-1] = 0
            # diagonal
            is_match = seq1[i-1] == seq2[j-1]
            value = scoring_matrix[i-1][j-1] + is_match*match - (not is_match)*mismatch
            if value > scoring_matrix[i][j]:
                scoring_matrix[i][j]=value
                path_matrix[i-1][j-1]=1
            # up
            value = scoring_matrix[i-1][j] - gapopen
            if value > scoring_matrix[i][j]:
                scoring_matrix[i][j]=value
                path_matrix[i-1][j-1]=2
            if max_score < scoring_matrix[i][j]:
                max_score = scoring_matrix[i][j]
                max_score_index = (i,j)
    # print(scoring_matrix)
    # print(max_score)
    # print(max_score_index)
    alnseq1, alnseq2 = traceback(scoring_matrix, path_matrix, seq1, seq2, max_score_index[0], max_score_index[1])
    
    return max_score, alnseq1, alnseq2

def traceback(score_matrix, path_matrix, seq1, seq2, i, j):
    if score_matrix[i][j] == 0:
        return ("", "")
    direction = path_matrix[i-1][j-1]
    # left
    if direction == 0:
        a,b = traceback(score_matrix, path_matrix, seq1, seq2, i, j-1)
        return(a+"-", b+seq2[j-1])
    # diagonal
    elif direction == 1:
        a,b = traceback(score_matrix, path_matrix, seq1, seq2, i-1, j-1)
        return(a+seq1[i-1], b+seq2[j-1])
    # up
    else:
        a,b = traceback(score_matrix, path_matrix, seq1, seq2, i-1, j)
        return (a+seq1[i-1], b+"-")



# smith_waterman("GGTTGACTA", "TGTTACGG", 3, 3, 2, 2)

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
