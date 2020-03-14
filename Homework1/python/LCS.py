'''
@Description: 
@Author: Zhaoxi Chen
@Github: https://github.com/FrozenBurning
@Date: 2020-03-13 23:29:13
@LastEditors: Zhaoxi Chen
@LastEditTime: 2020-03-14 10:41:09
'''
import sys
from argparse import ArgumentParser

def score_criterion(char1, char2,alpha,beta,gamma):
    '''
    @description: score computing between 2 matrix
    @param 
    @return: score
    '''
    if char1 == '-' or char2 == '-': #gap
        return gamma
    elif char1 != char2:   #mismatch
        return beta
    elif char1 == char2:   # match
        return alpha

def scoring_with_matrix(seq1, seq2,alpha,beta,gamma):
    '''
    @description: NW algorithms
    @param  
    @return: score matrix & trace matrix 
    '''
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}

    for i,char1 in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j,char2 in enumerate(seq2):
            if i == 0:                    # first row, gap in seq1
                score_mat[i][j] = j*gamma
                trace_mat[i][j] = [1]
                continue
            if j == 0:                    # first column, gap in seq2
                score_mat[i][j] = i*gamma
                trace_mat[i][j] = [2]
                continue
            ul = score_mat[i-1][j-1] + score_criterion(char1, char2,alpha,beta,gamma)     # from up-left, mark 0
            l  = score_mat[i][j-1]   + score_criterion('-', char2,alpha,beta,gamma)   # from left, mark 1, gap in seq1
            u  = score_mat[i-1][j]   + score_criterion(char1, '-',alpha,beta,gamma)   # from up, mark 2, gap in seq2
            chosen = max([ul,l,u])
            score_mat[i][j] = chosen
            trace_mat[i][j] = [ind for ind,s in enumerate([ul,l,u]) if s==chosen]
            # trace_mat[i][j] = [ul, l, u].index(chosen)   # record which direction
    return score_mat, trace_mat



pathcode=[]
def traceback(seq1, seq2, trace_mat,path_code=''):
    '''
    @description: traceback optimal path
    @param trace_mat: generated from scoring_with_matrix
    @return: path_code
    '''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    # path_code = ''
    if i>0 or j>0:
        direction = trace_mat[i][j] #supposed to be list
        for ind,direct in enumerate(direction):
            tmp_i=i
            tmp_j=j
            if direct == 0:                    # from up-left direction
                tmp_i = i-1
                tmp_j = j-1
                tmppath_code = '0' + path_code
            elif direct == 1:                  # from left
                tmp_j = j-1
                tmppath_code = '1' + path_code
            elif direct == 2:                  # from up
                tmp_i = i-1
                tmppath_code = '2' + path_code

            if tmp_i >0 or tmp_j>0:
                traceback(seq1[1:tmp_i+1],seq2[1:tmp_j+1],trace_mat,tmppath_code)
            else:
                pathcode.append(tmppath_code)            
    else:
        # pathcode.append(path_code)
        return 

def print_m(seq1, seq2, m):
    """print score matrix or trace matrix"""
    seq1 = '-' + seq1; seq2 = '-' + seq2
    print()
    print(' '.join(['%3s' % i for i in ' '+seq2]))
    for i, p in enumerate(seq1):
        line = [p] + [m[i][j] for j in range(len(seq2))]
        print(' '.join(['%3s' % i for i in line]))
    print()
    return

def pretty_print_align(seq1, seq2, path_code):
    '''
    return pair alignment result string from
    path code: 0 for match, 1 for gap in seq1, 2 for gap in seq2
    '''
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '0':
            align1 = align1 + seq1[0]
            align2 = align2 + seq2[0]
            if seq1[0] == seq2[0]:
                middle = middle + '|'
            else:
                middle = middle + ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        elif p == '1':
            align1 = align1 + '-'
            align2 = align2 + seq2[0]
            middle = middle + ' '
            seq2 = seq2[1:]
        elif p == '2':
            align1 = align1 + seq1[0]
            align2 = align2 + '-'
            middle = middle + ' '
            seq1 = seq1[1:]

    print('Alignment:\n\n   ' + align1 + '\n   ' + middle + '\n   ' + align2 + '\n')
    return


def main():
    parser = ArgumentParser(description='Needleman Wunsch Algorithm')
    # parser.add_argument('--seq1',type=str,required=True,help='input sequence 1')
    # parser.add_argument('--seq2',type=str,required=True,help='input sequence 2')
    parser.add_argument('--seq1',type=str,required=False,default='AAGC',help='input sequence 1')
    parser.add_argument('--seq2',type=str,required=False,default='AGT',help='input sequence 2')
    parser.add_argument('-a','--alpha',type=int,required=False,default=1)
    parser.add_argument('-b','--beta',type=int,required=False,default=-1)
    parser.add_argument('-g','--gamma',type=int,required=False,default=-1)
    parser.add_argument('-f','--file',type=str,required=False,default='',help='input test file including 2 strings')
    args = parser.parse_args()
    if not len(args.file) == 0: 
        with open(args.file, 'r') as f:
            strings = f.read().splitlines()
        seq1 = strings[0]
        seq2 = strings[1]
    else:
        try:
            seq1,seq2 = args.seq1,args.seq2
        except:
            print("Invalid Inputs")
    print('current params: alpha=%.1f,beta=%.1f,gamma=%.1f' % (args.alpha, args.beta,  args.gamma))
    print('sequence 1: %s' % seq1)
    print('sequence 2: %s' % seq2)
    # print('args:',args.gamma)
    score_mat, trace_mat = scoring_with_matrix(seq1, seq2,args.alpha,args.beta,args.gamma)
    print_m(seq1, seq2, score_mat)
    print_m(seq1, seq2, trace_mat)

    traceback(seq1, seq2, trace_mat)
    
    for i in range(0,len(pathcode)):
        pretty_print_align(seq1, seq2, pathcode[i])
    print(pathcode)
    # print('   '+path_code)

    print('total score:',score_mat[len(seq1)][len(seq2)])

if __name__ == '__main__':
    main()
