/*
 * @Description: 
 * @Author: Zhaoxi Chen
 * @Github: https://github.com/FrozenBurning
 * @Date: 2020-03-13 23:31:23
 * @LastEditors: Zhaoxi Chen
 * @LastEditTime: 2020-03-14 10:46:25
 */
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

int NeedlemanWunsch(string seq1, string seq2, int alpha, int beta, int gamma)
{
	char gap = 0;
	// Read data
	int len1 = seq1.size();
	int len2 = seq2.size();
    // trans string into vector
	char * seq1_vec = new char [len1+1];
	char * seq2_vec = new char [len2+1];
	for (int i = 0; i < len1; i++)
    {
        seq1_vec[i] = seq1[i];
    }
	for (int i = 0; i < len2; i++)
    {
        seq2_vec[i] = seq2[i];
    }

	// Initialize score matrix
	int ** score_mat = new int * [len1+1];

	// Initialize multiple auxiliary matrix
	int ** tracex = new int * [len1+1];
	int ** tracey = new int * [len1+1];
    int ** recorder1 = new int * [len1+1];
	int ** recorder2 = new int * [len1+1];

	for (int i = 0; i < len1+1; i++)
    {
		score_mat[i] = new int [len2+1];
		score_mat[i][0] = i*gamma;
		tracex[i] = new int [len2+1];
		tracey[i] = new int [len2+1];   
        recorder1[i] = new int [len2+1];	
		recorder2[i] = new int [len2+1];	
    }

	for (int i = 0; i < len2+1; i++) 
	{
		score_mat[0][i] = i*gamma;
	}

	for (int i = 1; i < len1+1; i++)
	{
		tracex[i][0] = i-1;
		tracey[i][0] = 0;
		recorder1[i][0] = seq1_vec[i-1];
		recorder2[i][0] = gap;
	}
	for (int i = 1; i < len2+1; i++)
	{
		tracex[0][i] = 0;
		tracey[0][i] = i-1;
		recorder1[0][i] = gap;
		recorder2[0][i] = seq2_vec[i-1];
	}
	tracex[0][0] = -1;
	tracey[0][0] = -1;



    // Main Solver
	for (int i = 1; i < len1+1; i++)
	{
		for (int j =1; j < len2+1; j++)
		{
			int diag, left, up,candidate;
			// compute scores in three directions
			if (seq1_vec[i-1] == seq2_vec[j-1]) 
			{
				diag = score_mat[i-1][j-1] + alpha;
			}
			else
			{
				diag = score_mat[i-1][j-1] + beta;
			}
			left = score_mat[i-1][j] + gamma;
			up = score_mat[i][j-1] + gamma;

            // update decision 
            candidate = max(diag,max(left,up));
            if (candidate == left)
            {
                score_mat[i][j] = left;
				tracex[i][j] = i-1;
				tracey[i][j] = j;
				recorder1[i][j] = seq1_vec[i-1];
				recorder2[i][j] = gap;
            }
            else if (candidate==up)
            {
				score_mat[i][j] = up;
				tracex[i][j] = i;
				tracey[i][j] = j-1;
				recorder1[i][j] = gap;
				recorder2[i][j] = seq2_vec[j-1];                
            }
            else
            {
				score_mat[i][j] = diag;
				tracex[i][j] = i-1;
				tracey[i][j] = j-1;
				recorder1[i][j] = seq1_vec[i-1];
				recorder2[i][j] = seq2_vec[j-1];
            }
		}
	}

	int row_ind = len1;
	int col_ind = len2;
	int * backtrace_row = new int [len1+len2];
	int * backtrace_col = new int [len1+len2];

	int counter = 0;
    int tmpx,tmpy;
	backtrace_row[counter] = row_ind;
	backtrace_col[counter] = col_ind;
	while(row_ind != 0 || col_ind != 0)
	{
		counter ++;
		tmpx = tracex[row_ind][col_ind];
		tmpy = tracey[row_ind][col_ind];
		row_ind = tmpx;
		col_ind = tmpy;
		backtrace_row[counter] = row_ind;
		backtrace_col[counter] = col_ind;
	}


	FILE * f = fopen("../Answer.txt", "w");
	for (int i = counter-1; i >= 0; i--)
	{
		char output1 = recorder1[backtrace_row[i]][backtrace_col[i]];
		char output2 = recorder2[backtrace_row[i]][backtrace_col[i]];
		int same = output1*output2;
		if (same != 0)
			fprintf(f,"%c",output1);
	}
	fclose(f);

	return(score_mat[len1][len2]);
}

int main()
{
    string seq1,seq2;
    // parameters
	int alpha = 1;
	int beta = 0;
	int gamma = 0;

	ifstream LCSfile("../LongestCommonSeq.txt");
	getline(LCSfile, seq1);
	getline(LCSfile, seq2);
	LCSfile.close();

	int score = NeedlemanWunsch(seq1,seq2, alpha, beta, gamma);
	printf("Score: %d\n",score);

	return 0;
}