#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

int get_score(char a, char b) {
    static const map<pair<char, char>, int> score_matrix = {
        {{'A', 'A'}, 3}, {{'A', 'T'}, -2}, {{'A', 'C'}, -1}, {{'A', 'G'}, 0},
        {{'T', 'A'}, -2}, {{'T', 'T'}, 3}, {{'T', 'C'}, 0}, {{'T', 'G'}, -1},
        {{'C', 'A'}, -1}, {{'C', 'T'}, 0}, {{'C', 'C'}, 3}, {{'C', 'G'}, -2},
        {{'G', 'A'}, 0}, {{'G', 'T'}, -1}, {{'G', 'C'}, -2}, {{'G', 'G'}, 3}
    };
    
    return score_matrix.at({a, b});
}




void needleman_wunsch(const string& seq1, const string& seq2, int gap_penalty) {
    int n = seq1.length();
    int m = seq2.length(); 

    vector<vector<int>> score(n + 1, vector<int>(m + 1));
    vector<vector<int>> path(n + 1, vector<int>(m + 1, 0));

    //initialize
    for (int i = 0; i <= n; ++i) {
        score[i][0] = i * gap_penalty;
        if (i > 0) path[i][0] = 2; 
    }
    for (int j = 0; j <= m; ++j) {
        score[0][j] = j * gap_penalty;
        if (j > 0) path[0][j] = 3; 
    }


    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int match_score = get_score(seq1[i - 1], seq2[j - 1]);

            int diagonal_score = score[i - 1][j - 1] + match_score;
            int up_score = score[i - 1][j] + gap_penalty; 
            int left_score = score[i][j - 1] + gap_penalty; 

            score[i][j] = max({diagonal_score, up_score, left_score});

            //save path
            if (score[i][j] == diagonal_score) {
                path[i][j] = 1; //diagonal
            } else if (score[i][j] == up_score) {
                path[i][j] = 2; //up
            } else { 
                path[i][j] = 3; // Left
            }
        }
    }


    //print
    cout << score[n][m] << endl;
    string aligned_seq1 = "";
    string aligned_seq2 = "";
    int i = n;
    int j = m;

    while (i > 0 || j > 0) {
        int direction = path[i][j];

        if (direction == 1) { 
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += seq2[j - 1];
            i--;
            j--;
        } else if (direction == 2) { 
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += '_'; 
            i--;
        } else { 
            aligned_seq1 += '_'; 
            aligned_seq2 += seq2[j - 1];
            j--;
        }
    }


    reverse(aligned_seq1.begin(), aligned_seq1.end());
    reverse(aligned_seq2.begin(), aligned_seq2.end());

    cout << aligned_seq1 << endl;
    cout << aligned_seq2 << endl;

}



int main() {
    string X, Y;
    if (!(cin >> X >> Y)) {
        return 1;
    }


    needleman_wunsch(X, Y, -1);

    return 0;
}

