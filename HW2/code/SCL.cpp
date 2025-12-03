#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

int get_score_strategic(char a, char b) {
    const int MATCH_SCORE = 0.1;      // m > 0
    const int MISMATCH_SCORE = -10000; // d < 0 
    
    if (a == b) {
        return MATCH_SCORE;
    } else {
        return MISMATCH_SCORE;
    }
}

void find_max_gaps(const string& seq1, const string& seq2) {
    int n = seq1.length(); 
    int m = seq2.length(); 

    const int GAP_PENALTY = -1; // g < 0

    vector<vector<int>> F(n + 1, vector<int>(m + 1));
    vector<vector<int>> P(n + 1, vector<int>(m + 1, 0)); 

    for (int i = 0; i <= n; ++i) {
        F[i][0] = i * GAP_PENALTY; 
        if (i > 0) P[i][0] = 2; // Up
    }
    for (int j = 0; j <= m; ++j) {
        F[0][j] = j * GAP_PENALTY; 
        if (j > 0) P[0][j] = 3; // Left
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int current_match_score = get_score_strategic(seq1[i - 1], seq2[j - 1]);

            int diagonal_score = F[i - 1][j - 1] + current_match_score;
            int up_score = F[i - 1][j] + GAP_PENALTY; 
            int left_score = F[i][j - 1] + GAP_PENALTY; 

            F[i][j] = max({diagonal_score, up_score, left_score});

            // ذخیره جهت
            if (F[i][j] == diagonal_score) {
                P[i][j] = 1; 
            } else if (F[i][j] == up_score) {
                P[i][j] = 2; 
            } else { 
                P[i][j] = 3; 
            }
        }
    }

    int max_gaps = 0;
    int i = n;
    int j = m;

    while (i > 0 || j > 0) {
        int direction = P[i][j];

        if (direction == 1) { 
            i--;
            j--;
        } else if (direction == 2) { 
            max_gaps++;
            i--;
        } else { 
            max_gaps++;
            j--;
        }
    }

    cout << max_gaps << endl;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string S, T;
    if (!(cin >> S >> T)) {
        return 1;
    }

    find_max_gaps(S, T);

    return 0;
}