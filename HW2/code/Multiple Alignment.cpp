#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>

using namespace std;

const int MATCH_SCORE = 1;
const int OTHER_SCORE = 0;


 
int get_column_score(char c1, char c2, char c3) {
    if (c1 != '_' && c1 == c2 && c2 == c3) {
        return MATCH_SCORE;
    }
    return OTHER_SCORE; 
}

void solve() {
    string X, Y, Z;
    if (!(cin >> X >> Y >> Z)) return;

    int N = X.length();
    int M = Y.length();
    int L = Z.length();

    vector<vector<vector<int>>> DP(
        N + 1, 
        vector<vector<int>>(
            M + 1, 
            vector<int>(L + 1, 0)
        )
    );
    vector<vector<vector<tuple<int, int, int>>>> Path(
        N + 1, 
        vector<vector<tuple<int, int, int>>>(
            M + 1, 
            vector<tuple<int, int, int>>(L + 1)
        )
    );
    
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= M; ++j) {
            for (int k = 0; k <= L; ++k) {
                if (i == 0 && j == 0 && k == 0) continue;

                int max_score = -1;
                tuple<int, int, int> best_move = {0, 0, 0};

                for (int dx = 0; dx <= 1; ++dx) {
                    for (int dy = 0; dy <= 1; ++dy) {
                        for (int dz = 0; dz <= 1; ++dz) {
                            if (dx == 0 && dy == 0 && dz == 0) continue; 

                            int prev_i = i - dx;
                            int prev_j = j - dy;
                            int prev_k = k - dz;

                            if (prev_i < 0 || prev_j < 0 || prev_k < 0) continue;

                            char c1 = (dx == 1) ? X[prev_i] : '_';
                            char c2 = (dy == 1) ? Y[prev_j] : '_';
                            char c3 = (dz == 1) ? Z[prev_k] : '_';

                            int current_score = DP[prev_i][prev_j][prev_k] + get_column_score(c1, c2, c3);

                            if (current_score > max_score) {
                                max_score = current_score;
                                best_move = {dx, dy, dz};
                            }
                        }
                    }
                }
                
                DP[i][j][k] = max_score;
                Path[i][j][k] = best_move;
            }
        }
    }

    cout << DP[N][M][L] << endl;

    string aligned_X = "", aligned_Y = "", aligned_Z = "";
    int i = N, j = M, k = L;

    while (i > 0 || j > 0 || k > 0) {
        tuple<int, int, int> move = Path[i][j][k];
        int dx = get<0>(move);
        int dy = get<1>(move);
        int dz = get<2>(move);
        
        if (dx == 1) aligned_X += X[i - 1]; else aligned_X += '_';
        if (dy == 1) aligned_Y += Y[j - 1]; else aligned_Y += '_';
        if (dz == 1) aligned_Z += Z[k - 1]; else aligned_Z += '_';
        
        i -= dx;
        j -= dy;
        k -= dz;
    }

    reverse(aligned_X.begin(), aligned_X.end());
    reverse(aligned_Y.begin(), aligned_Y.end());
    reverse(aligned_Z.begin(), aligned_Z.end());

    cout << aligned_X << endl;
    cout << aligned_Y << endl;
    cout << aligned_Z << endl;
}

int main() {
    solve();
    return 0;
}