#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <cmath>
#include <string>


using namespace std;
const string AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";

const int BLOSUM62[20][20] = {
    //   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
    { 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2}, // A
    { 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2}, // C
    {-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3}, // D
    {-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2}, // E
    {-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3}, // F
    { 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3}, // G
    {-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2}, // H
    {-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1}, // I
    {-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2}, // K
    {-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1}, // L
    {-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1}, // M
    {-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2}, // N
    {-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3}, // P
    {-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1}, // Q
    {-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2}, // R
    { 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2}, // S
    { 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2}, // T
    { 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1}, // V
    {-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2}, // W
    {-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7}  // Y
};

int get_substitution_score(char a, char b) {
    size_t idx_a = AMINO_ACIDS.find(a);
    size_t idx_b = AMINO_ACIDS.find(b);

    if (idx_a == string::npos || idx_b == string::npos) {
        return -9999;
    }
    return BLOSUM62[idx_a][idx_b];
}



struct MaxScorePos {
    int score = 0;
    int i = 0;
    int j = 0;
};

void local_alignment_affine(const string& seq_a, const string& seq_b) {
    int n = seq_a.length();
    int m = seq_b.length();

    const int ALPHA = -11; // Gap Opening Penalty
    const int BETA = -1;   // Gap Extension Penalty
    

    vector<vector<int>> M(n + 1, vector<int>(m + 1));
    vector<vector<int>> I(n + 1, vector<int>(m + 1));
    vector<vector<int>> D(n + 1, vector<int>(m + 1));


    vector<vector<char>> PM(n + 1, vector<char>(m + 1));
    vector<vector<char>> PI(n + 1, vector<char>(m + 1));
    vector<vector<char>> PD(n + 1, vector<char>(m + 1));


    //initialize
    vector<vector<int>> P(n + 1, vector<int>(m + 1, 0));
    MaxScorePos max_score_pos = {0, 0, 0};
    for (int i = 0; i <= n; ++i) {
        M[i][0] = 0;
        I[i][0] = 0; 
        D[i][0] = 0;
        
        PM[i][0] = '\0';
        PI[i][0] = '\0'; 
        PD[i][0] = '\0';
         
    }
    for (int j = 0; j <= m; ++j) {
        M[0][j] = 0;
        I[0][j] = 0;
        D[0][j] = 0; 

        PM[0][j] = '\0';
        PI[0][j] = '\0'; 
        PD[0][j] = '\0';
    }


    //fil matrix
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int sub_score = get_substitution_score(seq_a[i - 1], seq_b[j - 1]);
            
           
            D[i][j] = max({D[i - 1][j] + BETA, 
                              M[i - 1][j] + ALPHA , 0});

            if(D[i][j] == D[i-1][j] + BETA) PD[i][j] = 'D';
            else if (D[i][j] == M[i-1][j] + ALPHA ) PD[i][j] = 'M';
            else PD[i][j] = 'R';


        
            I[i][j] = max({I[i][j - 1] + BETA,
                              M[i][j - 1] + ALPHA , 0});

            if(I[i][j] == I[i][j-1] + BETA) PI[i][j] = 'I';
            else if (I[i][j] == M[i][j-1] + ALPHA ) PI[i][j] = 'M';
            else PI[i][j] = 'R';                  

           
            M[i][j] = max({0, // Smith-Waterman MUST include 0
                                M[i - 1][j - 1] + sub_score, 
                                D[i][j], 
                                I[i][j]});
            if(M[i][j] == M[i-1][j-1]+sub_score)PM[i][j] = 'M';
            else if (M[i][j] == D[i][j])PM[i][j] = 'D';
            else if (M[i][j] == I[i][j])PM[i][j] = 'I';
            else PM[i][j] = 'R';

            if (M[i][j] > max_score_pos.score) {
                max_score_pos.score = M[i][j];
                max_score_pos.i = i;
                max_score_pos.j = j;
            }
        }
    }


    cout << max_score_pos.score << endl;
    
    
    // Start traceback
    int i = max_score_pos.i;
    int j = max_score_pos.j;

    int current_matrix = 1;
    char move;

    while (i > 0  && j > 0) {

    if(current_matrix == 1)move = PM[i][j];
    else if(current_matrix == 2)move = PD[i][j];
    else move = PI[i][j];

    if (move == 'M'){
        if (current_matrix != 3) i--;
        if(current_matrix !=2)  j --;
        current_matrix = 1; 
    }
    else if(move == 'D'){
        if(current_matrix != 1)i--;
        current_matrix = 2;
    }
    else if(move == 'I'){
        if(current_matrix != 1) j--;
        current_matrix = 3;
    }else{
        break;
    }

    }

    for(int x = i ; x < max_score_pos.i ; x++){
        cout << seq_a[x] ;
    }
    cout << endl;

    for(int y = j ; y < max_score_pos.j ; y++){
        cout << seq_b[y] ;
    }
    cout << endl;





}



int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string A, B;
    if (!(cin >> A >> B)) {
        return 1;
    }

    local_alignment_affine(A, B);

    return 0;
}