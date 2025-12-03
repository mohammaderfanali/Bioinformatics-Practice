#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>

using namespace std;

const int INF = 1e9;
const int GAP = -1;
const int MATCH = 1;
const int MISMATCH = -1;

int get_score(char a, char b) {
    return (a == b) ? MATCH : MISMATCH;
}


vector<int> get_alignment_scores(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    
 
    vector<int> prev(n + 1);
    vector<int> curr(n + 1);

    for (int j = 0; j <= n; ++j) {
        prev[j] = j * GAP;
    }

    for (int i = 1; i <= m; ++i) {
        curr[0] = i * GAP; 
        for (int j = 1; j <= n; ++j) {
            int match_val = prev[j - 1] + get_score(pattern[i - 1], text[j - 1]);
            int delete_val = prev[j] + GAP;   
            int insert_val = curr[j - 1] + GAP; 
            
            curr[j] = max({match_val, delete_val, insert_val});
        }
        prev = curr; 
    }
    return prev; 
}

void solve() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string S;
    int N;
    if (!(cin >> S >> N)) return;

    vector<string> W(N);
    for (int i = 0; i < N; ++i) {
        cin >> W[i];
    }

    int s_len = S.length();


    vector<int> best_prefix(s_len + 1, -INF);

    for (const string& w : W) {
        vector<int> scores = get_alignment_scores(S, w);
        for (int k = 0; k <= s_len; ++k) {
            best_prefix[k] = max(best_prefix[k], scores[k]);
        }
    }

    string rev_S = S;
    reverse(rev_S.begin(), rev_S.end());
    
    vector<int> best_suffix(s_len + 1, -INF);

    for (string w : W) {
        reverse(w.begin(), w.end());
        vector<int> scores = get_alignment_scores(rev_S, w);
        for (int k = 0; k <= s_len; ++k) {
            best_suffix[k] = max(best_suffix[k], scores[k]);
        }
    }

   
    
    long long max_total_score = -1e18; 

    for (int k = 0; k <= s_len; ++k) {
        
        if (best_prefix[k] != -INF && best_suffix[s_len - k] != -INF) {
            long long current_score = (long long)best_prefix[k] + best_suffix[s_len - k];
            if (current_score > max_total_score) {
                max_total_score = current_score;
            }
        }
    }

    cout << max_total_score << endl;
}

int main() {
    solve();
    return 0;
}