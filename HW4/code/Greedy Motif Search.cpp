#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

int getScore(const vector<string>& motifs) {
    if (motifs.empty()) return 0;
    int k = motifs[0].length();
    int t = motifs.size();
    int totalScore = 0;

    for (int j = 0; j < k; ++j) {
        map<char, int> counts;
        counts['A'] = 0; counts['C'] = 0; counts['G'] = 0; counts['T'] = 0;
        
        for (int i = 0; i < t; ++i) {
            counts[motifs[i][j]]++;
        }
        
        int maxCount = 0;
        for (map<char, int>::iterator it = counts.begin(); it != counts.end(); ++it) {
            if (it->second > maxCount) {
                maxCount = it->second;
            }
        }
        totalScore += maxCount; 
    }
    return totalScore;
}

map<char, vector<double> > createProfile(const vector<string>& motifs) {
    int t = motifs.size();
    int k = motifs[0].length();
    map<char, vector<double> > profile;
    
    string bases = "ACGT";
    for (int b = 0; b < 4; ++b) {
        profile[bases[b]] = vector<double>(k, 1.0);
    }

    for (int j = 0; j < k; ++j) {
        for (int i = 0; i < t; ++i) {
            profile[motifs[i][j]][j] += 1.0;
        }
        for (int b = 0; b < 4; ++b) {
            profile[bases[b]][j] /= (t + 4.0);
        }
    }
    return profile;
}

string profileMostProbable(string text, int k, map<char, vector<double> >& profile) {
    double maxProb = -1.0;
    string bestKmer = text.substr(0, k);

    for (int i = 0; i <= (int)text.length() - k; ++i) {
        string kmer = text.substr(i, k);
        double prob = 1.0;
        for (int j = 0; j < k; ++j) {
            prob *= profile[kmer[j]][j];
        }

        if (prob > maxProb) {
            maxProb = prob;
            bestKmer = kmer;
        }
    }
    return bestKmer;
}

vector<string> greedyMotifSearch(vector<string> dna, int k, int t) {
    vector<string> bestMotifs;
    for (int i = 0; i < t; ++i) {
        bestMotifs.push_back(dna[i].substr(0, k));
    }
    int bestScore = getScore(bestMotifs);

    string firstSeq = dna[0];
    for (int i = 0; i <= (int)firstSeq.length() - k; ++i) {
        vector<string> currentMotifs;
        currentMotifs.push_back(firstSeq.substr(i, k)); 

        for (int j = 1; j < t; ++j) {
            map<char, vector<double> > profile = createProfile(currentMotifs);
            currentMotifs.push_back(profileMostProbable(dna[j], k, profile));
        }

        int currentScore = getScore(currentMotifs);
        if (currentScore > bestScore) {
            bestScore = currentScore;
            bestMotifs = currentMotifs;
        }
    }
    return bestMotifs;
}

int main() {
    int k, t;
    if (!(cin >> k >> t)) return 0;

    vector<string> dna_sequences;
    for (int i = 0; i < t; ++i) {
        string temp;
        cin >> temp;
        dna_sequences.push_back(temp);
    }

    vector<string> result = greedyMotifSearch(dna_sequences, k, t);

    for (int i = 0; i < (int)result.size(); ++i) {
        cout << result[i] << endl;
    }

    return 0;
}