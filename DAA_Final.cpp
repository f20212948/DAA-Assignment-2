#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>

using namespace std;

// Function to check if two characters form a valid RNA base pair
bool isValidPair(char a, char b){
    return (a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
}

// Function to calculate the maximum number of base pairs in the RNA sequence
vector<vector<int>> numBasePairs(string rna){
    int n = rna.size();
    vector<vector<int>> dp(n , vector<int>(n , 0)); // Dynamic programming table to store intermediate results
    
    // Initialization of dynamic programming table
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
           int length = j - i;
            if(length <= 4){
                dp[i][j] = 0;
            }
        }
    }

    // Filling the dynamic programming table using bottom-up approach
    for(int k = 1; k < n; k++){
        for(int i = 0; i < n - k; i++){ // Loop over all possible substrings of the RNA sequence
            int j = i + k;
            if(k <= 4){
                dp[i][j] = 0;
                continue;
            }
            int temp = dp[i][j - 1]; // Store the maximum number of base pairs excluding the current nucleotide
            int maxi = 0;
            
            for(int t = i; t < j - 4; t++){ // Loop over all possible positions to split the substring
                if(isValidPair(rna[t], rna[j])){ // Check if the current pair of nucleotides forms a valid base pair
                    if(t == 0){
                        maxi = max(maxi, 0 + dp[t + 1][j - 1] + 1); // Update the maximum number of base pairs
                        continue;
                    }
                    maxi = max(maxi , dp[i][t - 1] + dp[t + 1][j - 1] + 1); // Update the maximum number of base pairs
                }
            }
            dp[i][j] = max(maxi,temp); // Update the value in the dynamic programming table
        }
    }

    return dp; // Return the dynamic programming table
}

// Function to print the indices of base pairs in the RNA sequence
vector<pair<int,int>> printPairs(string rna , vector<vector<int>> dp){
    stack<pair<int, int>> s;
    vector<pair<int, int>> pairs;
    s.push({0, rna.size() - 1});
    while(!s.empty()){
        pair<int, int> p = s.top();
        s.pop();
        int i = p.first;
        int j = p.second;
        if(i >= j){
            continue;
        }
        if(j > 0 && dp[i][j] == dp[i][j - 1]){
            s.push({i, j - 1});
        }
        else if(i < dp.size() - 1 && dp[i][j] == dp[i + 1][j]){
            s.push({i + 1, j});
        }
        else if(i < dp.size() - 1 && j > 0 && dp[i][j] == dp[i + 1][j - 1] +1 && isValidPair(rna[i], rna[j])){
            pairs.push_back({i, j}); // Store the indices of the base pairs
            s.push({i + 1, j - 1});
        }
        else{
            for(int k = i + 1; k < j; k++){
                if(k < dp.size() - 1 && (dp[i][k] + dp[k + 1][j] == dp[i][j])){
                    s.push({i, k});
                    s.push({k + 1, j});
                    break;
                }
            }
        }
    }

    for(auto it : pairs){
       cout << it.first << " " << it.second << endl; // Print the indices of base pairs
   }

    return pairs; // Return the indices of base pairs
}

int main(){

    string rna;

    ifstream in;
    in.open("sample.txt");

    in >> rna;

    in.close();
    cout<< rna.length() << endl;
    auto start = chrono::high_resolution_clock::now();  // Start time
    // cout << rna << endl;
    vector<vector<int>> ans = numBasePairs(rna); // Calculate the maximum number of base pairs
    auto end = chrono::high_resolution_clock::now();  // End time
    cout << ans[0][rna.size()-1] << " pairs" << endl;
    cout << rna << endl;
    
    // Calculate and print the time taken by the function
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
    vector<pair<int , int>> pairs = printPairs(rna , ans); // Print the indices of base pairs

    string dotBracket = "";
    for(int i = 0; i < rna.size(); i++){
        dotBracket += ".";
    }

    // Construct the dot-bracket notation representing the secondary structure of the RNA
    for(auto it : pairs){
        dotBracket[it.first] = '(';
        dotBracket[it.second] = ')';
    }

    cout << dotBracket << endl;

    ofstream out;
    out.open("output.txt");

    out << dotBracket;

    execlp("python", "python", "visual.py", NULL); // Execute the Python script for visualization

    return 0;
}
