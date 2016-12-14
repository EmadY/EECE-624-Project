#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <random>
#include <stack>

using namespace std;

// needed structs
struct TrieNode{
    TrieNode* child[2] = {nullptr, nullptr};
    int count = 0;
};

// function signatures
inline bool can_equate(string& s1, string& s2);
inline string create_equal(string& s1, string& s2);

void reduce_to_unique(vector<string>&);
void join2_reduce(vector<string>&);
void trie_reduce(vector<string>&);


void delete_trie(TrieNode*);

int main()
{
    std::srand ( unsigned ( std::time(0) ) );
    int c; cin >> c;
    vector<string> V(c);
    for(int i = 0; i < c; i++) cin >> V[i];


    reduce_to_unique(V);
//    for(int i = 0; i < V.size(); i++) cout << V[i] << endl;
    cout << endl;

    join2_reduce(V);

    cout << endl << "Reduced List by join2: " << endl;
    for(int i = 0; i < V.size(); i++) cout << V[i] << endl;


    return 0;
}

inline bool can_equate(string& s1, string& s2){
    if(s1.size() != s2.size()) return false;
    for(int cei = 0; cei < s1.size(); cei++)
        if((s1[cei] == '1' && s2[cei] == '0')
         || (s1[cei] == '0' && s2[cei] == '1')) return false;
     return true;
}

inline string create_equal(string& s1, string& s2){
    string create_equal_string = s1;
    for(int i = 0; i < s1.size(); i++)
        if(create_equal_string[i] == 'X' && s2[i] != 'X')
            create_equal_string[i] = s2[i];
    return create_equal_string;
}

void reduce_to_unique(vector<string>& V){
    set<string> S;
    for(int i = 0; i< V.size(); i++) S.insert(V[i]);
    V = vector<string>(S.size());
    int ind = 0;
    for(auto it = S.begin(); it != S.end(); ++it) V[ind++] = *it;
}

void join2_reduce(vector<string>& V){
    /*
     * Works by taking all N elemnts, and trying to join any two together.
     * To avoid unintentional bias, everytime we find an answer, we shuffle
     * the now reduced vector and go aagin.
     */

    vector<string> next = V;
    vector<string> curr;
    while(true){
        // current is equal to previous next
        curr = vector<string>(next);
        // move through all possibilities while we haven't removed one
        for(int i = 0; i < curr.size() && (next.size() == curr.size()); i++){
            for(int j = i+1; j < curr.size(); j++){
                if(can_equate(curr[i], curr[j])){
                    cout << "Detected " << curr[i] << " " << curr[j] << " equates to " << create_equal(curr[i], curr[j]) << endl;
                    next[i] = create_equal(curr[i], curr[j]);
                    next.erase(next.begin() + j);
                    break;
                }
            }
        }

        // if we couldn't remove anything, we're done
        if(curr.size() == next.size()) break;

        // reduced by one. shuffle next and move.
        random_shuffle(next.begin(), next.end());
    }
    V = curr;
}

void trie_reduce(vector<string>&){
    mt19937 engine(time(0));
    uniform_int_distribution<> dist;

    //TODO(emad): do this.

}




void delete_trie(TrieNode* head){
    stack<TrieNode*> rem; rem.push(head);
    while(!rem.empty()){
        TrieNode* tn = rem.top(); rem.pop();
        if(tn->child[0] != nullptr) rem.push(tn->child[0]);
        if(tn->child[1] != nullptr) rem.push(tn->child[1]);
        delete tn;
    }
}
