#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <random>
#include <stack>
#include <functional>
#define vs std::vector<string>

using namespace std;
using std::placeholders::_1;

// needed structs
struct TrieNode{
    TrieNode* child[2] = {nullptr, nullptr};
    int count = 0;
};

// function signatures
inline bool can_equate(const string& s1, const string& s2);
inline string create_equal(const string& s1, const string& s2);

void brute_force_reduce(vs&);
void join2_reduce(vs&);
void trie_reduce(vs&); // TODO(emad): implement

void reduce_to_unique(vs&);
void seperate_vectors(const vs&, vs&, vs&);
bool check_correct_reduce(const vs&, const vs&);

void print_vector(const vs&);
void delete_trie(TrieNode*);

int main() {
    std::srand ( unsigned ( std::time(0) ) );

    int c; cin >> c; vs V(c);
    for(int i = 0; i < c; i++) cin >> V[i];
    reduce_to_unique(V);

    vector< function<void(vs&)> > functions;
    vector< vs > vectors;
    vector<bool> check_reduce;

    // For each method, add a copy of V to vectors (vs(V)), and a boolean that
    // controls whether this method will be tested. Finally add the function to
    // test using the same syntax as below.
    // Set the boolean of the function you want to test to true and rest to false,
    // unless you want to compare.

    // join2_reduce
    vectors.push_back(vs(V)); check_reduce.push_back(true);
    functions.push_back(bind(&join2_reduce, _1));

    // trie_reduce
    vectors.push_back(vs(V)); check_reduce.push_back(false); // not implemented. don't set to true.
    functions.push_back(bind(&trie_reduce, _1));



    bool display_results = true; // set to true if you want to see the reduced set.
    for(int i = 0; i < functions.size(); i++){
        if(!check_reduce[i]) continue;
        clock_t begin_clock = clock();
        functions[i](vectors[i]);
        clock_t end_clock = clock();
        bool correct_answer = check_correct_reduce(vectors[i], V);
        cout << "Function " << i+1 << " gave a " <<  (correct_answer ? "correct" : "wrong")
             << " reduced set. The reduced set had a size of " << vectors[i].size() << "."
             <<" The code ran in " << double(end_clock - begin_clock) / CLOCKS_PER_SEC
             << "s." << endl;
        if(display_results){
            cout << endl << "Reduced set is:" << endl;
            print_vector(vectors[i]);
        }
    }

    return 0;
}

/*
Returns if two vectors can be equal by setting don't cares correctly.
*/
inline bool can_equate(const string& s1, const string& s2){
    if(s1.size() != s2.size()) return false;
    for(int cei = 0; cei < s1.size(); cei++)
        if((s1[cei] == '1' && s2[cei] == '0')
         || (s1[cei] == '0' && s2[cei] == '1')) return false;
     return true;
}

/*
Returns a vector that is in fact equal to both vectors, but with more
don't cares set. The vector returned has the minimum number of don't cares removed.

Assumes that can_equate(s1, s2) is true, otherwise there is a risk of undefined behaviour.
*/
inline string create_equal(const string& s1, const string& s2){
    string create_equal_string = s1;
    for(int i = 0; i < s1.size(); i++)
        if(create_equal_string[i] == 'X' && s2[i] != 'X')
            create_equal_string[i] = s2[i];
    return create_equal_string;
}

/*
Works by taking all N elemnts, and trying to join any two together.
To avoid unintentional bias, everytime we find an answer, we shuffle
the now reduced vector and go aagin.
*/
void join2_reduce(vs& V){


    vs next = V;
    vs curr;
    while(true){
        // current is equal to previous next
        curr = vs(next);
        // move through all possibilities while we haven't removed one
        for(int i = 0; i < curr.size() && (next.size() == curr.size()); i++){
            for(int j = i+1; j < curr.size(); j++){
                if(can_equate(curr[i], curr[j])){
                    //cout << "Detected " << curr[i] << " " << curr[j] << " equates to " << create_equal(curr[i], curr[j]) << endl;
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

/*
Start with all elements that do not contain don't cares. Add them to a trie.
Go over all the rest, moving through the trie. When we have a don't care,
select a direction randomly based on how many there are in each one.
when finished, go over the trie and get the values that remained.
*/
void trie_reduce(vs&){
    mt19937 engine(time(0));
    uniform_real_distribution<> dist; //dist(engine) gives a random int




    //TODO(emad): do this.

}

/*
Takes a list of vectors and removes duplicate entries. The returned vector is
also sorted lexicographically (side-effect).
*/
void reduce_to_unique(vs& V){
    set<string> S;
    for(int i = 0; i< V.size(); i++) S.insert(V[i]);
    V = vs(S.size());
    int ind = 0;
    for(auto it = S.begin(); it != S.end(); ++it) V[ind++] = *it;
}

/*
Start with all elements that do not contain don't cares. Add them to a trie.
Go over all the rest, moving through the trie. When we have a don't care,
select a direction randomly based on how many there are in each one.
when finished, go over the trie and get the values that remained.
*/
void seperate_vectors(const vs& data, vs& known, vs& unknown) {
    known = vs();
    unknown = vs();
    for(int i = 0; i < data.size(); i++){
        for(int j = 0; j < data[i].size(); j++) {
            if(data[i][j] == 'X')
            {
                unknown.push_back(data[i]);
                break;
            }
        }
        if((i+1) != known.size() + unknown.size()) known.push_back(data[i]);
    }
}

/*
Given two lists of vectors, this method returns whether each of the vectors
in the first list is in fact a vector in the second with some don't cares set.

Used for testing purposes. O(n^2) complexity
*/
bool check_correct_reduce(const vs& initial, const vs& reduced){
    for(int i = 0; i < initial.size(); i++){
        for(int j = 0; j < reduced.size()+1; j++){
            if(j == reduced.size()) return false;
            if(can_equate(initial[i], reduced[j]) && (reduced[j] == create_equal(initial[i], reduced[j])))
                break;
        }
    }
    return true;
}

/*
Used for testing, mainly.
*/
void print_vector(const vs& V){
    for(int i = 0; i < V.size(); i++) cout << V[i] << endl;
    cout << endl;
}

/*
Deletes a trie given the trie head and avoids dangling pointers.
Note that this assumes that the trie is valid and will crash otherwise.
*/
void delete_trie(TrieNode* head){
    stack<TrieNode*> rem; rem.push(head);
    while(!rem.empty()){
        TrieNode* tn = rem.top(); rem.pop();
        if(tn->child[0] != nullptr) rem.push(tn->child[0]);
        if(tn->child[1] != nullptr) rem.push(tn->child[1]);
        delete tn;
    }
}
