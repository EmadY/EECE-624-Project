#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <random>
#include <stack>
#include <functional>
#include <fstream>
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
void mult_join2_reduce(vs&);
void trie_reduce(vs&); // TODO(emad): implement

void reduce_to_unique(vs&);
void seperate_vectors(const vs&, vs&, vs&);
bool check_correct_reduce(const vs&, const vs&);
int get_min_cover(const vs&, vs&, int);
void greddy_reduce(vs&);
int get_min(std::vector<int>&);
void print_vector(const vs&);
void delete_trie(TrieNode*);
inline bool bit(int, int);

const bool USE_BRUTE_FORCE = true; // minimize start time when brute force not required
const int DP_LEN = USE_BRUTE_FORCE ? 33554432 : 10;
int DP[DP_LEN]; // dp for brute force method
vs DPvs[DP_LEN];

int main() {
    // Note: to reduce startup time when you DO NOT want to test brute force
    // set USE_BRUTE_FORCE to false above.

    std::srand ( unsigned ( std::time(0) ) );

    bool custom_in = false; // change if you want to test on a specific vector

    vector< vs > V;

    if(custom_in){
        int c; cin >> c; V.push_back(vs(c));
        for(int i = 0; i < c; i++) cin >> V[0][i];
    } else {
        ifstream in; in.open("data.txt");
        while(true){
            int c; in >> c; if(c == 0) break;
            V.push_back(vs(c));
            for(int i = 0; i < c; i++) in >> V[V.size()-1][i];
        }
    }

    for(int vi = 0; vi < V.size(); vi++){

        reduce_to_unique(V[vi]);

        vector< function<void(vs&)> > functions;
        vector< vs > vectors;
        vector<bool> check_reduce;
        vs function_name;

        // For each method, add a copy of V to vectors (vs(V)), and a boolean that
        // controls whether this method will be tested. Finally add the function to
        // test using the same syntax as below.
        // Set the boolean of the function you want to test to true and rest to false,
        // unless you want to compare.

        // brute_force_reduce
        vectors.push_back(vs(V[vi])); check_reduce.push_back(true);
        functions.push_back(bind(&brute_force_reduce, _1));
        function_name.push_back("Brute Force");

        // join2_reduce
        vectors.push_back(vs(V[vi])); check_reduce.push_back(true);
        functions.push_back(bind(&join2_reduce, _1));
        function_name.push_back("Join2 Reduce");

        //mult_join2_reduce
        vectors.push_back(vs(V[vi])); check_reduce.push_back(true);
        functions.push_back(bind(&mult_join2_reduce, _1));
        function_name.push_back("Multiple Join2");


        // trie_reduce
        vectors.push_back(vs(V[vi])); check_reduce.push_back(false); // not implemented. don't set to true.
        functions.push_back(bind(&trie_reduce, _1));
        function_name.push_back("Trie Reduce");


        bool display_results = false; // set to true if you want to see the reduced set.
        for(int i = 0; i < functions.size(); i++){
            if(!check_reduce[i]) continue;
            clock_t begin_clock = clock();
            functions[i](vectors[i]);
            clock_t end_clock = clock();
            bool correct_answer = check_correct_reduce(V[vi], vectors[i]);
            cout << "Function :'" << function_name[i] << "'. Verdict : " <<  (correct_answer ? "Yes" : "No")
                 << ". Size : " << vectors[i].size() << ". Time: "
                 << double(end_clock - begin_clock) / CLOCKS_PER_SEC << "s." << endl;
            if(display_results){
                cout << endl << "Reduced set is:" << endl;
                print_vector(vectors[i]);
            }
        }

        cout << endl;
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
Brute force method. Will find optimal set, but will take a long time.
complexity is n^2*2^n.
Requires (and is mainly implemented in) helper function get_min_cover()
*/
void brute_force_reduce(vs& V){
    // more than 20 will take too much time, so just skip.
    if(V.size() > 25) {V = vs(); return;}

    int start_dpv = 0;
    for(int i = 0; i < V.size(); i++) start_dpv = start_dpv*2 + 1;

    fill(&DP[0], &DP[0] + start_dpv + 1, -1);


    vs Vc(V);
    get_min_cover(Vc, V, start_dpv);
}

/*
Gets the minimum number of vectors required to cover all the input vectors.
The input vectors are determined by V and dpv:
    if bit i of dpv is 1, then V[i] is in the input, otherwise it's not.
the cover is returned in cover.
Uses memoization. memoized in DP.
*/
int get_min_cover(const vs& V, vs& cover, int dpv){
    //cout << "called: " << dpv << endl;
    if(dpv == 0) {cover = vs(); return 0;}
    if(DP[dpv] != -1) {cover = DPvs[dpv]; return DP[dpv]; }
    int rem_vs = 0;
    cover = vs();
    for(int i = 0; i < V.size() && rem_vs < 2; i++)
    {
        if(!bit(dpv, i)) continue;
        rem_vs++;
        cover.push_back(V[i]);
    }
    if(rem_vs == 1) return 1;
    cover = vs();

    // Get all pairs of vectors that 'can_equate'. Try all other vectors and see
    // if they can be added to the two we found above. When done, remove from dpv
    // and call with new dpv.
    //
    // Use the set that returned the minimum dpv when returning (vs = created set +
    // vs returned by next call).
    //
    // If we can not find any pairs, then we are done, return ans = number of input
    // vectors, and vs = input vectors.
    //
    // Base case was handled above: 0 or 1 in input only remaining.

    bool found = false;
    vs best = vs(); int min_v = 1000; string best_eq_str;

    for(int i = 0; i < V.size(); i++) {
        if(!bit(dpv,i)) continue;
        for(int j = i+1; j < V.size(); j++) {
            if(!bit(dpv,j)) continue;
            string eq_string;
            if(!can_equate(V[i], V[j])) continue; // we only want compatible strings
            eq_string = create_equal(V[i], V[j]);

            int to_sub = (1 << i) + (1 << j); // to be able to change dpv

            // loop over all other numbers, see if we can find any that want to join.
            for(int k = 0; k < V.size(); k++){
                if(k == i || k == j) continue;
                if(!bit(dpv, k)) continue;
                if(!can_equate(V[k], eq_string)) continue;
                eq_string = create_equal(V[k], eq_string);
                to_sub += (1 << k);
            }

            vs ret_cover;
            int ret_cover_size = get_min_cover(V, ret_cover, dpv-to_sub);
            if(!found) found = true;
            if(ret_cover_size < min_v){
                best = ret_cover; min_v = ret_cover_size; best_eq_str = eq_string;
            }
        }
    }

    if(!found) {
        for(int i = 0; i < V.size(); i++) if(bit(dpv, i)) best.push_back(V[i]);
    } else {
        best.push_back(best_eq_str);
    }

    cover = best;
    DPvs[dpv] = cover;
    DP[dpv] = cover.size();
    //cout << "returned: " << dpv << " w/ result: " << cover.size() << endl;

    return DP[dpv];
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
                    //cout << "Detected " << curr[i] << " " << curr[j] << " equates to " 
                    //<< create_equal(curr[i], curr[j]) << endl;
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
* similiar concept to join2 reduce, except that
* it is done multiuple times with different starting permutations, as the position
* and identity of starting elemnts directly affects the output. Maintains which
* permutation lead to best output and uses that.
*/
void mult_join2_reduce(vs& V){
    vs best(V.size());
    for(int i = 0; i < V.size(); i++){
        random_shuffle(V.begin(), V.end());
        vs Vc(V);
        join2_reduce(Vc);
        if(Vc.size() < best.size()) best = Vc;
    }
    V = best;
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
            //if(j == reduced.size()) cout << "the breaker " << initial[i] << endl;
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

void greddy_reduce(vs& V) {
    std::vector<std::vector<int>> adj(V.size());
    for (int i = 0; i < V.size(); i++) {
        std::vector<int> temp;
        adj[i] = temp;
    }
    for (int i = 0; i < V.size(); i++) {
        for (int j = i+1; j < V.size(); j++) {
            if (can_equate(V[i], V[j])) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }
    std::vector<int> degrees(V.size());
    for (int i = 0; i < V.size(); i++) {
        degrees[i] = adj[i].size();
    }
    int n = V.size();
    unordered_set<int> removed;
    vs result;
    while (true) {
        int min = get_min(degrees);
        if (min == 1e9)
            break;
        removed.insert(min);
        if (degrees[min] != 0) {
            // NEEDS WORK
            int node_to_remove;//pick some node
            removed.insert(node_to_remove);
            result.push_back(create_equal(V[min], V[node_to_remove]));
        } else {
            result.push_back(V[min]);
        }
        degrees[min] = 1e9;
        for (int i = 0; i < V.size(); i++) {
            if (degrees[i] != 1e9) {
                int count = 0;
                for (int j = 0; j < adj[i].size(); j++) {
                    if (removed.find(j) == removed.end()) {
                        count++;
                    }
                }
                degrees[i] = count;
            }
        }
    }
    V = result;
}

int get_min(std::vector<int>& degrees) {
    int min = degrees[0];
    int ind = 0;
    for (int i = 1; i < degrees.size(); i++) {
        if (degrees[i] < min) {
            min = degrees[i];
            ind = i;
        }
    }
    if (min == 1e9)
        return -1;
    return ind;
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

/*
Gets the bit at position i of n.
*/
inline bool bit(int n, int i){
    return (n >> i) & 1;
}
