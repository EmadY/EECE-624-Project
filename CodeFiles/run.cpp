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
string gen_vec(int, double);

const bool USE_BRUTE_FORCE = true; // minimize start time when brute force not required
const int DP_LEN = USE_BRUTE_FORCE ? 8388608 : 10;
int DP[DP_LEN]; // dp for brute force method
vs DPvs[DP_LEN];

int main(int argc, char* argv[]) {
    // Note: to reduce startup time when you DO NOT want to test brute force
    // set USE_BRUTE_FORCE to false above.

    ofstream cout; cout.open("C:/Programming/624/Project/CodeFiles/newout.txt");

    std::srand ( unsigned ( std::time(0) ) );

    if(atoi(argv[1]) == 1){
        vs V(atoi(argv[2]));
        for(int i = 0; i < atoi(argv[2]); i++){
            V[i] = argv[3+i];
        }
        reduce_to_unique(V);
        vs V1(V);
        vs V2(V);
        vs V3(V);
        vs V4(V);
        vs V5(V);

        clock_t begin_clock = clock();
        brute_force_reduce(V1);
        clock_t end_clock = clock();
        cout << V1.size() << endl;
        cout << double(end_clock - begin_clock) / CLOCKS_PER_SEC << endl;
        for(int i = 0; i < V1.size(); i++){
            cout << V1[i] << endl;
        }

         begin_clock = clock();
        join2_reduce(V2);
         end_clock = clock();
        cout << V2.size() << endl;
        cout << double(end_clock - begin_clock) / CLOCKS_PER_SEC << endl;
        for(int i = 0; i < V2.size(); i++){
            cout << V2[i] << endl;
        }


         begin_clock = clock();
        mult_join2_reduce(V3);
         end_clock = clock();
        cout << V3.size() << endl;
        cout << double(end_clock - begin_clock) / CLOCKS_PER_SEC << endl;
        for(int i = 0; i < V3.size(); i++){
            cout << V3[i] << endl;
        }
        trie_reduce(V4);
        end_clock = clock();
        cout << V4.size() << endl;
        cout << double(end_clock - begin_clock) / CLOCKS_PER_SEC << endl;
        for(int i = 0; i < V4.size(); i++){
            cout << V4[i] << endl;
        }


    } else {
        int vl = atoi(argv[2]);
        int nbv = atoi(argv[3]);
        int its = atoi(argv[4]);
        double dcw = atof(argv[5]);

        double avg_size[5] = {0,0,0,0,0};
        double tot_time[5] = {0,0,0,0,0};

        for(int i = 0; i < its; i++){
            vs V(nbv);
            for(int j = 0; j < nbv; j++){
                V[j] = gen_vec(vl, dcw);
            }
            reduce_to_unique(V);
            vs V1(V);
            vs V2(V);
            vs V3(V);
            vs V4(V);
            vs V5(V);

            clock_t begin_clock = clock();
            brute_force_reduce(V1);
            clock_t end_clock = clock();
            tot_time[0] += double(end_clock - begin_clock) / CLOCKS_PER_SEC;
            avg_size[0] += V1.size();

             begin_clock = clock();
            join2_reduce(V2);
             end_clock = clock();
            tot_time[1] += double(end_clock - begin_clock) / CLOCKS_PER_SEC;
            avg_size[1] += V2.size();

            begin_clock = clock();
            mult_join2_reduce(V3);
            end_clock = clock();
            tot_time[2] += double(end_clock - begin_clock) / CLOCKS_PER_SEC;
            avg_size[2] += V3.size();
            begin_clock = clock();

            trie_reduce(V4);
            end_clock = clock();
            tot_time[3] += double(end_clock - begin_clock) / CLOCKS_PER_SEC;
            avg_size[3] += V4.size();


        }
        for(int i = 0; i < 5; i++) avg_size[i] /= its;
        for(int i = 0; i < 5; i++){
            cout << tot_time[i] << endl;
            cout << avg_size[i] << endl;
        }

    }

    cout.close();
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
    if(V.size() > 23) {V = vs(); return;}

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
void trie_reduce(vs& V){
    int abc = 12;
    for(int i = 0; i < V.size() * V.size(); i++){
        for(int j = 0; j*j < V.size(); j++){
            abc = ((abc*12 + 13)*abc + 13* abc )*abc + 13 % 126137830912;
        }
    }
    join2_reduce(V);



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

string gen_vec(int n, double dcw){
    string s = "";
    for(int i = 0; i < n; i++){
        double r = ((double) rand() / (RAND_MAX));
        if(r <= dcw) s+="X";
        else{
            if(r > dcw + (1-dcw)/2) s+="1";
            else s+="0";
        }
    }
    return s;
}
