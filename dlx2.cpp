#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;

const int m = 12;
const int n = (1 << m) - 1;
const size_t WORDS = (n + 64) / 64;

// Global state
uint64_t adj[n + 1][WORDS] = {{0}};
int L[n + 2], R[n + 2];

void dlx_remove(int i) {
    R[L[i]] = R[i];
    L[R[i]] = L[i];
}

void dlx_restore(int i) {
    L[R[i]] = i;
    R[L[i]] = i;
}

class FastHybrid {
    uint64_t nodes = 0, cliques = 0;
    uint64_t p_stack[110][WORDS], x_stack[110][WORDS];

public:
    int get_pivot(const uint64_t* P, const uint64_t* X) {
        int best_u = -1;
        int max_intersect = -1;
        auto search = [&](const uint64_t* set) {
            for (size_t k = 0; k < WORDS; ++k) {
                uint64_t bits = set[k];
                while (bits) {
                    int u = (int)(k << 6) + __builtin_ctzll(bits);
                    if (u <= n) {
                        int count = 0;
                        for (size_t i = 0; i < WORDS; ++i)
                            count += __builtin_popcountll(adj[u][i] & P[i]);
                        if (count > max_intersect) {
                            max_intersect = count;
                            best_u = u;
                        }
                    }
                    bits &= (bits - 1);
                }
            }
        };
        search(P); search(X);
        return best_u;
    }

    void expand(int depth, uint64_t* P, uint64_t* X) {
        nodes++;
        bool p_empty = true;
        for (size_t i = 0; i < WORDS; i++) if (P[i]) { p_empty = false; break; }
        if (p_empty) {
            for (size_t i = 0; i < WORDS; i++) if (X[i]) return;
            cliques++; return;
        }

        int u = get_pivot(P, X);
        vector<int> cand;
        // Use DLX to find candidates
        for (int v = R[0]; v != 0; v = R[v]) {
            size_t v_idx = (size_t)v;
            if ((P[v_idx >> 6] & (1ULL << (v_idx & 63))) && !(adj[u][v_idx >> 6] & (1ULL << (v_idx & 63)))) {
                cand.push_back(v);
            }
        }

        for (int v : cand) {
            dlx_remove(v); // The "Dancing" part
            size_t v_idx = (size_t)v;
            for (size_t k = 0; k < WORDS; k++) {
                p_stack[depth + 1][k] = P[k] & adj[v][k];
                x_stack[depth + 1][k] = X[k] & adj[v][k];
            }
            expand(depth + 1, p_stack[depth + 1], x_stack[depth + 1]);
            
            P[v_idx >> 6] &= ~(1ULL << (v_idx & 63));
            X[v_idx >> 6] |= (1ULL << (v_idx & 63));
        }
        for (int v : cand) dlx_restore(v);
    }

    void run() {
        R[0] = 1; L[0] = n;
        for (int i = 1; i <= n; i++) {
            L[i] = i - 1; R[i] = (i == n ? 0 : i + 1);
        }
        uint64_t P[WORDS] = {0}, X[WORDS] = {0};
        for (int i = 1; i <= n; i++) P[(size_t)i >> 6] |= (1ULL << ((size_t)i & 63));
        expand(0, P, X);
        cout << "Maximal Cliques: " << cliques << "\nNodes: " << nodes << endl;
    }
};

int main() {
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i != j && (i & j) == 0) {
                adj[i][(size_t)j >> 6] |= (1ULL << ((size_t)j & 63));
            }
        }
    }
    FastHybrid solver;
    auto s = chrono::high_resolution_clock::now();
    solver.run();
    auto e = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration<double>(e - s).count() << "s" << endl;
    return 0;
}