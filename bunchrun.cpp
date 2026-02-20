#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>

using namespace std;

// Maximum possible N for m=12
const int MAX_N = (1 << 12);
const int MAX_WORDS = (MAX_N + 64) / 64;

// Global storage to avoid re-allocation
uint64_t adj_bitset[MAX_N][MAX_WORDS];
int degrees[MAX_N];

// ---------------------------------------------------------
// VERSION 1: STANDARD TOMITA (Dynamic Allocation)
// ---------------------------------------------------------
struct Tomita {
    uint64_t nodes, cliques;
    int words, n_val;

    int choosePivot(const vector<uint64_t>& P, const vector<uint64_t>& X) {
        int best_u = -1, max_intersect = -1;
        for (int i = 0; i < words; i++) {
            uint64_t candidates = P[i] | X[i];
            while (candidates) {
                int u = (i << 6) + __builtin_ctzll(candidates);
                int intersect = 0;
                for (int k = 0; k < words; k++) intersect += __builtin_popcountll(adj_bitset[u][k] & P[k]);
                if (intersect > max_intersect) { max_intersect = intersect; best_u = u; }
                candidates &= candidates - 1;
            }
        }
        return best_u;
    }

    void expand(vector<uint64_t> P, vector<uint64_t> X) {
        nodes++;
        bool p_empty = true;
        for (int i = 0; i < words; i++) if (P[i]) { p_empty = false; break; }
        if (p_empty) {
            for (int i = 0; i < words; i++) if (X[i]) return;
            cliques++; return;
        }

        int u = choosePivot(P, X);
        for (int i = 0; i < words; i++) {
            uint64_t mask = P[i] & ~adj_bitset[u][i];
            while (mask) {
                int v = (i << 6) + __builtin_ctzll(mask);
                vector<uint64_t> P_next(words), X_next(words);
                for (int k = 0; k < words; k++) {
                    P_next[k] = P[k] & adj_bitset[v][k];
                    X_next[k] = X[k] & adj_bitset[v][k];
                }
                expand(P_next, X_next);
                P[i] &= ~(1ULL << (v & 63));
                X[i] |= (1ULL << (v & 63));
                mask &= mask - 1;
            }
        }
    }
};

// ---------------------------------------------------------
// VERSION 2: HYBRID DLX (Stack-based + Short-circuit)
// ---------------------------------------------------------
struct HybridDLX {
    uint64_t nodes, cliques;
    int words, n_val;
    uint64_t p_stack[200][MAX_WORDS];
    uint64_t x_stack[200][MAX_WORDS];

    int choosePivot(const uint64_t* P, const uint64_t* X) {
        int best_u = -1, max_p = -1, p_size = 0;
        for (int i = 0; i < words; i++) p_size += __builtin_popcountll(P[i]);
        if (p_size <= 2) {
            for (int i = 0; i < words; i++) if (P[i]) return (i << 6) + __builtin_ctzll(P[i]);
        }
        const int threshold = (p_size * 9) / 10;

        auto check = [&](const uint64_t* set) {
            for (int k = 0; k < words; k++) {
                uint64_t bits = set[k];
                while (bits) {
                    int u = (k << 6) + __builtin_ctzll(bits);
                    int cnt = 0;
                    for (int i = 0; i < words; i++) cnt += __builtin_popcountll(adj_bitset[u][i] & P[i]);
                    if (cnt > max_p) { max_p = cnt; best_u = u; }
                    if (max_p >= threshold) return true;
                    bits &= bits - 1;
                }
            }
            return false;
        };
        if (!check(P)) check(X);
        return best_u;
    }

    void expand(int depth, uint64_t* P, uint64_t* X) {
        nodes++;
        bool p_empty = true;
        for (int i = 0; i < words; i++) if (P[i]) { p_empty = false; break; }
        if (p_empty) {
            for (int i = 0; i < words; i++) if (X[i]) return;
            cliques++; return;
        }

        int u = choosePivot(P, X);
        for (int i = 0; i < words; i++) {
            uint64_t mask = P[i] & ~adj_bitset[u][i];
            while (mask) {
                int v = (i << 6) + __builtin_ctzll(mask);
                for (int k = 0; k < words; k++) {
                    p_stack[depth+1][k] = P[k] & adj_bitset[v][k];
                    x_stack[depth+1][k] = X[k] & adj_bitset[v][k];
                }
                expand(depth + 1, p_stack[depth+1], x_stack[depth+1]);
                P[i] &= ~(1ULL << (v & 63));
                X[i] |= (1ULL << (v & 63));
                mask &= mask - 1;
            }
        }
    }
};

void run_experiment(int m) {
    int n = (1 << m) - 1;
    int words = (n + 64) / 64;

    // Reset Globals
    for (int i = 0; i <= n; i++) {
        degrees[i] = 0;
        for (int j = 0; j < words; j++) adj_bitset[i][j] = 0;
    }

    // Build Graph
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i != j && (i & j) == 0) {
                adj_bitset[i][j >> 6] |= (1ULL << (j & 63));
                degrees[i]++;
            }
        }
    }

    // Prepare P
    vector<uint64_t> p_init(words, 0);
    for (int i = 1; i <= n; i++) p_init[i >> 6] |= (1ULL << (i & 63));

    // Run Tomita
    Tomita t = {0, 0, words, n};
    auto s1 = chrono::high_resolution_clock::now();
    t.expand(p_init, vector<uint64_t>(words, 0));
    auto e1 = chrono::high_resolution_clock::now();
    double time_t = chrono::duration<double>(e1 - s1).count();

    // Run Hybrid DLX
    HybridDLX h = {0, 0, words, n};
    uint64_t p_hybrid[MAX_WORDS] = {0}, x_hybrid[MAX_WORDS] = {0};
    for(int i=0; i<words; i++) p_hybrid[i] = p_init[i];
    auto s2 = chrono::high_resolution_clock::now();
    h.expand(0, p_hybrid, x_hybrid);
    auto e2 = chrono::high_resolution_clock::now();
    double time_h = chrono::duration<double>(e2 - s2).count();

    // Print Row
    cout << "|" << setw(3) << m << " |" << setw(10) << h.cliques << " |" << setw(10) << h.nodes 
         << " |" << setw(12) << fixed << setprecision(5) << time_t 
         << " |" << setw(12) << time_h << " |" << setw(8) << setprecision(2) << time_t/time_h << "x |" << endl;
}

int main() {
    cout << "| m  | Cliques    | Nodes      | Tomita (s)   | Hybrid (s)   | Speedup  |" << endl;
    cout << "|----|------------|------------|--------------|--------------|----------|" << endl;
    for (int m_val = 4; m_val <= 12; m_val++) {
        run_experiment(m_val);
    }
    return 0;
}