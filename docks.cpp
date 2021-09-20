// DOCKS algorithm for binary alphabets
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

using std::runtime_error;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::max;

struct DeBruijnGraph {
  explicit DeBruijnGraph(const size_t _k, const size_t _w);
  double estimate_ram();
  void exclude_decycling(const string &infile);
  size_t count_paths_exclude_top();
  void print_excluded_nodes();
  void get_nodes_from(const size_t mer, size_t &from_1, size_t &from_2);
  void get_nodes_to(const size_t mer, size_t &to_1, size_t &to_2);
  void print_binary(size_t mer);


  size_t sz;
  size_t hash_mask;
  size_t top_node;

  vector<bool> V;
  vector< vector<size_t> > D;
  vector< vector<size_t> > F;
  vector< vector<size_t> > T;
  static size_t k;
  static size_t w;
};

size_t DeBruijnGraph::k = 0;
size_t DeBruijnGraph::w = 0;

DeBruijnGraph::DeBruijnGraph(const size_t _k, const size_t _w) {
  DeBruijnGraph::w = _w;
  DeBruijnGraph::k = _k;
  sz = (1 << k);
  hash_mask = (1ull << k) - 1;

  cerr << "[required RAM: " <<  estimate_ram() << "]\n";
  V = vector<bool>(sz, true);
  D = vector<vector<size_t>>(sz, vector<size_t>(w+1, 0));
  F = vector<vector<size_t>>(sz, vector<size_t>(w+1, 0));
  T = vector<vector<size_t>>(sz, vector<size_t>(w+1, 0));
};

double
DeBruijnGraph::estimate_ram() {
  const size_t num_bytes = 
    sizeof(bool)*sz + //V
    sizeof(size_t)*(3*sz*(w+1));

  return static_cast<double>(num_bytes)/(1024.0*1024.0);
}

void
DeBruijnGraph::get_nodes_from(const size_t mer, size_t &from_1, size_t &from_2) {
  from_1 = (mer >> 1) & hash_mask;
  from_2 = ((mer >> 1) | (1ull << DeBruijnGraph::k)) & hash_mask;
}

void
DeBruijnGraph::get_nodes_to(const size_t mer, size_t &to_1, size_t &to_2) {
  to_1 = (mer << 1) & hash_mask;
  to_2 = ((mer << 1) | (1)) & hash_mask;
}

size_t
DeBruijnGraph::count_paths_exclude_top() {
  // get edges from binary hash
  size_t from_1, from_2, to_1, to_2;

  size_t ans = 0;
  size_t top_node = 0;
  for (size_t j = 1; j <= w; ++j) {
    for (size_t i = 0; i < sz; ++i) {
      if (j == 1) {
        D[i][1] = F[i][1] = T[i][1] = V[i];
      }
      else if (V[i]) {
        get_nodes_from(i, from_1, from_2);
        get_nodes_from(i, to_1, to_2);
        D[i][j] = V[from_1]*D[from_1][j-1] + V[from_2]*D[from_2][j-1];
        F[i][j] = V[to_1]*F[to_1][j-1] + V[to_2]*F[to_2][j-1];
        T[i][j] = D[i][j]*F[i][j];
      }
    }
  }

  // find node with most paths
  for (size_t i = 0; i < sz; ++i)
    ans = (V[i]) ? (max(ans, T[i][w])) : ans;

  // find the top node
  for (size_t i = 0; i < sz; ++i)
    top_node = ((V[i]) && (T[i][w] == ans)) ? (i) : (top_node);

  // exclude node with most paths
  print_binary(top_node);
  V[top_node] = false;
  return ans;
}

static size_t
to_binary(const string &s) {
  size_t ans = 0ull;
  const size_t lim = s.size();

  for (size_t i = 0; i < lim; ++i)
    ans = (ans << 1) | (s[i] == '1');

  return ans;
}

void
DeBruijnGraph::exclude_decycling(const string &infile) {
  ifstream in(infile);
  if (!in)
    throw runtime_error("bad decycling file: " + infile);

  string mer;
  while (in >> mer)
    V[to_binary(mer)] = false;

}

void
DeBruijnGraph::print_binary(size_t mer) {
  for (size_t i = 0; i < k; ++i) {
    cout.put((mer & 1) ? '1' : '0');
    mer >>= 1;
  }
  cout.put('\n');
}

void
DeBruijnGraph::print_excluded_nodes() {
  for (size_t i = 0; i < sz; ++i)
    if (!V[i])
      print_binary(i);
}

void
run_docks(const size_t key_weight, const size_t window_size,
          const string &decycling_file) {
  cerr << "[constructing graph]\n";
  DeBruijnGraph G(key_weight, window_size);

  cerr << "[excluding decycling]\n";
  G.exclude_decycling(decycling_file);

  cerr << "[main loop]\n";
  while (G.count_paths_exclude_top() > 0);

  //cerr << "[printing nodes to cout]\n";
  //G.print_excluded_nodes();

  cerr << "[done]\n";
}

int main(int argc, const char **argv) {
  if (argc != 4) {
    cerr << "Usage: docks <k> <w> <decycling-k.txt>" << endl;
    return EXIT_FAILURE;
  }

  const size_t k = atoi(argv[1]);
  const size_t w = atoi(argv[2]);
  const string decycling_file = string(argv[3]);

  cerr << "[running docks with k = " << k << " and w = " << w << "]\n";

  run_docks(k, w, decycling_file);
  return EXIT_SUCCESS;
}
