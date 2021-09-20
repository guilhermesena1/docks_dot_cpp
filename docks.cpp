// DOCKS algorithm for binary alphabets
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <bitset>

using std::bitset;
using std::runtime_error;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::max;

typedef uint32_t node_t;

struct DeBruijnGraph {
  explicit DeBruijnGraph(const size_t _k, const size_t _w);
  double estimate_ram() const;
  void exclude_decycling(const string &infile);
  void get_nodes_from(const node_t mer, node_t &from_1, node_t &from_2);
  void get_nodes_to(const node_t mer, node_t &to_1, node_t &to_2);
  void print_binary(node_t mer) const;

  bool count_paths_exclude_top();

  node_t sz; // the number of nodes in the graph
  size_t hash_mask; // the mask of all 1s in kmers

  vector<bool> V; // whether nodes are still in the graph
  vector<size_t> D;
  vector<size_t> F;
  vector<size_t> T;

  vector<node_t> from_1, from_2, to_1, to_2;
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

  cerr << "[required RAM: " << estimate_ram() << " MB]\n";
  V = vector<bool>(sz, true);
  
  const size_t dp_size = sz*(w+1);
  D = vector<size_t>(dp_size, 0);
  F = vector<size_t>(dp_size, 0);
  T = vector<size_t>(sz, 0);

  cerr << "[building edges]\n";
  from_1 = vector<node_t>(sz);
  from_2 = vector<node_t>(sz);
  to_1 = vector<node_t>(sz);
  to_2 = vector<node_t>(sz);
  for (node_t i = 0; i < sz; ++i) {
    get_nodes_from(i, from_1[i], from_2[i]);
    get_nodes_to(i, to_1[i], to_2[i]);
  }
};

double
DeBruijnGraph::estimate_ram() const {
  const size_t num_bytes = 
    sizeof(bool)*sz +
    sizeof(size_t)*(3*sz*(w+1)) +
    sizeof(node_t)*4;

  return static_cast<double>(num_bytes)/(1024.0*1024.0);
}

void
DeBruijnGraph::get_nodes_from(const node_t mer, node_t &f1, node_t &f2) {
  f1 = (mer >> 1) & hash_mask;
  f2 = ((mer >> 1) | (1ull << (DeBruijnGraph::k - 1))) & hash_mask;
}

void
DeBruijnGraph::get_nodes_to(const node_t mer, node_t &t1, node_t &t2) {
  t1 = (mer << 1) & hash_mask;
  t2 = ((mer << 1) | (1)) & hash_mask;
}

bool
DeBruijnGraph::count_paths_exclude_top() {
  const size_t lim = w + 1;
  const size_t start = 0;

#pragma omp parallel for
  for (node_t i = 0; i < sz; ++i)
    D[i*lim + start] = F[i*lim + start] = V[i];

  for (size_t j = start + 1; j < lim; ++j) {
#pragma omp parallel for
    for (node_t i = 0; i < sz; ++i) {
      if (V[i]) {
        // D: paths starting at i with length j
        D[i*lim + j] = V[to_1[i]]*D[to_1[i]*lim + j - 1] +
                       V[to_2[i]]*D[to_2[i]*lim + j - 1];

        // F: paths ending at i with length j
        F[i*lim + j] = V[from_1[i]]*F[from_1[i]*lim + j - 1] +
                       V[from_2[i]]*F[from_2[i]*lim + j - 1];


      }
    }
  }

  // count paths through each node
#pragma omp parallel for
  for (node_t i = 0; i < sz; ++i) {
    T[i] = 0;
    if (V[i]) {
      for (size_t j = 0; j < lim; ++j) {
        T[i] += D[i*lim + j]*F[i*lim + lim - j - 1];
      }
    }
  }

  // get which node to report
  const node_t the_best = max_element(begin(T), end(T)) - begin(T);

  // no more paths of length w
  if (!V[the_best] || !T[the_best]) return false;

  // report the k-mer with most paths
  print_binary(the_best);

  // remove it from the graph
  V[the_best] = false;
  return true;
}

static node_t
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
DeBruijnGraph::print_binary(const node_t mer) const {
  for (size_t i = 0; i < k; ++i)
    cout.put((mer & (1ull << (k - i - 1))) ? '1' : '0');
  cout.put('\n');
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

  cerr << "[done]\n";
}

int main(int argc, const char **argv) {
  try {
    if (argc != 4 && argc != 5) {
      cerr << "Usage: docks <k> <w> <decycling-k.txt> [num-threads]" << endl;
      return EXIT_FAILURE;
    }

    const size_t k = atoi(argv[1]);
    const size_t w = atoi(argv[2]);
    const string decycling_file = string(argv[3]);
    if (argc == 5)
      omp_set_num_threads(atoi(argv[4]));

    cerr << "[running docks with k = " << k << " and w = " << w << "]\n";
    cerr << "[using " << omp_get_num_threads() << " threads]\n";

    run_docks(k, w, decycling_file);
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  catch (std::bad_alloc &a) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
