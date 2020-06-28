#pragma once
#include <algorithm>
#include <cmath> 
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class HungaryAlgorithm {
 public:
  static std::pair<std::vector<int>, double> solve(const std::vector<std::vector<double>>& A) {
    // step 1
    auto l = HungaryAlgorithm::trivial_labeling(A);  // labelling
    auto M = HungaryAlgorithm::match(A, l);
    while (true) {
      // step 2
      if (M.size() == 2 * A.size()) break;
      auto u = HungaryAlgorithm::free_vertices(M, A.size()).front();
      std::unordered_set<int> S = {u};
      std::unordered_set<int> T = {};
      std::unordered_map<int, std::unordered_set<int>> alternating_tree;
      alternating_tree[u] = {};               // root
      std::unordered_map<int, double> slack;  // for each column
      for (size_t j = 0; j < A[0].size(); j++) slack[A.size() + j] = l[u] + l[A.size() + j] - A[u][j];

      // step 3
      std::unordered_set<int> N;
      while (true) {
        if (N.size() == T.size()) {
          auto it = HungaryAlgorithm::min_slack(slack);
          double alpha = it.second;
          N.insert(it.first);
          slack.erase(it.first);
          for (auto& p : slack) p.second -= alpha;
          for (size_t i = 0; i < l.size(); i++) {
            if (S.count(i)) l[i] -= alpha;
            if (T.count(i)) l[i] += alpha;
          }
        }
        // step 4
        auto y = HungaryAlgorithm::set_difference(N, T).front();
        auto s = HungaryAlgorithm::neighbor_reverse(A, S, l, y).front();
        alternating_tree[s].insert(y);
        if (M.count(y) == 0) {  // y is free
          auto path = HungaryAlgorithm::find_path(alternating_tree, u, y);
          for (size_t i = 0; i + 1 < path.size(); i++) {
            if (i % 2 == 0) {
              M[path[i]] = path[i + 1];
              M[path[i + 1]] = path[i];
            }
          }
          break;
        } else {  // y is matched
          int z = M[y];
          S.insert(z);
          T.insert(y);
          N.insert(y);
          alternating_tree[y].insert(z);
          // update slack
          for (auto& it : slack) it.second = std::min(it.second, l[z] + l[it.first] - A[z][it.first - A.size()]);
        }
      }
      if (M.size() == 2 * A.size()) break;
    }
    std::vector<int> res(A.size(), 0);
    double c = 0;
    for (size_t i = 0; i < res.size(); i++) {
      res[i] = M[i] - A.size();
      c += A[i][res[i]];
    }
    return {res, c};
  }

 private:
  static std::vector<double> trivial_labeling(const std::vector<std::vector<double>>& A) {
    std::vector<double> label(A.size() + A[0].size(), 0);
    for (size_t i = 0; i < A.size(); i++) label[i] = *max_element(A[i].begin(), A[i].end());
    return label;
  }

  static std::unordered_map<int, int> match(const std::vector<std::vector<double>>& A, const std::vector<double>& l) {
    std::unordered_map<int, int> M;
    for (size_t i = 0; i < A.size(); i++) {
      for (size_t j = 0; j < A[0].size(); j++) {
        if (M.count(A.size() + j) > 0 || M.count(i) > 0) continue;
        if (A[i][j] == l[i] + l[A.size() + j]) {
          M[i] = A.size() + j;
          M[A.size() + j] = i;
        }
      }
    }
    return M;
  }

  static std::unordered_set<int> neighbor(const std::vector<std::vector<double>>& A, const std::vector<double>& l,
                                          int i) {
    std::unordered_set<int> res;
    for (size_t j = 0; j < A[0].size(); j++)
      if (A[i][j] == l[i] + l[A.size() + j]) res.insert(A.size() + j);
    return res;
  }

  static std::unordered_set<int> neighbor(const std::vector<std::vector<double>>& A, const std::vector<double>& l,
                                          std::unordered_set<int>& S) {
    std::unordered_set<int> res;
    for (auto i : S) {
      auto p = neighbor(A, l, i);
      res.insert(p.begin(), p.end());
    }
    return res;
  }

  static std::vector<int> neighbor_reverse(const std::vector<std::vector<double>>& A, const std::unordered_set<int>& S,
                                           const std::vector<double>& l, int j) {
    std::vector<int> res;
    std::pair<int, double> candidate = {-1, std::numeric_limits<double>::max()};
    for (auto i : S) {
      // if (A[i][j - A.size()] == l[i] + l[j]) res.push_back(i); // ideal but not stable
      double d = fabs(A[i][j - A.size()] - l[i] - l[j]);
      if (d < candidate.second) {
        candidate.second = d;
        candidate.first = i;
      }
    }
    res.push_back(candidate.first);
    return res;
  }

  static std::vector<int> set_difference(const std::unordered_set<int>& a, const std::unordered_set<int>& b) {
    std::vector<int> res;
    for (auto p : a)
      if (b.count(p) == 0) res.push_back(p);
    return res;
  }

  static std::vector<int> free_vertices(const std::unordered_map<int, int>& M, int m) {
    std::vector<int> res;
    for (int i = 0; i < m; i++)
      if (M.count(i) == 0) res.push_back(i);
    return res;
  }

  static std::vector<int> find_path(const std::unordered_map<int, std::unordered_set<int>>& tree, int s, int d) {
    if (s == d) return {d};
    if (tree.count(s) != 0) {
      for (auto p : tree.at(s)) {
        auto t = find_path(tree, p, d);
        if (!t.empty()) {
          t.insert(t.begin(), s);
          return t;
        }
      }
    }
    return {};
  }

  static std::pair<int, double> min_slack(const std::unordered_map<int, double>& slack) {
    std::pair<int, double> candidate = {-1, std::numeric_limits<double>::max()};
    for (auto& p : slack) {
      if (p.second < candidate.second) {
        candidate.first = p.first;
        candidate.second = p.second;
      }
    }
    return candidate;
  }
};
