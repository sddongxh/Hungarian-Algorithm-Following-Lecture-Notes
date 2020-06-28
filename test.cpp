#include <ctime>
#include <iostream>
#include <vector>

#include "hungary_algorithm.hpp"

using namespace std;

vector<vector<double>> reverse(const vector<vector<double>>& A) {
  auto B = A;
  for (auto& p : B)
    for (auto& q : p) q = -q;
  return B;
}

int main(int argc, char** argv) {
  // Input matrix must be square
  vector<vector<double>> A = {{1, 6, 0}, {0, 8, 6}, {4, 0, 1}};
  std::srand(std::time(nullptr));  // use current time as seed for random generator
  int k = 10;
  A = vector<vector<double>>(k, vector<double>(k, 0));
  for (size_t i = 0; i < A.size(); i++) {
    for (size_t j = 0; j < A[0].size(); j++) {
      A[i][j] = (std::rand() % 100) / 100.0;
    }
  }
  auto [a, c] = HungaryAlgorithm::solve(A);
  cout << "reward = " << c << endl;
  cout << "assignment: ";
  for (auto p : a) cout << p << " ";
  cout << endl;
  return 0;
}