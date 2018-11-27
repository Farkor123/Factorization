#include "../lib/number_theorem.h"
#include <unistd.h>
int main() {
  /*system("clear");
  mpz_class num, comparator;
  std::cout << "Input number: ";
  std::cin >> num;
  if (number_theorem::is_prime(num)) {
    std::cout << "Number is prime!\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator.get_mpz_t(), 2, 16);
  if (num < comparator) {
    std::cout << "Factorization with lookup table:\n";
    auto x = number_theorem::lookup_table(num);
    for (unsigned int i = 0; i < x.size()-1; i++) {
      std::cout << x[i] << "*";
    }
    std::cout << x.back() << "\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator.get_mpz_t(), 2, 70);
  if (num < comparator) {
    std::cout << "Brent's modification of Pollard's rho algorithm:\n";
    auto x = number_theorem::pollard_rho(num);
    for (unsigned int i = 0; i < x.size()-1; i++) {
      std::cout << x[i] << "*";
    }
    std::cout << x.back() << "\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator.get_mpz_t(), 10, 50);
  if (num < comparator) {
    std::cout << "Lenstra's elliptic curve factorization.\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator.get_mpz_t(), 10, 100);
  if (num < comparator) {
    std::cout << "Quadratic Sieve.\n";
    return 0;
  }
  std::cout << "General Number Field Sieve.\n";*/
  int j;
  std::clock_t start;
  mpz_class i, y;
  for (j = 0; j < 5; j++) {
    start = std::clock();
    for (i = 40000; i < 60001; i++) {
      y = i;
      auto w = number_theorem::pollard_rho(y);
    }
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";
  }
  return 0;
}
