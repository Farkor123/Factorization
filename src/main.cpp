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
  gmp_randclass r(gmp_randinit_default);
  while (true) {
    mpz_class num;
    mpz_ui_pow_ui(num.get_mpz_t(), 2, 70);
    mpz_class y, z, c = 1;
    y = z = r.get_z_range(num);
    std::cout << y << "\n";
    auto w = number_theorem::simpler_pollard_rho(y);
    for (auto i : w) {
      c *= i;
    }
    if (c == z) {
      std::cout << " Yes: ";
      for (unsigned int i = 0; i < w.size()-1; i++) {
        std::cout << w[i] << "*";
      }
      std::cout << w.back() << "\n";
    }
    else {
      std::cout << " No\n";
      return -1;
    }
  }
  return 0;
}
