#include "number_theorem.h"

int main() {
  system("clear");
  mpz_t num, comparator;
  std::cout << "Input number: ";
  utility::mpz_init_set_cin(num, 10);
  mpz_init(comparator);
  mpz_ui_pow_ui(comparator, 2, 16);
  if (!(mpz_cmp(num, comparator) > 0)) {
    mpz_clear(comparator);
    unsigned long int x = mpz_get_ui(num);
    mpz_clear(num);
    std::cout << "Factorization with lookup table:\n";
    number_theorem::lookup_table(x);
    return 0;
  }
  mpz_ui_pow_ui(comparator, 2, 70);
  if (!(mpz_cmp(num, comparator) > 0)) {
    mpz_clear(comparator);
    std::cout << "Brent's modification of Pollard's rho algorithm:\n";
    number_theorem::pollard_rho(num);
    return 0;
  }
  mpz_ui_pow_ui(comparator, 10, 50);
  if (!(mpz_cmp(num, comparator) > 0)) {
    mpz_clear(comparator);
    std::cout << "Lenstra's elliptic curve factorization.\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator, 10, 100);
  if (!(mpz_cmp(num, comparator) > 0)) {
    mpz_clear(comparator);
    std::cout << "Quadratic Sieve.\n";
    return 0;
  }
  mpz_clear(comparator);
  std::cout << "General Number Field Sieve.\n";
  return 0;
}
