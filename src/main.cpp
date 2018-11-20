#include "../lib/number_theorem.h"

int main() {
  system("clear");
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
    for (int i = 0; i < number_theorem::lookup_table(num).size()-1; i++) {
      std::cout << number_theorem::lookup_table(num)[i] << "*";
    }
    std::cout << number_theorem::lookup_table(num)[number_theorem::lookup_table(num).size()-1] << "\n";
    return 0;
  }
  mpz_ui_pow_ui(comparator.get_mpz_t(), 2, 70);
  if (num < comparator) {
    std::cout << "Brent's modification of Pollard's rho algorithm:\n";
    //number_theorem::pollard_rho(num);
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
  std::cout << "General Number Field Sieve.\n";
  return 0;
}
