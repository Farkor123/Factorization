#pragma once
#include "utility.h"

namespace number_theorem {
  int least_divisor(const int& num) {
    for(int i = 3; i <= sqrt(num); i += 2) {
      if (num % i == 0) {
        return i;
      }
    }
    return num;
  }

  void lookup_table (unsigned long int& num) {
    int arr[65537] = {2};
    for (int i = 2; i < 65537; i++) {
      if (i % 2 == 0) {
        arr[i] = 2;
      }
      else {
        arr[i] = least_divisor(i);
      }
    }
    while (num != 1) {
      std::cout << arr[num];
      num /= arr[num];
      if(num != 1) {
        std::cout << "*";
      }
    }
    std::cout << "\n";
  }

  bool miller_rabin_test (const mpz_t& _d, const mpz_t& _n) {
    mpz_t a, d, n, temp;
    mpz_init(temp);
    mpz_init_set(d, _d);
    mpz_init_set(n, _n);
    mpz_init_set_ui(a, rand());
    mpz_sub_ui(temp, n, 4);
    mpz_mod(a, a, temp);
    mpz_add_ui(a, a, 2);
    mpz_sub_ui(temp, n, 1);
    utility::mpz_pow_mod(a, a, d, n);
    if (mpz_cmp_ui(a, 1) == 0 || mpz_cmp(a, temp) == 0) {
      return true;
    }
    while (mpz_cmp(d, temp) != 0) {
      mpz_mul(a, a, a);
      mpz_mod(a, a, n);
      mpz_mul_ui(d, d, 2);
      if (mpz_cmp_ui(a, 1) == 0) {
        mpz_clear(a);
        mpz_clear(d);
        mpz_clear(n);
        mpz_clear(temp);
        return false;
      }
      if (mpz_cmp(a, temp) == 0) {
        mpz_clear(a);
        mpz_clear(d);
        mpz_clear(n);
        mpz_clear(temp);
        return true;
      }
    }
    mpz_clear(a);
    mpz_clear(d);
    mpz_clear(n);
    mpz_clear(temp);
    return false;
  }

  bool is_prime(const mpz_t &num, int k = 256) {
    if (mpz_cmp_ui(num, 1) <= 0 || mpz_cmp_ui(num, 4) == 0) {
      return false;
    }
    if (mpz_cmp_ui(num, 3) == 0 || mpz_cmp_ui(num, 2) == 0) {
      return true;
    }
    mpz_t d;
    mpz_init(d);
    mpz_sub_ui(d, num, 1);
    while (mpz_even_p(d) != 0) {
      mpz_fdiv_q_ui(d, d, 2);
    }
    for (int i = 0; i < k; i++) {
      if(!miller_rabin_test(d, num)) {
        return false;
      }
    }
    return true;
  }

  void pollard_rho (mpz_t& num) {
    mpz_t x[100001];
    mpz_t q, d, diff;
    mpz_init_set_ui(q, 1);
    mpz_inits(d, diff);
    mpz_init_set_ui(x[1], 3);
    mpz_init_set_ui(x[2], 8);
    for (int i = 1; ; i++) {
      mpz_sub(diff, x[2*i], x[i]);
      mpz_mul(q, q, diff);
      mpz_gcd(d, q, num);
      if (mpz_cmp_ui(d, 1) > 0) {
        std::cout << d << "*";
      }
    }
  }
}
