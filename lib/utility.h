#pragma once
#include <bits/stdc++.h>
#include <gmp.h>
#include <gmpxx.h>

namespace utility {

  void mpz_init_set_cin(mpz_t &rop, int base) {
    std::string str;
    std::cin >> str;
    char arr[1024];
    std::strcpy(arr, str.c_str());
    mpz_init_set_str(rop, arr, base);
  }

  void mpz_gcd (mpz_t &rop, const mpz_t& _a, const mpz_t& _b) {
    mpz_t temp, a, b;
    mpz_init_set(a, _a);
    mpz_init_set(b, _b);
    mpz_init_set_ui(temp, 0);
    while (mpz_cmp_ui(b, 0) != 0) {
      mpz_mod(temp, a, b);
      mpz_set(a, b);
      mpz_set(b, temp);
    }
    mpz_set(rop, a);
    mpz_clear(temp);
    mpz_clear(a);
    mpz_clear(b);
  }

  void mpz_pow_mod(mpz_t& rop, const mpz_t& _base, const mpz_t& _exp, const mpz_t& mod) {
    mpz_t base, exp;
    mpz_init_set(base, _base);
    mpz_init_set(exp, _exp);
    mpz_set_ui(rop, 1);
    mpz_mod(base, base, mod);
    while (mpz_cmp_ui(exp, 0) > 0) {
      if(mpz_odd_p(exp) != 0) {
        mpz_mul(rop, rop, base);
        mpz_mod(rop, rop, mod);
      }
      mpz_fdiv_q_ui(exp, exp, 2);
      mpz_mul(base, base, base);
      mpz_mod(base, base, mod);
    }
    mpz_clear(base);
    mpz_clear(exp);
  }
}
