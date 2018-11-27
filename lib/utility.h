#pragma once
#include <bits/stdc++.h>
#include <gmp.h>
#include <gmpxx.h>

namespace utility {
  mpz_class gcd (mpz_class a, mpz_class b) {
    mpz_class c;
    while (b != 0) {
      c = a % b;
      a = b;
      b = c;
    }
    return a;
  }

  mpz_class pow_mod(mpz_class base, mpz_class exponent, const mpz_class& mod) {
    mpz_class ret = 1;
    base %= mod;
    while (exponent > 0) {
      if (exponent % 2 != 0) {
        ret = (ret * base) % mod;
      }
      exponent = exponent>>1;
      base = (base * base) % mod;
    }
    return ret;
  }

  mpz_class x_squared_minus_one(const mpz_class& x) {
    return x*x-1;
  }

  bool miller_rabin_test (mpz_class d, const mpz_class& n) {
    mpz_class x = 2 + rand() % (n-4);
    x = utility::pow_mod(x, d, n);
    if (x == 1 || x == n-1) {
      return true;
    }
    while (d != n-1) {
      x = (x * x) % n;
      d *= 2;
      if (x == 1) {
        return false;
      }
      if (x == n-1) {
        return true;
      }
    }
    return false;
  }

  bool is_prime(const mpz_class &num, int k = 256) {
    if (num <= 1 || num == 4) {
      return false;
    }
    if (num == 2 || num == 3) {
      return true;
    }
    mpz_class d = num - 1;
    while (d % 2 == 0) {
      d /= 2;
    }
    for (int i = 0; i < k; i++) {
      if(!miller_rabin_test(d, num)) {
        return false;
      }
    }
    return true;
  }
}
