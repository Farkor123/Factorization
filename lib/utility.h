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
}
