#pragma once
#include "lookupTable.h"
#include "utility.h"

namespace number_theorem {
  std::vector<mpz_class> brent_pollard_rho (mpz_class& num) {
  /*
   * An implementation from Brent's paper "An improved Monte Carlo factorization algorithm."
   */
    // A vector to hold all factors.
    std::vector<mpz_class> ret;
    while (num % 2 == 0) {
      ret.push_back(2);
      num /= 2;
    }
    if (num == 1) {
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    if (utility::is_prime(num)) {
      ret.push_back(num);
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    if (num < 65537) {
      auto w = lookup_table(num);
      ret.insert(ret.end(), w.begin(), w.end());
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    //  Innitialize random number generator.
    gmp_randclass rng(gmp_randinit_default);
    mpz_class y, m, g=1, r=1, q=1, ys, x;
    y = rng.get_z_range(num-1);
    m = rng.get_z_range(num-1);
    while (g == 1) {
      x = y;
      for (mpz_class i = 0; i < r; i++) {
        y = utility::x_squared_minus_one(y) % num;
      }
      mpz_class k = 0;
      while (k < r && g == 1) {
        ys = y;
        mpz_class min;
        m <= r-k ? min = m : min = r - k;
        for (mpz_class i = 0; i < min; i++) {
          y = utility::x_squared_minus_one(y) % num;
          q = (q * (abs(x - y))) % num;
        }
        k += m;
      }
      r *= 2;
    }
    //  If d==num, we can as well start from beginning,
    //  with different parameters.
    if (g == num) {
      while (true) {
        ys = utility::x_squared_minus_one(ys) % num;
        g = gcd(abs(x - ys), num);
        break;
      }
    }
    num /= g;
    auto w = brent_pollard_rho(g);
    ret.insert(ret.end(), w.begin(), w.end());
    w.clear();
    w = brent_pollard_rho(num);
    ret.insert(ret.end(), w.begin(), w.end());
    std::sort(ret.begin(), ret.end());
    return ret;
  }
}
