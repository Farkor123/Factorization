#pragma once
#include "lookupTable.h"
#include "utility.h"

namespace number_theorem {
  std::vector<mpz_class> pollard_rho (mpz_class& num) {
    /*
     *  Pollard's Monte Carlo as stated in his paper
     *  "A Monte Carlo method for factorization".
     *  With few addidtions.
     */
    //vector holding all factors
    std::vector<mpz_class> ret;
    //  Checking %2 is cheap, so we do that.
    while (num % 2 == 0) {
      ret.push_back(2);
      num /= 2;
    }
    //  If num was power of two, we can return.
    if (num == 1) {
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    //  If now num is lesser than 2**16, rho is overkill.
    //  We use lookup table instead.
    if (num < 65537) {
      auto w = lookup_table(num);
      ret.insert(ret.end(), w.begin(), w.end());
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    //  Rho doesn't like primes. We have to check it now.
    if (utility::is_prime(num)) {
      ret.push_back(num);
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    //  Innitialize random number generator.
    gmp_randclass r(gmp_randinit_default);
    //  Let's parametrize algorithm.
    mpz_class x = 2, y = 2, d = 1;
    //  I know, goto.
    //  But it's the simplest method, and doesn't crash.
    FLAG:
    //  The hearth of Pollard's Rho.
    while (d == 1) {
      x = utility::x_squared_minus_one(x) % num;
      y = utility::x_squared_minus_one(utility::x_squared_minus_one(y)) % num;
      d = gcd(abs(x - y), num);
    }
    //  If d==num, we can as well start from beginning,
    //  with different parameters.
    if (d == num) {
      x = y = r.get_z_range(num-3)+3;
      d = 1;
      goto FLAG;
    }
    num /= d;
    auto w = pollard_rho(d);
    ret.insert(ret.end(), w.begin(), w.end());
    w.clear();
    w = pollard_rho(num);
    ret.insert(ret.end(), w.begin(), w.end());
    std::sort(ret.begin(), ret.end());
    return ret;
  }
}
