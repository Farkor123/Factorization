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
    //  Innitialize random number generator.
    gmp_randclass r(gmp_randinit_default);
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
    else {
      //  Else, if d is prime it is a valid factor.
      if (utility::is_prime(d)) {
        num /= d;
        ret.push_back(d);
        //  If num is now prime, there is no point in further fctorization.
        if(utility::is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        //  If now num is less than 2**16, rho is overkill.
        //  We will factor it with lookup table.
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
      //  If d ain't prime, but is less than 65537, rho is overkill.
      //  We will use this complex d, factoring it with lookup table.
      else if(d < 65537) {
        num /= d;
        auto w = lookup_table(d);
        ret.insert(ret.end(), w.begin(), w.end());
        //  If num is now prime, there is no point in further fctorization.
        if (utility::is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        //  If now num is less than 2**16, rho is overkill.
        //  We will factor it with lookup table.
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
      //  However, if d > 2**16, we will just use Rho on it.
      else {
        num /= d;
        auto w = pollard_rho(d);
        ret.insert(ret.end(), w.begin(), w.end());
        //  If num is now prime, there is no point in further fctorization.
        if (utility::is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        //  If now num is less than 2**16, rho is overkill.
        //  We will factor it with lookup table.
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
    }
  }
}
